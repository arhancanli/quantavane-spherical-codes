"""The entries above the size limit whose cosine the table gives exactly (census version 8): Proposition 3.19 in exact
arithmetic, with a sketched rank test over F_p that fits in memory.

  nice -n 10 python3 -I run_above.py OUTDIR [--only d_n,...] [--max-cols 20000] [--mem-mb 1500]

For every entry of the archived table with 3 <= d <= 32, N <= 1024, no verdict because of size, and a minimal polynomial
of degree at most 2 (K = Q or Q(sqrt D), D the squarefree part of the discriminant):
 1. certify/v5/core.py (released, unchanged) builds the exact core: the core's Gram matrix recognised over K and checked
    exactly (unit diagonal, rank r, G = G_.B G_BB^-1 G_B.), the contacts E at the exact largest cosine c with every other
    core pair strictly below c, the stored rattlers proved below c (balls), and the certified distance to the stored
    code (an exactly orthogonal rational frame). If r < d the core is flat and Proposition 2.7 applies.
 2. If r = d: bar-rigidity by Proposition 3.19(c,d): the rank over F_p (p < 2^16 prime, sqrt D -> a root mod p) of
    S M, with M the rigidity matrix in the coordinates of the basis B (the rows of v5/rigidity.py R_rows_modp) and
    S = M^T, or, if that falls short, S uniformly random with |V| d rows and another prime (any S gives a lower bound
    for the rank over K). The rank is computed by blocked elimination in double precision on a matrix stored as 16-bit
    integers, exact because every intermediate value is an integer below 2^53 (entries below 2^16, products below 2^32,
    dot products of at most 256 terms). The infinitesimal rotations of the whole core and, if E is disconnected, of
    each component, are checked to be exact kernel vectors (over K), and independent mod p. Bar-rigid (modulo these)
    when the rank equals |V| d minus their number.
 3. A strictly positive equilibrium stress, checked exactly over K: the uniform one; else one constant on classes of
    contacts (by common neighbours; by the values of the stress closest to the uniform one, a projection computed by
    sparse least squares; or, up to 4,000 columns, by the values of an interior-point linear-programming stress; the
    last two are constant on the orbits of the symmetries): an exact basis of the stresses constant on the classes
    (the nullspace over Q of the equilibrium equations, by FLINT; over Q(sqrt D) as a doubled system over Q), and the
    combination closest to the proposal's values, its coefficients rounded to rationals; else (at most 800 columns) the
    exact solution of a square system near a linear-programming stress (v5/rigidity.py). If none is found and the stress
    closest to the uniform one vanishes on some contacts, (iii) and (iv) are tried on its support E2 (Theorem 2.4 allows
    a subset of the contacts): E2 must touch every core point, be connected and pass the rank test itself, and carry an
    exact strictly positive stress constant on classes. The classes and coefficients are untrusted; only the exact
    checks count.
 0. If the cosine is rational and the core's Gram matrix is not recognised over Q, K = Q(sqrt D) is tried for
    D = 2, 3, 5, 6, 7, 10, 11, 13, 14, 15.
Verdicts: STRICT (Theorem 2.4), LOCAL, NOT STRICT (Proposition 2.8), NOT LOCAL, FLAT CORE (Proposition 2.7), or the reason
it is undecided. Each entry runs in a child process whose resident memory is checked every second (stopped above
--mem-mb). The rank test is skipped above --max-cols columns.
"""
import hashlib, json, math, os, sys, time, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, os.path.join(CERT, "v5"))
import warnings
warnings.filterwarnings("ignore")
import numpy as np

SHA = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
NB, NBI = 256, 16


def _unit_lower_inv(L, p):
    """The inverse mod p of the unit lower triangular matrix whose strictly lower part is that of L (float64, entries
    in [0, p)); every dot product has at most 256 terms below 2^32."""
    k = L.shape[0]
    X = np.zeros((k, k))
    for i in range(k):
        X[i, i] = 1.0
        if i:
            X[i, :i] = np.mod(-(L[i, :i] @ X[:i, :i]), p)
    return X


def _check(B, p):
    if not np.all((B >= 0) & (B < p)):
        raise ArithmeticError("an entry left [0, p): the exactness bound was violated")


def rank_mod_p(A, p):
    """Rank over F_p, p < 2^16 prime, of the matrix A (unsigned 16-bit integers in [0, p); modified in place). Blocked
    LU with partial pivoting: panels of 256 columns are factored in double precision (inner blocks of 16 columns whose
    updates are reduced mod p only before an entry is used as a pivot, a multiplier or in a pivot row: at most 256
    pending products below 2^32 each), and the rest of the matrix is updated by products whose dot products have at
    most 256 terms below 2^32, then reduced mod p and stored again as 16-bit integers. Every intermediate value is an
    integer of absolute value below 2^53, so the double-precision arithmetic is exact; every reduced block is checked
    to lie in [0, p)."""
    n, m = A.shape
    assert p < (1 << 16)
    r = 0
    for j0 in range(0, m, NB):
        if r >= n:
            break
        j1 = min(m, j0 + NB)
        r0 = r
        P = A[r0:, j0:j1].astype(np.float64)                 # the panel, rows r0.. (row i of P is row r0 + i of A)
        pc = []
        for i0 in range(0, j1 - j0, NBI):
            if r >= n:
                break
            i1 = min(j1 - j0, i0 + NBI)
            rs = r
            spc = []
            for j in range(i0, i1):
                if r >= n:
                    break
                col = P[r - r0:, j]
                np.mod(col, p, out=col)
                _check(col, p)
                nz = np.flatnonzero(col)
                if len(nz) == 0:
                    continue
                k = r + int(nz[0])
                if k != r:
                    P[[r - r0, k - r0]] = P[[k - r0, r - r0]]
                    if j1 < m:
                        A[[r, k], j1:] = A[[k, r], j1:]
                inv = pow(int(P[r - r0, j]), p - 2, p)
                u = None
                if j + 1 < i1:
                    u = np.mod(P[r - r0, j + 1:i1], p)
                    P[r - r0, j + 1:i1] = u
                if r + 1 < n:
                    f = np.mod(P[r - r0 + 1:, j] * inv, p)
                    P[r - r0 + 1:, j] = f
                    if u is not None:
                        P[r - r0 + 1:, j + 1:i1] -= np.outer(f, u)
                spc.append(j)
                r += 1
            if spc and i1 < j1 - j0:
                X = _unit_lower_inv(P[rs - r0:r - r0, spc], p)
                U = np.mod(X @ np.mod(P[rs - r0:r - r0, i1:], p), p)
                P[rs - r0:r - r0, i1:] = U
                if r < n:
                    P[r - r0:, i1:] -= P[r - r0:, spc] @ U
            pc.extend(spc)
        if pc and j1 < m:
            k = len(pc)
            X = _unit_lower_inv(P[:k, pc], p)
            L21 = P[k:, pc] if r < n else None
            del P
            for c0 in range(j1, m, NB):
                c1 = min(m, c0 + NB)
                U = np.mod(X @ A[r0:r, c0:c1].astype(np.float64), p)
                _check(U, p)
                A[r0:r, c0:c1] = U
                if L21 is not None:
                    T = A[r:, c0:c1].astype(np.float64)
                    T -= L21 @ U
                    np.mod(T, p, out=T)
                    _check(T, p)
                    A[r:, c0:c1] = T
            del L21
    return r


def modp_csr(G, B, E, d, nV, p, rt):
    """The rigidity matrix mod p as a scipy CSR matrix of int64: the rows of certify/v5/rigidity.py R_rows_modp
    (released, unchanged), computed for chunks of the contacts to bound the memory."""
    import scipy.sparse as sps
    import rigidity as RG
    parts = []

    def add(rows, off):
        nnz = sum(len(row) for row in rows)
        ri, ci, vi = np.empty(nnz, np.int64), np.empty(nnz, np.int64), np.empty(nnz, np.int64)
        t = 0
        for i, row in enumerate(rows):
            for j, v in row.items():
                ri[t], ci[t], vi[t] = off + i, j, v
                t += 1
        parts.append((ri, ci, vi))
    add(RG.R_rows_modp(G, B, [], d, nV, p, rt), 0)
    for e0 in range(0, len(E), 4096):
        add(RG.R_rows_modp(G, B, E[e0:e0 + 4096], d, nV, p, rt)[nV:], nV + e0)
    ri, ci, vi = (np.concatenate([q[k] for q in parts]) for k in range(3))
    return sps.csr_matrix((vi, (ri, ci)), shape=(nV + len(E), nV * d))


def kernel_csr(kv, ncols):
    import scipy.sparse as sps
    ri = [i for i, row in enumerate(kv) for _ in row]
    ci = [j for row in kv for j in row]
    vi = [v for row in kv for v in row.values()]
    return sps.csr_matrix((np.array(vi, np.int64), (np.array(ri, np.int64), np.array(ci, np.int64))), shape=(len(kv), ncols))


def gram_modp(M, p):
    """M^T M mod p (the choice S = M^T in Proposition 3.19(c)) for the sparse M (CSR, int64 entries in [0, p)), as
    unsigned 16-bit integers, built in blocks of 256 columns; every entry is a sum of at most (nonzeros of a column of
    M) products below 2^32: exact in int64."""
    ncols = M.shape[1]
    Mt = M.T.tocsr()
    Mc = M.tocsc()
    A = np.zeros((ncols, ncols), dtype=np.uint16)
    for k0 in range(0, ncols, 256):
        k1 = min(ncols, k0 + 256)
        A[:, k0:k1] = np.mod((Mt @ Mc[:, k0:k1]).toarray(), p)
    return A


def sketch(M, ncols, p, seed):
    """S M mod p for the sparse M (CSR, int64 entries in [0, p)), S uniformly random over F_p with ncols rows (or M
    itself when it has at most ncols rows), as a matrix of unsigned 16-bit integers; stored transposed, (S M)^T, which
    has the same rank. Every entry is a sum of at most (nonzeros of a column of M) products below 2^32: exact in int64."""
    nr = M.shape[0]
    if nr <= ncols:
        return np.mod(M.toarray(), p).astype(np.uint16)
    rng = np.random.default_rng(seed)
    Mt = M.T.tocsr()                                        # ncols x nr
    blocks = [Mt[:, i0:min(nr, i0 + 1024)] for i0 in range(0, nr, 1024)]
    del Mt
    A = np.zeros((ncols, ncols), dtype=np.uint16)
    for k0 in range(0, ncols, 256):
        k1 = min(ncols, k0 + 256)
        acc = np.zeros((ncols, k1 - k0), dtype=np.int64)
        for Mb in blocks:
            acc += Mb @ rng.integers(0, p, size=(Mb.shape[1], k1 - k0), dtype=np.int64)
        A[:, k0:k1] = np.mod(acc, p)
        del acc
    return A


def run_entry(d, n, disc, max_cols, log):
    import core
    import rigidity as RG
    from kfield import K, ONE, ksum, prime_for
    t0 = time.time()
    C = core.build(SNAP, d, n, disc)
    if disc == 1 and str(C.get("status", "")).startswith("core Gram matrix not recognised"):
        # a rational cosine with irrational inner products: try K = Q(sqrt D) for small D
        tried = []
        for D2 in (2, 3, 5, 6, 7, 10, 11, 13, 14, 15):
            tried.append(D2)
            C2 = core.build(SNAP, d, n, D2)
            if not str(C2.get("status", "")).startswith("core Gram matrix not recognised"):
                C, disc = C2, D2
                break
        log.update(field_fallback_tried=tried, disc=disc)
    log["seconds_core"] = round(time.time() - t0, 1)
    print("core done", d, n, C.get("status"), log["seconds_core"], flush=True)
    if C["status"] != "ok":
        log.update(status="core: " + C["status"], **{k: v for k, v in C.items() if k != "status"})
        return
    log.update(C["log"])
    G, V, E, c, r, B, W, GBB, GBBinv = (C[k] for k in ("G", "V", "E", "c", "r", "B", "W", "GBB", "GBBinv"))
    nV = len(V)
    if r < d:
        log["status"] = "NOT LOCAL, FLAT CORE (the exact core spans a proper subspace; Proposition 2.7)"
        return
    g = d * (d - 1) // 2
    cols = nV * d
    log["cols"] = cols
    if cols > max_cols:
        log["status"] = f"undecided: rank test above {max_cols} columns"
        return
    comps = RG.components(nV, E)
    kvecs, nflex = [], 0
    if len(comps) > 1:
        rot = RG.rotations(nV, W, GBBinv, d)
        flex = []
        for comp in comps[1:]:
            flex.extend(RG.rotations(nV, W, GBBinv, d, comp=comp))
        for u in rot + flex:
            if any(not x.iszero() for x in RG.apply_R(G, B, E, u, d, nV)):
                raise ValueError("a claimed kernel vector is not in the kernel")
        kvecs, nflex = rot + flex, len(flex)
    # the infinitesimal rotations (of the core, and of each further component) form a space of flexes of dimension
    # g + nflex (the core spans R^d, Lemma 2.2(a)), so the rank over K is at most cols - g - nflex; a rank mod p equal
    # to it proves equality (Proposition 3.19(b,c)): bar-rigid modulo these flexes.
    want = cols - g - nflex
    tries, start = [], 1 << 16
    for seed in (None, 2027):                                   # S = M^T; if that falls short, a random S, another p
        p, rt = prime_for(disc, start=start)
        while True:
            try:
                M = modp_csr(G, B, E, d, nV, p, rt)
                break
            except ValueError:                                  # a denominator divisible by p
                p, rt = prime_for(disc, start=p)
        A = gram_modp(M, p) if seed is None else sketch(M, cols, p, seed)
        del M
        rk = rank_mod_p(A, p)
        del A
        tries.append({"S": "M^T" if seed is None else "random", "prime": p, "seed": seed, "rank": int(rk)})
        p_last, rt_last, start = p, rt, p
        if rk >= want:
            break
    rk = max(t["rank"] for t in tries)
    log["seconds_rank"] = round(time.time() - t0, 1)
    print("rank done", d, n, rk, cols, log["seconds_rank"], flush=True)
    if kvecs:
        kv = [{j: x.modp(p_last, rt_last) for j, x in enumerate(u) if not x.iszero()} for u in kvecs]
        ind = rank_mod_p(sketch(kernel_csr(kv, cols), cols, p_last, None), p_last) if len(kv) <= cols else 0
    else:
        ind = g
    span_exact = (ind == g + nflex) and (rk == want)
    log.update(rank_tries=tries, rank_mod_p=int(rk), kernel_dim_upper=cols - int(rk), components=[len(s_) for s_ in comps],
               kernel_vectors=g + nflex, kernel_vectors_independent=bool(ind == g + nflex), flex_space_exact=bool(span_exact))
    om = None
    unif = [ONE] * len(E)
    if RG.check_stress(W, E, d, nV, unif):
        om, log["stress"] = unif, "uniform"
    elif class_stress(W, E, d, nV, log, neighbour_classes(E, nV), "neighbours") is not None:
        om = log.pop("_class_stress")
        log["stress"] = "constant on classes of contacts (common-neighbour counts), solved exactly"
    elif (lambda kt: kt is not None and class_stress(W, E, d, nV, log, kt[0], "projection", kt[1]) is not None)(proj_classes(C["Xf"][V], E, d, log)):
        om = log.pop("_class_stress")
        log["stress"] = "constant on the classes of the stress closest to the uniform one, solved exactly"
    elif cols <= 4000 and (lambda kt: kt is not None and class_stress(W, E, d, nV, log, kt[0], "lp", kt[1]) is not None)(lp_classes(C["Xf"][V], E, d, log)):
        om = log.pop("_class_stress")
        log["stress"] = "constant on the classes of a linear-programming stress, solved exactly"
    elif cols > 800:
        log["stress"] = "no balanced stress constant on the classes (the exact square solve is not tried above 800 columns)"
    else:
        XV = C["Xf"][V]
        wf, tmax = RG.lp_stress(XV, E, d)
        log["lp_maxmin_stress"] = tmax
        if wf is not None and tmax > 1e-9:
            cand = RG.exact_stress_from_lp(G, B, W, E, d, nV, XV, wf, disc, log)
            if cand is not None and RG.check_stress(W, E, d, nV, cand) and all(x.sign() > 0 for x in cand):
                om = cand
                log["stress"] = "exact (square solve near the linear program's)"
            else:
                log["stress"] = "no exact strictly positive stress found"
        else:
            log["stress"] = "no strictly positive stress (linear program)"
    w = log.pop("_proj_w", None)
    sub = False
    if om is None and span_exact and nflex == 0 and w is not None and w.min() > -1e-9 * float(w.mean()):
        # the stresses vanish on some contacts: (iii) and (iv) on the support E2 of the stress closest to the uniform
        # one (Theorem 2.4 allows a subset of the contacts); E2 must touch every core point and be bar-rigid itself
        mw = float(w.mean())
        Es = [e for e in range(len(E)) if w[e] > 1e-6 * mw]
        E2 = [E[e] for e in Es]
        log["support_contacts"] = len(E2)
        log["support_points_untouched"] = nV - len({a for e_ in E2 for a in e_})
        if 0 < len(E2) < len(E) and log["support_points_untouched"] == 0 and len(RG.components(nV, E2)) == 1:
            rk2 = 0
            for seed in (None, 2027):
                M2 = modp_csr(G, B, E2, d, nV, p_last, rt_last)
                A2 = gram_modp(M2, p_last) if seed is None else sketch(M2, cols, p_last, seed)
                del M2
                rk2 = max(rk2, int(rank_mod_p(A2, p_last)))
                del A2
                if rk2 >= want:
                    break
            log["support_rank_mod_p"] = rk2
            if rk2 == want:
                key2, target2 = cluster(w[Es], 1e-6)
                om2 = class_stress(W, E2, d, nV, log, key2, "support", target2)
                if om2 is not None:
                    log.pop("_class_stress", None)
                    om, sub = om2, True
                    log["stress"] = ("on the support of the stress closest to the uniform one (" + str(len(E2)) + " of "
                                     + str(len(E)) + " contacts), constant on its classes, solved exactly")
    if om is not None and sub:
        log["status"] = ("STRICT (bar-rigid with a strictly positive stress on a subset of the contacts; Theorem 2.4 "
                         "with (iii) and (iv) on the subset)")
    elif om is not None and span_exact and nflex == 0:
        log["status"] = "STRICT (bar-rigid core with a strictly positive stress; Theorem 2.4)"
    elif om is not None and span_exact and len(comps) > 1:
        log["status"] = "LOCAL, NOT STRICT (components each bar-rigid with a strictly positive stress; Proposition 2.8)"
    elif om is None:
        log["status"] = "undecided: " + log["stress"]
    else:
        log["status"] = "undecided: not shown bar-rigid (rank mod p below the bound)"
    log["seconds"] = round(time.time() - t0, 1)


def neighbour_classes(E, nV):
    """UNTRUSTED: classes of contacts by their number of common neighbours."""
    import numpy as np
    A = np.zeros((nV, nV), dtype=np.int32)
    for a, b in E:
        A[a, b] = A[b, a] = 1
    com = A @ A
    return [int(com[a, b]) for a, b in E]


def proj_classes(XV, E, d, log):
    """UNTRUSTED: the equilibrium stress closest to the uniform one (the projection of (1,...,1) onto the stresses, by
    sparse least squares), which is invariant under the symmetries of the code, so constant on the orbits of contacts;
    its values clustered to relative 1e-6. None if it is not strictly positive."""
    import numpy as np
    import scipy.sparse as sps
    from scipy.sparse.linalg import lsqr
    nV = len(XV)
    Q = []
    for a in range(nV):
        q_, _ = np.linalg.qr(np.c_[XV[a], np.eye(d)])
        Q.append(q_[:, 1:d])                                  # an orthonormal basis of the tangent space at x_a
    rows, cols, vals = [], [], []
    for e, (a, b) in enumerate(E):
        for i_, (u, v) in enumerate(((a, b), (b, a))):
            t = Q[u].T @ XV[v]                                # tangent part of x_v at x_u
            for k in range(d - 1):
                rows.append(u * (d - 1) + k)
                cols.append(e)
                vals.append(t[k])
    T = sps.csr_matrix((vals, (rows, cols)), shape=(nV * (d - 1), len(E)))
    z = lsqr(T.T, np.ones(len(E)), atol=1e-14, btol=1e-14, iter_lim=20000)[0]
    w = np.ones(len(E)) - T.T @ z
    log["proj_stress_min_over_mean"] = float(w.min() / w.mean())
    log["proj_stress_residual"] = float(np.abs(T @ w).max())
    log["_proj_w"] = w                                       # for the support of a stress (run_entry); not written
    if not w.min() > 1e-6 * w.mean():
        return None
    return cluster(w, 1e-6)


def cluster(w, rel):
    """UNTRUSTED: classes of the entries of w, consecutive sorted values within relative rel joined; with the class
    means (over the mean of w) as targets."""
    import numpy as np
    order = np.argsort(w)
    key = [0] * len(w)
    c, prev = 0, None
    for e in order:
        if prev is not None and w[e] - prev > rel * max(abs(prev), 1e-12):
            c += 1
        key[int(e)] = c
        prev = w[e]
    sums = {}
    for e in range(len(w)):
        sums.setdefault(key[e], []).append(float(w[e]))
    mean = float(np.mean(w))
    return key, {c_: (sum(v_) / len(v_)) / mean for c_, v_ in sums.items()}


def lp_classes(XV, E, d, log):
    """UNTRUSTED: the max-min stress of the linear program (interior point, so constant on the orbits of the symmetries
    of the problem), its values clustered to relative 1e-7; with the class means as targets."""
    import rigidity as RG
    wf, tmax = RG.lp_stress(XV, E, d)
    log["lp_maxmin_stress"] = tmax
    if wf is None or not tmax > 1e-9:
        return None
    order = sorted(range(len(E)), key=lambda e: wf[e])
    key = [0] * len(E)
    c, prev = 0, None
    for e in order:
        if prev is not None and wf[e] - prev > 1e-7 * max(abs(prev), 1e-12):
            c += 1
        key[e] = c
        prev = wf[e]
    sums = {}
    for e in range(len(E)):
        sums.setdefault(key[e], []).append(float(wf[e]))
    mean = sum(map(sum, sums.values())) / len(E)
    return key, {c_: (sum(v_) / len(v_)) / mean for c_, v_ in sums.items()}


def class_stress(W, E, d, nV, log, key, tag, target=None, max_classes=300):
    """A stress constant on given classes of contacts. UNTRUSTED: the classes, and the coefficients: an exact basis of
    the stresses constant on the classes (the nullspace over Q, by FLINT, of the equilibrium equations with the
    multipliers eliminated; over Q(sqrt D) the system for (x, y), omega = x + y sqrt D, is a doubled system over Q), and
    the combination closest to the target values (default 1), by least squares in double precision with the
    coefficients rounded to rationals. TRUSTED: the stress is then checked exactly (signs over K, equilibrium by
    v5/rigidity.py check_stress). Returns the stress (list over E, K elements) or None; it is also stored in
    log['_class_stress']."""
    import flint
    from fractions import Fraction
    import rigidity as RG
    from kfield import K, ZERO, fq
    classes = sorted(set(key))
    nc = len(classes)
    log["contact_classes_" + tag] = nc
    if nc > max_classes:
        return None
    ci = {c: i for i, c in enumerate(classes)}
    D = K.D
    two = D != 1
    nu = 2 * nc if two else nc
    u = [dict() for _ in range(nV)]               # u[a][c][k] = sum of W[k][b] over the contacts (a, b) of class c
    for e, (a, b) in enumerate(E):
        c = ci[key[e]]
        for s_, t_ in ((a, b), (b, a)):
            v = u[s_].get(c)
            if v is None:
                v = u[s_][c] = [ZERO] * d
            for k in range(d):
                w = W[k][t_]
                if not w.iszero():
                    v[k] = v[k] + w
    rows = []                                     # integer rows over the unknowns (x_c) or (x_c, y_c)
    for a in range(nV):
        wa = [W[k][a] for k in range(d)]
        k0 = next(k for k in range(d) if not wa[k].iszero())
        for k in range(d):
            if k == k0:
                continue
            row = {}
            for c, v in u[a].items():
                x = v[k] * wa[k0] - v[k0] * wa[k]
                if not x.iszero():
                    row[c] = x
            if not row:
                continue
            den = 1
            for x in row.values():
                den = math.lcm(den, int(x.a.q), int(x.b.q))
            if not two:
                rows.append({c: int(x.a.p) * (den // int(x.a.q)) for c, x in row.items()})
            else:
                ra, rb = {}, {}
                for c, x in row.items():
                    A_ = int(x.a.p) * (den // int(x.a.q))
                    B_ = int(x.b.p) * (den // int(x.b.q))
                    # (A + B sqrt D)(X + Y sqrt D) = (A X + D B Y) + (B X + A Y) sqrt D
                    ra[c], ra[nc + c] = A_, D * B_
                    rb[c], rb[nc + c] = B_, A_
                rows.append(ra)
                rows.append(rb)
    if not rows:
        return None
    Bm = flint.fmpz_mat(len(rows), nu)
    for i, row in enumerate(rows):
        for c, v in row.items():
            if v:
                Bm[i, c] = v
    del rows, u
    N, k = (Bm.transpose() * Bm).nullspace()      # over Q, the nullspace of B^T B is that of B
    log["class_stress_space_dim_" + tag] = int(k)
    if k == 0:
        return None
    sq = math.sqrt(D)
    basis = [[int(N[i, j]) for i in range(nu)] for j in range(k)]

    def fl(v, e):                                 # approximately v / 2^e, for integers of any size
        sh = max(0, abs(v).bit_length() - 62)
        return math.ldexp(float(abs(v) >> sh) * (1 if v >= 0 else -1), sh - e)
    ex = [max(abs(v).bit_length() for v in bv) for bv in basis]   # basis vector j is used as bv / 2^ex[j]
    R = np.array([[fl(bv[c], e) + (fl(bv[nc + c], e) * sq if two else 0.0) for bv, e in zip(basis, ex)]
                  for c in range(nc)])
    t = np.array([float(target.get(c, 1.0)) if target else 1.0 for c in classes])
    y = np.linalg.lstsq(R, t, rcond=None)[0]
    coef = [Fraction(float(y_)).limit_denominator(10 ** 12) / (1 << e) for y_, e in zip(y, ex)]
    f = []
    for c in range(nc):
        X = sum((q_ * bv[c] for q_, bv in zip(coef, basis)), Fraction(0))
        Y = sum((q_ * bv[nc + c] for q_, bv in zip(coef, basis)), Fraction(0)) if two else Fraction(0)
        f.append(K(fq(X.numerator, X.denominator), fq(Y.numerator, Y.denominator)))
    if all(x.sign() < 0 for x in f):
        f = [-x for x in f]
    if not all(x.sign() > 0 for x in f):
        log["class_stress_" + tag] = "the closest class-constant stress has entries of both signs or zeros"
        return None
    om = [f[ci[key[e]]] for e in range(len(E))]
    if not RG.check_stress(W, E, d, nV, om):
        log["class_stress_" + tag] = "the exact check failed"
        return None
    reps = [repr(x) for x in f]
    if sum(map(len, reps)) <= 4000:
        log["class_stress_values"] = reps
    else:
        log["class_stress_values_float"] = [int(x.a.p) / int(x.a.q) + (int(x.b.p) / int(x.b.q) * sq if two else 0.0)
                                            for x in f][:60]
        log["class_stress_values_sha256"] = hashlib.sha256("\n".join(reps).encode()).hexdigest()
        log["class_stress_values_chars"] = sum(map(len, reps))
    log["class_stress_from"] = tag
    log["_class_stress"] = om
    return om


def _child(conn, d, n, disc, max_cols):
    log = {}
    try:
        run_entry(d, n, disc, max_cols, log)
    except MemoryError:
        log["status"] = "stopped: MemoryError"
    except Exception:
        log["status"] = "error"
        log["error"] = traceback.format_exc()[-1500:]
    log.pop("_proj_w", None)
    log.pop("_class_stress", None)
    conn.send(json.dumps(log, default=str))
    conn.close()


def parse_poly(s):
    """Coefficients {degree: integer} of an integer polynomial in x written like '111*x^2+2*x-1'."""
    import re
    out = {}
    for sign, body in re.findall(r"([+-]?)([^+-]+)", s.replace(" ", "")):
        sg = -1 if sign == "-" else 1
        if "x" in body:
            cpart, _, epart = body.partition("x")
            cf = int(cpart.rstrip("*")) if cpart.rstrip("*") else 1
            e = int(epart.lstrip("^")) if epart.startswith("^") else 1
        else:
            cf, e = int(body), 0
        out[e] = out.get(e, 0) + sg * cf
    return out


def squarefree(D):
    s = 1 if D > 0 else -1
    D = abs(D)
    out, f = 1, 2
    while f * f <= D:
        while D % (f * f) == 0:
            D //= f * f
        if D % f == 0:
            out *= f
            D //= f
        f += 1
    return s * out * D


def main():
    import multiprocessing, subprocess
    import snapshot
    outd = sys.argv[1]
    args = sys.argv[2:]
    only = set(args[args.index("--only") + 1].split(",")) if "--only" in args else None
    max_cols = int(args[args.index("--max-cols") + 1]) if "--max-cols" in args else 20000
    mem = int(args[args.index("--mem-mb") + 1]) if "--mem-mb" in args else 1500
    os.makedirs(outd, exist_ok=True)
    s5 = json.load(open(os.path.join(CERT, "v5", "out", "summary_v5.json")))["class"]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index() if e["listed"]}
    todo = []
    for k, v in s5.items():
        d, n = map(int, k.split(","))
        if v[0] == "none_size" and idx[(d, n)]["poly"]:
            todo.append((d * n, d, n, idx[(d, n)]["poly"]))
    todo.sort()
    for _, d, n, poly in todo:
        if only and f"{d}_{n}" not in only:
            continue
        js = os.path.join(outd, f"d{d}_n{n}.json")
        if os.path.exists(js):
            continue
        coef = parse_poly(poly)
        deg = max(coef)
        log = {"d": d, "n": n, "poly": poly, "degree": deg, "program_sha256": SHA}
        t0 = time.time()
        if deg > 2:
            log["status"] = "not attempted: cosine of degree above 2"
        else:
            if deg == 1:
                disc = 1
            else:
                a_, b_, c_ = coef.get(2, 0), coef.get(1, 0), coef.get(0, 0)
                disc = squarefree(b_ * b_ - 4 * a_ * c_)
            log["disc"] = disc
            ctx = multiprocessing.get_context("spawn")
            rd, wr = ctx.Pipe(duplex=False)
            pr = ctx.Process(target=_child, args=(wr, d, n, disc, max_cols))
            pr.start()
            wr.close()
            peak = 0
            got = None
            while True:
                if rd.poll(1.0):
                    got = json.loads(rd.recv())
                    break
                if not pr.is_alive():
                    break
                o = subprocess.run(["ps", "-o", "rss=", "-p", str(pr.pid)], capture_output=True, text=True).stdout.strip()
                rss = int(o) // 1024 if o else 0
                peak = max(peak, rss)
                if rss > mem:
                    pr.kill()
                    break
            pr.join()
            if got is None:
                got = {"status": f"stopped: memory above {mem} MB" if peak > mem else "stopped: child process ended"}
            log.update(got)
            log["peak_rss_mb_sampled"] = peak
        log["wall"] = round(time.time() - t0, 1)
        with open(js + ".part", "w") as f:
            json.dump(log, f, default=str)
        os.replace(js + ".part", js)
        print(d, n, log.get("status"), log["wall"], flush=True)


if __name__ == "__main__":
    main()
