"""Census certifier, version 4: the certificates of the paper (version 4), Sections 3 and 4.

  python3 -I certify/census/cq.py CODE.txt D N [--tol 1e-10] [--json OUT.json] [--flex]

ISOSTATIC contact graph (|E| = |V|(d-1) - d(d-1)/2 + 1): the existence theorem and its quantitative form. An exact
configuration P* exists within the printed distance of the code, its contacts are exactly E, its largest cosine c* is
enclosed, its stress (unique up to scale) is strictly positive and its contact framework bar-rigid, so P* is a strict
local minimum of the largest cosine modulo rotations and rattlers (status PASS); with a sharpness constant kappa and an
explicit radius rho, and the radius r_u of epsilon-inflation. If the existence, span and gap tests pass and the
stress has a provably negative entry, P* is not a local minimum (status NONOPT).

HYPERSTATIC contact graph (|E| larger): the hyperstatic theorem. A subset E' of the isostatic size is chosen
(untrusted), P' is
proved to exist for E', the other contacts are enclosed within tau' of c', and a strictly positive stress of the whole
of E at P' is enclosed: a local minimum lies within an explicit distance of P' (status LOCALMIN). Whether it is strict
is decided, when it can be, by strict.py.

FLEXIBILITY of a contact set (option --flex): the flexibility proposition. E = all pairs within the tolerance of the
largest cosine
(no rattlers removed). It is proved that the rows of the norm and contact equations are linearly independent at every
configuration within an explicit distance h of the code (|q_i - x_i|_1 < h for every point): no configuration there
carries a nonzero equilibrium stress on E, so none whose contact graph lies in E is a local minimum (status FLEX).

The stored code is the stored decimal coordinates normalised exactly; every distance to it that is printed includes
the error of reading and normalising it in double precision (stored_norm_error).

TRUSTED / UNTRUSTED as in sphcert.py. The trusted part uses only exact integer arithmetic, elementwise numpy
operations, numpy sums and scipy's sparse products (Lemmas R and S of sphcert.py; the rounding lemma of the paper); no
dense
BLAS product enters a decision. The approximate inverse may be stored in single precision for large systems (it is any
matrix; every block that enters a decision is converted exactly to double precision). Every result carries the SHA-256
of cq.py and sphcert.py (code_sha256).
"""
import argparse, hashlib, json, math, os, sys, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.dirname(HERE))  # certify/ (sphcert.py)
import numpy as np
import scipy.sparse as sps
from scipy.linalg import lapack

import sphcert as S
from sphcert import up, gam, frac_up, frac_down, U, ETA

BLOCK = 128
DEN = 1.0 - 2.0 ** -50   # (1 - beta) rounded down: multiply by this
VERSION = 4


def _sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


CODE_SHA256 = {"cq.py": _sha256(os.path.abspath(__file__)), "sphcert.py": _sha256(S.__file__)}


# ============================================================================================ reading
def read_code(path, d, n):
    """Points as doubles, normalised in double precision (the stored code X of the paper)."""
    rows = []
    for line in open(path):
        s = line.strip()
        parts = s.split(",")
        if len(parts) != d:
            continue
        try:
            rows.append([float(t) for t in parts])
        except ValueError:
            continue
    X = np.array(rows, dtype=float)
    if X.shape != (n, d):
        raise SystemExit(f"{path}: {len(rows)} points of dimension {d}, expected {n}")
    return X / np.linalg.norm(X, axis=1, keepdims=True)


def read_decimals(path, d, n):
    """The stored coordinates as strings (the lines of d comma-separated numbers)."""
    rows = []
    for line in open(path):
        parts = [t.strip() for t in line.strip().split(",")]
        if len(parts) != d:
            continue
        try:
            [float(t) for t in parts]
        except ValueError:
            continue
        rows.append(parts)
    assert len(rows) == n
    return rows


def _sqrt_bounds(q, bits=160):
    """TRUSTED. Rationals lo <= sqrt(q) <= hi for a rational q > 0, from integer square roots."""
    a, b = q.numerator, q.denominator
    s = math.isqrt(a * b << (2 * bits))
    return Fraction(s, b << bits), Fraction(s + 1, b << bits)


def stored_norm_error(dec, xs):
    """TRUSTED. dec: stored coordinates as decimal strings (k x d); xs: the same points as the certifier uses them
    (rounded to double precision and normalised in double precision). Returns floats e (k x d) with
    |xhat_ik - xs_ik| <= e_ik, where xhat_i = (stored coordinates of point i) / |.|, normalised exactly: the stored code
    of the paper, an exact point of the sphere."""
    out = np.empty(xs.shape)
    for i, row in enumerate(dec):
        q = [Fraction(t) for t in row]
        s2 = sum(x * x for x in q)
        lo, hi = _sqrt_bounds(s2)
        for k, x in enumerate(q):
            f = Fraction(float(xs[i, k]))
            a, b = x / hi, x / lo            # xhat_ik lies between a and b
            out[i, k] = frac_up(max(abs(a - f), abs(b - f)))
    return out


# ============================================================================================ UNTRUSTED
def approx_inverse(Jc, dtype):
    """UNTRUSTED. Inverse by LAPACK getrf + getri in place (one dense array of the given dtype)."""
    n = Jc.shape[0]
    A = np.zeros((n, n), dtype=dtype, order="F")
    co = Jc.tocoo()
    A[co.row, co.col] = co.data
    if dtype == np.float64:
        getrf, getri, lw = lapack.dgetrf, lapack.dgetri, lapack.dgetri_lwork
    else:
        getrf, getri, lw = lapack.sgetrf, lapack.sgetri, lapack.sgetri_lwork
    lu, piv, info = getrf(A, overwrite_a=True)
    if info != 0:
        raise np.linalg.LinAlgError("singular")
    lwork, _ = lw(n)
    inv, info = getri(lu, piv, lwork=int(lwork), overwrite_lu=True)
    if info != 0:
        raise np.linalg.LinAlgError("getri")
    for i0 in range(0, n, 4096):
        if not np.all(np.isfinite(inv[:, i0:i0 + 4096])):
            raise np.linalg.LinAlgError("not finite")
    return inv


def rows64(R, i0, i1):
    return np.asarray(R[i0:i1], dtype=np.float64)


def cols64(R, j0, j1):
    return np.asarray(R[:, j0:j1], dtype=np.float64)


def matvec(R, v):
    """UNTRUSTED. R v in double precision, by row blocks (R may be single precision)."""
    out = np.empty(R.shape[0])
    for i0 in range(0, R.shape[0], 2048):
        out[i0:i0 + 2048] = rows64(R, i0, i0 + 2048) @ v
    return out


def rmatvec(R, v):
    """UNTRUSTED. R^T v in double precision."""
    out = np.empty(R.shape[1])
    for j0 in range(0, R.shape[1], 2048):
        out[j0:j0 + 2048] = cols64(R, j0, j0 + 2048).T @ v
    return out


def refine(ex, Z, C, R, K, iters=14, log=None):
    """UNTRUSTED. Newton's chord method with exact residuals, rounded to the grid 2^-K."""
    s = ex.s
    hist = []
    for _ in range(iters):
        Ff, _ = ex.residual(Z, C)
        res = float(np.abs(Ff).max())
        hist.append(res)
        if res < 2.0 ** -(K - 6) or (len(hist) > 2 and res > 0.5 * hist[-2]):
            break
        dz = np.rint(-matvec(R, Ff) * 2.0 ** K)
        for a in range(s.n):
            row = Z[a]
            for k in range(s.d):
                row[k] += int(dz[a * s.d + k])
        C += int(dz[s.n * s.d])
    if log is not None:
        log["newton_residuals"] = [float(f"{h:.3e}") for h in hist]
    return Z, C, hist[-1]


def choose_square_subset_qr(sysm_full, X0, d, n, m_iso):
    """UNTRUSTED. Rank-revealing choice of E': the contact rows of the Jacobian are projected onto the orthogonal
    complement of the row space of the norm and slice rows, and QR with column pivoting (LAPACK geqp3, single
    precision) on the transposed projected rows picks m_iso of them."""
    import scipy.linalg as sl
    J = sysm_full.jac(X0).tocsr()
    m = sysm_full.m
    B = sps.vstack([J[:n], J[n + m:]]).toarray()           # norm and slice rows
    Qb, _ = np.linalg.qr(B.T)                              # orthonormal basis of their row space (N x (n+g))
    del B
    Cr = J[n:n + m].toarray().astype(np.float32)           # contact rows, m x N
    for i0 in range(0, m, 1024):                           # project in blocks (no m x N temporaries)
        blk = Cr[i0:i0 + 1024].astype(np.float64)
        Cr[i0:i0 + 1024] = blk - (blk @ Qb) @ Qb.T
    _, piv = sl.qr(Cr.T, mode="r", pivoting=True, overwrite_a=True, check_finite=False)
    return sorted(int(e) for e in piv[:m_iso])


# ============================================================================================ TRUSTED
def krawczyk(sysm, Ff, Ferr, Jc, R, dfx, log):
    """TRUSTED. Rump 2010, Thm 13.3 with X = [-r, r]^N, through the componentwise bound Ybar + r rs(r) < r, where
    Ybar >= |R F(z~)| and rs(r) >= the row sums of |I - R J| over J in J(Xf) +- Jr(dfx + r). The row sums are split
    as rs(r) <= A + r B (A: |I - R Jc| + gamma |R||Jc| + |R| Jr(dfx); B: |R| Jr(1)), so the test can be repeated for
    any r (epsilon-inflation). Also returns upper bounds of the row sums of |R| (all columns; the columns of the C and
    G rows) for the bound on the inverse of the Jacobian."""
    n, d, m = sysm.n, sysm.d, sysm.m
    N = sysm.nu
    g = gam(N + 1)
    absJc = abs(Jc)
    uJ = up(np.asarray(absJc.sum(1)).ravel())            # row sums of |Jc|
    a = np.zeros(N)                                      # row sums of Jr(dfx)
    a[:n] = up(2.0 * up(dfx.sum(1)))
    sd = up(dfx.sum(1))
    a[n:n + m] = up(sd[sysm.Ea] + sd[sysm.Eb])
    jr1 = np.zeros(N)                                    # row sums of Jr(1): 2d on the S and C rows
    jr1[:n + m] = 2.0 * d
    Ybar = np.empty(N)
    rowA = np.empty(N)
    rowB = np.empty(N)
    nR = np.empty(N)
    nRcg = np.empty(N)
    for i0 in range(0, N, BLOCK):
        Rb = rows64(R, i0, i0 + BLOCK)
        aR = np.abs(Rb)
        RF = (Rb * Ff[None, :]).sum(1)
        A_ = (aR * np.abs(Ff)[None, :]).sum(1)
        B_ = (aR * Ferr[None, :]).sum(1)
        Ybar[i0:i0 + BLOCK] = up(np.abs(RF) + up(g * up(A_)) + up(B_))
        P = np.asarray(Rb @ Jc)
        idx = np.arange(P.shape[0])
        Eh = -P
        Eh[idx, i0 + idx] = 1.0 - P[idx, i0 + idx]
        t1 = up(up(np.abs(Eh).sum(1)))
        t2 = up(g * up((aR * uJ[None, :]).sum(1)))
        t3 = up((aR * a[None, :]).sum(1))
        rowA[i0:i0 + BLOCK] = up(t1 + t2 + t3)
        rowB[i0:i0 + BLOCK] = up(jr1.max() * up(aR[:, :n + m].sum(1)))
        nR[i0:i0 + BLOCK] = up(aR.sum(1))
        nRcg[i0:i0 + BLOCK] = up(aR[:, n:].sum(1))
    out = {"Ybar": Ybar, "rowA": rowA, "rowB": rowB, "nR": nR, "nRcg": nRcg}

    def test(r):
        rs = up(rowA + up(r * rowB))
        err = up(Ybar + up(rs * r))
        return bool(np.all(err < r)), rs, err

    r = float(up(2.0 * Ybar.max()))
    ok, rs, err = test(r)
    att = [{"r": r, "max_rowsum": float(rs.max()), "max_err_over_r": float((err / r).max())}]
    k = 0
    while not ok and rs.max() < 1.0 and k < 3:
        r = float(up(4.0 * Ybar.max() / (1.0 - float(rs.max()))))
        ok, rs, err = test(r)
        att.append({"r": r, "max_rowsum": float(rs.max()), "max_err_over_r": float((err / r).max())})
        k += 1
    log["krawczyk_attempts"] = att
    out.update(ok=ok, r=r, rs=rs, err=err)
    if ok:  # epsilon-inflation: the largest radius (within a factor 2) for which the test still passes
        q = float(np.min((1.0 - rowA) / np.maximum(rowB, 1e-300)))
        ru, best = 0.5 * q, None
        for _ in range(60):
            if ru <= r:
                break
            ok2, _, _ = test(ru)
            if ok2:
                best = ru
                break
            ru *= 0.5
        out["r_unique"] = best if best is not None else r
    return out


def interval_jac_radius(sysm, rad):
    """TRUSTED. The radius matrix of the Jacobian when every coordinate (a, k) is within rad[a, k] of Xf."""
    return sysm.jac_radius(rad)


def transposed_enclosure(sysm, Jc, Jr, R, b, brad, log, key):
    """TRUSTED (y~ untrusted). Encloses y = J^-T b for every J in Jc +- Jr and every right-hand side within brad of b:
    |y - y~| <= w + cs ||y - y~||_inf, with w >= |R|^T (|b - J^T y~| + brad) and cs >= the column sums of |I - J R|.
    Returns (y~, bound) or (None, None) when beta = max cs >= 1."""
    N = sysm.nu
    g = gam(N + 1)
    yt = rmatvec(R, b)
    for _ in range(3):
        yt = yt + rmatvec(R, b - Jc.T @ yt)
    absJc = abs(Jc)
    rh = b - Jc.T @ yt
    S1 = absJc.T @ np.abs(yt)
    S2 = Jr.T @ np.abs(yt)
    Resbar = up(up(up(np.abs(rh)) + up(g * up(S1)) + up(S2)) + brad)
    vJ = up(np.asarray(absJc.sum(0)).ravel())   # column sums of |Jc|
    vJr = up(np.asarray(Jr.sum(0)).ravel())     # column sums of Jr
    w = np.empty(N)
    cs = np.empty(N)
    for j0 in range(0, N, BLOCK):
        Rc = cols64(R, j0, j0 + BLOCK)
        aR = np.abs(Rc)
        w[j0:j0 + BLOCK] = up((aR * Resbar[:, None]).sum(0))
        Q = np.asarray(Jc @ Rc)
        idx = np.arange(Q.shape[1])
        Eh = -Q
        Eh[j0 + idx, idx] = 1.0 - Q[j0 + idx, idx]
        t1 = up(up(np.abs(Eh).sum(0)))
        t2 = up(g * up((aR * vJ[:, None]).sum(0)))
        t3 = up((aR * vJr[:, None]).sum(0))
        cs[j0:j0 + BLOCK] = up(t1 + t2 + t3)
    beta = float(cs.max())
    log[key + "_beta"] = beta
    if not beta < 1.0:
        return None, None
    dinf = up(up(w.max()) / ((1.0 - beta) * DEN))
    return yt, up(w + up(cs * dinf))


def exact_pair_cos(Z, K, a, b):
    return Fraction(sum(x * y for x, y in zip(Z[a], Z[b])), 1 << (2 * K))


# ============================================================================================ driver
def certify(X, d, n, tol=1e-10, K=100, log=None, dtype=None, contacts=None, dump=None, structure_out=None, dec=None,
            keep=None):
    """X: the code (n x d doubles, normalised in double precision); dec: the stored coordinates as decimal strings
    (then every distance printed is to the exactly normalised stored code). Returns the log dict with 'status' among
    PASS, NONOPT, LOCALMIN, and failure strings. keep: a dict that receives the data the strictness stage needs."""
    log = {} if log is None else log
    log["code_sha256"] = CODE_SHA256
    cpu0, t0 = time.process_time(), time.time()
    try:
        return _certify(X, d, n, tol, K, log, dtype, contacts, dump, structure_out, dec, keep)
    finally:
        log["cpu_seconds"] = round(time.process_time() - cpu0, 2)
        log["wall_seconds"] = round(time.time() - t0, 2)


def choose_beta(Esq, extra, Xf, R, N, nV, m, d, log):
    """UNTRUSTED. Coefficients beta_f >= 0 of the extra contacts of a hyperstatic graph: from approximate solutions
    y^f = J^-T grad F_f of the transposed systems, the combined stress closest to uniform (least squares in beta), or,
    if that is not positive, a linear program that maximises its smallest entry. Returns beta or None."""
    from scipy.optimize import linprog
    b = np.zeros(N)
    b[N - 1] = -1.0
    rr, cc, vv = [], [], []
    for f, (a_, b_) in enumerate(extra):
        rr += [f] * (2 * d + 1)
        cc += list(range(a_ * d, (a_ + 1) * d)) + list(range(b_ * d, (b_ + 1) * d)) + [N - 1]
        vv += list(Xf[b_]) + list(Xf[a_]) + [-1.0]
    Jf = sps.csr_matrix((vv, (rr, cc)), shape=(len(extra), N))
    y0 = rmatvec(R, b)
    Y = np.empty((len(extra), N))                       # rows J_f R = (R^T J_f^T)^T, by column blocks of R
    for j0 in range(0, N, 2048):
        Y[:, j0:j0 + 2048] = np.asarray(Jf @ cols64(R, j0, j0 + 2048))
    om0 = y0[nV:nV + len(Esq)]
    Ysub = Y[:, nV:nV + len(Esq)]
    kx = len(extra)
    u = 1.0 / m
    Amat = Ysub @ Ysub.T + np.eye(kx)
    rhs = Ysub @ (om0 - u) + u
    try:
        beta_ls = np.linalg.solve(Amat, rhs)
    except np.linalg.LinAlgError:
        beta_ls = None
    if beta_ls is not None and beta_ls.min() > 0 and (om0 - Ysub.T @ beta_ls).min() > 0:
        log["stress_choice"] = "least squares"
        log["stress_lp_margin"] = float(min(beta_ls.min(), (om0 - Ysub.T @ beta_ls).min()))
        return beta_ls
    top = sps.csr_matrix(np.hstack([Ysub.T, np.ones((len(Esq), 1))]))
    bot = sps.hstack([-sps.identity(kx, format="csr"), sps.csr_matrix(np.ones((kx, 1)))], format="csr")
    A_ub = sps.vstack([top, bot], format="csr")
    b_ub = np.concatenate([om0, np.zeros(kx)])
    lp = linprog(np.r_[np.zeros(kx), -1.0], A_ub=A_ub, b_ub=b_ub, bounds=[(0, None)] * kx + [(None, 1.0)],
                 method="highs-ipm", options={"time_limit": 120.0})
    log["stress_lp_status"] = int(lp.status)
    log["stress_choice"] = "linear program"
    if lp.status != 0 or lp.x is None or lp.x[-1] <= 0:
        log["stress_lp_margin"] = float(lp.x[-1]) if (lp.status == 0 and lp.x is not None) else None
        return None
    log["stress_lp_margin"] = float(lp.x[-1])
    return lp.x[:kx].copy()


def extra_rhs(extra, beta_f, Xf, radz, N, d):
    """TRUSTED. The right-hand side b = -e_gamma - sum_f beta_f grad F_f(z) of the combined stress, at Xf, with a
    componentwise radius brad that covers the rounding of b and every z whose points are within radz of Xf."""
    Bv = np.zeros(N)
    Babs = np.zeros(N)
    Brad = np.zeros(N)
    for f, (a_, b_) in enumerate(extra):
        Bv[a_ * d:(a_ + 1) * d] -= beta_f[f] * Xf[b_]
        Bv[b_ * d:(b_ + 1) * d] -= beta_f[f] * Xf[a_]
        Babs[a_ * d:(a_ + 1) * d] += beta_f[f] * np.abs(Xf[b_])
        Babs[b_ * d:(b_ + 1) * d] += beta_f[f] * np.abs(Xf[a_])
        Brad[a_ * d:(a_ + 1) * d] += beta_f[f] * radz[b_]
        Brad[b_ * d:(b_ + 1) * d] += beta_f[f] * radz[a_]
    b = Bv
    b[N - 1] = -1.0 + float(np.sum(beta_f))
    Babs[N - 1] = 1.0 + float(np.sum(beta_f))
    kmax = max(1, len(extra) + 1)
    brad = up(up(gam(min(kmax, 2 ** 20)) * up(Babs)) + up(Brad))
    return b, brad


def _certify(X, d, n, tol, K, log, dtype, contacts, dump=None, structure_out=None, dec=None, keep=None):
    V, E, M, _ = S.contact_structure(X, d, tol)
    nV, m = len(V), len(E)
    iso = nV * (d - 1) - d * (d - 1) // 2 + 1
    log.update(d=d, n=n, tol=tol, bits=K, n_nonrattlers=int(nV), contacts=int(m), isostatic_count=int(iso),
               rattlers=int(n - nV), max_cos_double=M)
    if nV == 0 or m < iso:
        log["status"] = "no contact graph" if nV == 0 else "hypostatic: contacts - isostatic count = %d" % (m - iso)
        return log
    xs = X[V]
    hyper = m > iso
    Eall = E
    if hyper:
        log["kind"] = "hyperstatic"
        full = S.System(xs, E, d)
        sub = choose_square_subset_qr(full, xs, d, nV, iso)
        if sub is None:
            log["status"] = "hyperstatic: no square subsystem chosen"
            return log
        Esq = [E[e] for e in sub]
        extra = [E[e] for e in range(m) if e not in set(sub)]
    else:
        log["kind"] = "isostatic"
        Esq, extra = E, []
    sysm = S.System(xs, Esq, d)
    N = sysm.nu
    if structure_out:
        json.dump({"V": [int(v) for v in V], "E": [[int(V[a_]), int(V[b_])] for a_, b_ in Esq]},
                  open(structure_out, "w"))
    log["unknowns"] = int(N)
    if dtype is None:
        dtype = np.float64 if N <= 12000 else np.float32
    log["inverse_dtype"] = "float64" if dtype == np.float64 else "float32"
    c0 = float(np.mean([xs[a] @ xs[b] for a, b in Esq]))
    Z = [[int(v) for v in row] for row in (xs * 2.0 ** K)]
    C = int(c0 * 2.0 ** K)
    ex = S.Exact(sysm, K)
    t1 = time.time()
    try:
        R = approx_inverse(sysm.jac(xs), dtype)
    except (np.linalg.LinAlgError, ValueError):
        log["status"] = "Jacobian numerically singular"
        return log
    log["seconds_inverse"] = round(time.time() - t1, 2)
    Z, C, res = refine(ex, Z, C, R, K, log=log)
    # ---------------------------------------------------------------- trusted part
    t2 = time.time()
    Ff, Ferr = ex.residual(Z, C)
    Xf = np.array([[z / (1 << K) for z in row] for row in Z])
    dfx = np.abs(Xf) * 2.0 ** -52 + ETA
    Jc = sysm.jac(Xf)
    kr = krawczyk(sysm, Ff, Ferr, Jc, R, dfx, log)
    log["krawczyk"] = kr["ok"]
    log["residual_max"] = float(np.abs(Ff).max())
    log["seconds_krawczyk"] = round(time.time() - t2, 2)
    if not kr["ok"]:
        log["status"] = "FAIL existence (Krawczyk)"
        return log
    r, err, rs = kr["r"], kr["err"], kr["rs"]
    beta_k = float(rs.max())
    err_pts = up(dfx + err[:nV * d].reshape(nV, d))
    c_err = float(err[-1])
    cq = Fraction(C, 1 << K)
    c_lo, c_hi = frac_down(cq - Fraction(c_err)), frac_up(cq + Fraction(c_err))
    log.update(c_mid=str(C * 10 ** 40 // (1 << K)), c_err=c_err, coord_err_max=float(err[:-1].max()),
               krawczyk_r=r, krawczyk_beta=beta_k, r_unique=kr["r_unique"],
               krawczyk_margin=float((err / r).max()),
               R_norm_inf=float(kr["nR"].max()), J_norm_inf=float(up(np.asarray(abs(Jc).sum(1)).max())))
    log["cond_inf_upper"] = float(up(log["R_norm_inf"] * log["J_norm_inf"] / ((1 - beta_k) * DEN)))
    # ---- stress (and, for hyperstatic codes, the combination with the extra contacts)
    rad = up(dfx + r)
    Jr = interval_jac_radius(sysm, rad)
    b = np.zeros(N)
    b[N - 1] = -1.0
    brad = np.zeros(N)
    beta_f = None
    if hyper:
        beta_f = choose_beta(Esq, extra, Xf, R, N, nV, m, d, log)
        if beta_f is None:
            log["status"] = "hyperstatic: no positive stress found"
            return log
        # b = -e_c - sum_f beta_f J_f(z')^T, with its radius (rounding and the enclosure of z')
        radz = up(dfx + err[:nV * d].reshape(nV, d))
        b, brad = extra_rhs(extra, beta_f, Xf, radz, N, d)
    t3 = time.time()
    yt, ybound = transposed_enclosure(sysm, Jc, Jr, R, b, brad, log, "stress")
    log["seconds_stress"] = round(time.time() - t3, 2)
    if yt is None:
        log["status"] = "FAIL stress enclosure"
        return log
    msq = len(Esq)
    om = yt[nV:nV + msq]
    ob = ybound[nV:nV + msq]
    om_lo = om - ob  # rounding: om - ob is a float difference; take one ulp down to be safe
    om_lo = np.nextafter(om_lo, -np.inf)
    om_hi = np.nextafter(om + ob, np.inf)
    if hyper:
        wmin = float(min(om_lo.min(), beta_f.min()))
        stress_pos = bool(om_lo.min() > 0 and beta_f.min() > 0)
        stress_neg = False
        mtot = m
    else:
        wmin = float(om_lo.min())
        stress_pos = bool(om_lo.min() > 0)
        stress_neg = bool(om_hi.min() < 0)
        mtot = msq
    log.update(stress_min_lower=wmin, stress_min_rel=wmin * mtot, stress_positive=stress_pos,
               stress_provably_negative_entry=stress_neg, stress_max=float(om.max()),
               stress_err_max=float(ob.max()),
               stress_quantiles_over_mean=[float(x) for x in np.quantile(
                   (np.concatenate([om, beta_f]) if hyper else om) * mtot, [0, 0.01, 0.1, 0.5, 0.9, 0.99, 1])])
    # ---- span
    ok_span, span_rs, _ = S.span_check(Xf, err_pts, d)
    log.update(span=ok_span, span_rowsum=span_rs)
    # ---- the extra contacts of a hyperstatic code: |<p'_a, p'_b> - c'| <= tau'
    e_max = Fraction(float(up(err.max())))   # |p' - z~| <= err (z~ = Z / 2^K exactly; Xf plays no role here)
    tau = Fraction(0)
    if hyper:
        for (a_, b_) in extra:
            s_ = exact_pair_cos(Z, K, a_, b_)
            l1 = Fraction(sum(abs(x) for x in Z[a_]) + sum(abs(y) for y in Z[b_]), 1 << K)
            dev = abs(s_ - cq) + e_max * l1 + d * e_max * e_max + Fraction(c_err)
            tau = max(tau, dev)
        log["tau_extra"] = float(frac_up(tau))
    # ---- rattlers and gaps (the extra contacts count as contacts; the margin covers tau' and the move of Q0)
    Y, rinfo = S.place_rattlers(X, V, Xf, float(cq), d)
    log["rattler_placement"] = rinfo
    E_glob = sorted((int(V[a_]), int(V[b_])) for a_, b_ in Eall)
    # ---- sharpness constant and radius (the rigidity lemma and the quantitative theorem of the paper)
    mu = float(up(log["R_norm_inf"] / ((1.0 - beta_k) * DEN)))
    rows_p = slice(0, nV * d)
    mu_p = float(up(np.max(up(kr["nRcg"][rows_p] + up(rs[rows_p] * mu)))))
    dxi = float(up(np.abs(Xf - xs).max() + err_pts.max()))   # max |xi - p*| (xi = the stored code, the slice reference)
    slice_term = float(up(2.0 * nV * dxi * mu_p))
    log.update(mu_inf=mu, mu_p=mu_p, slice_term=slice_term)
    kappa = rho = None
    if stress_pos and slice_term < 1.0:
        s_rig = (1.0 - slice_term) * DEN / mu_p
        s_rig = float(np.nextafter(s_rig, 0))
        kappa = float(np.nextafter(s_rig * wmin / (1.0 - wmin) * DEN, 0))
        lam = float(up(up(1.0 + max(abs(float(c_lo)), abs(float(c_hi)))) + float(frac_up(tau))))
        rho = float(np.nextafter(min(kappa / (4.0 * lam * math.sqrt(d) * (1 + 2 ** -40)), 0.5), 0))
        log.update(rigidity_s=s_rig, kappa=kappa, kappa_point=float(np.nextafter(kappa / (math.sqrt(d) * (1 + 2 ** -40)), 0)),
                   rho=rho, lam=lam)
    # distance of the stored code from P*: max coordinate and least squares. With dec, the stored code is the stored
    # decimal coordinates normalised exactly (nerr bounds its distance from xs); otherwise xs itself.
    nerr = stored_norm_error([dec[v] for v in V], xs) if dec is not None else np.zeros_like(xs)
    log["stored_norm_err_max"] = float(nerr.max())
    dV = float(up(up(np.abs(Xf - xs) + nerr).max() + err_pts.max()))
    D2 = float(up(math.sqrt(float(up((up(up(np.abs(Xf - xs) + nerr) + err_pts) ** 2).sum()))) * (1 + 2 ** -40)))
    log.update(dist_nonrattlers=dV, dist_ls=D2)
    if rho is not None:
        log["code_inside_rho"] = bool(D2 < rho)
    # gap margin: for a hyperstatic code also tau' and the displacement 4 sqrt(n d) tau'/kappa of the local minimum
    extra_margin = Fraction(0)
    if hyper and kappa is not None:
        extra_margin = 2 * tau + Fraction(float(up(4.0 * math.sqrt(d) * float(frac_up(tau)) / kappa)))
    c_err_g = float(frac_up(Fraction(c_err) + extra_margin))
    c_lo_g = frac_down(cq - Fraction(c_err_g))
    ok_g, ginfo = S.gap_check(sysm, V, E_glob, Y, Xf, np.sqrt(d) * err_pts.max(axis=1) * S.UP_F, Z, C, K,
                              float(err.max() + dfx.max()) * S.UP_F, c_lo_g, c_err_g, d)
    log["gaps"] = ginfo
    log["gaps_ok"] = ok_g
    rat = np.setdiff1d(np.arange(n), V)
    dR = 0.0
    if len(rat):
        e_r = max(frac_up(abs(1 - sum(Fraction(float(x)) ** 2 for x in Y[r_]))) for r_ in rat)
        nr = stored_norm_error([dec[r_] for r_ in rat], X[rat]) if dec is not None else np.full((len(rat), d), 2.0 ** -52)
        dR = float(up(up(np.abs(Y[rat] - X[rat]) + nr).max() + e_r))
    log.update(dist_rattlers=dR, c_lo=repr(c_lo), c_hi=repr(c_hi))
    lo = (cq - Fraction(c_err)) * 10 ** 30
    hi = (cq + Fraction(c_err)) * 10 ** 30
    log["c_interval_30"] = [f"{lo.numerator // lo.denominator}", f"{-(-hi.numerator // hi.denominator)}"]
    ok_base = kr["ok"] and ok_span and ok_g
    if not hyper:
        if ok_base and stress_pos:
            log["status"] = "PASS"
        elif ok_base and stress_neg:
            log["status"] = "NONOPT"
        else:
            log["status"] = "FAIL " + ",".join(x for x, f in (("stress", stress_pos), ("span", ok_span),
                                                                ("gaps", ok_g)) if not f)
    else:
        # existence of a local minimum: need D(Q0) <= 4 sqrt(n d) tau'/kappa < rho
        if ok_base and stress_pos and kappa is not None:
            disp = float(up(4.0 * math.sqrt(nV * d) * float(frac_up(tau)) / kappa * (1 + 2 ** -40)))
            log["localmin_displacement"] = disp
            log["status"] = "LOCALMIN" if disp < rho else "FAIL hyperstatic: tau too large for kappa"
        else:
            log["status"] = "FAIL " + ",".join(x for x, f in (("stress", stress_pos), ("span", ok_span),
                                                                ("gaps", ok_g), ("kappa", kappa is not None)) if not f)
    if keep is not None:
        keep.update(V=V, E=Eall, Esq=Esq, extra=extra, Z=Z, C=C, K=K, R=R, sysm=sysm, xs=xs, Xf=Xf, err=err,
                    Y=Y, beta_f=beta_f, c_err=c_err, tau=tau)
    # the exact configuration (rounded to double precision, rattlers as placed) when its certified value, rounded up
    # to 12 decimals as the table displays it, is below the displayed value (a candidate improvement; checked apart)
    if dump is not None and log["status"] in ("PASS", "LOCALMIN"):
        from decimal import Decimal, ROUND_CEILING
        path, shown = dump
        top = Fraction(int(log["c_interval_30"][1]), 10 ** 30) + (tau if hyper else 0)
        top12 = Decimal(-(-top.numerator * 10 ** 12 // top.denominator)) / Decimal(10 ** 12)
        log["value_ceil12"] = str(top12)
        if top12 < Decimal(shown):
            with open(path, "w") as f:
                for p_ in Y:
                    f.write(",".join(repr(float(x)) for x in p_) + "\n")
            log["dumped"] = path
    return log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("code")
    ap.add_argument("d", type=int)
    ap.add_argument("n", type=int)
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--json")
    ap.add_argument("--single", action="store_true")
    ap.add_argument("--flex", action="store_true")
    a = ap.parse_args()
    X = read_code(a.code, a.d, a.n)
    dec = read_decimals(a.code, a.d, a.n)
    if a.flex:
        log = certify_flex(X, a.d, a.n, a.tol, dec=dec)
        print(" ".join(f"{k}={log[k]}" for k in ("status", "contacts", "points_in_contacts", "flex_beta",
                                                  "flex_radius_l1", "flex_radius_l1_stored", "forest", "cpu_seconds")
                       if k in log), flush=True)
        if a.json:
            json.dump(log, open(a.json, "w"), indent=1, default=str)
        return
    log = certify(X, a.d, a.n, a.tol, dtype=np.float32 if a.single else None, dec=dec)
    keys = ["status", "kind", "n_nonrattlers", "contacts", "isostatic_count", "unknowns", "krawczyk_beta", "r_unique",
            "stress_min_rel", "kappa", "rho", "dist_ls", "code_inside_rho", "tau_extra", "localmin_displacement",
            "cpu_seconds", "wall_seconds"]
    print(" ".join(f"{k}={log[k]}" for k in keys if k in log), flush=True)
    if a.json:
        with open(a.json, "w") as f:
            json.dump(log, f, indent=1, default=str)



# ============================================================================================ flexibility of a contact set
def certify_flex(X, d, n, tol=1e-7, log=None, dec=None, max_rows=None, qr_max_entries=12 * 10 ** 6):
    """The flexibility proposition of the paper. E = all pairs within tol of the largest cosine of X (double precision; no rattlers
    removed). L(Q) is the matrix of the norm rows (2 q_i in block i) and contact rows (q_b in block a, q_a in block b)
    of E at a configuration Q, restricted to the points that lie in a pair of E (the norm rows of the other points are
    independent of everything else). UNTRUSTED: a right inverse C of L(X) (from (L L^T)^-1, or from a QR factorisation
    of L^T when that is not accurate enough). TRUSTED: beta >= max row sum of |I - L(X) C| < 1, so L(X) C is
    nonsingular, and then every Q with |q_i - x_i|_1 < h = (1 - beta) / (2 ||C||_inf) for all i has L(Q) of full row
    rank: (Q, E) has no nonzero equilibrium stress, so no such Q (and no rotation of one) whose contact graph lies in E
    is a local minimum of the largest cosine. With dec, h_stored = h - max_i |xhat_i - x_i|_1 is the radius about the
    exactly normalised stored code (status FLEX)."""
    log = {} if log is None else log
    log["code_sha256"] = CODE_SHA256
    cpu0, t0 = time.process_time(), time.time()
    try:
        G = np.empty((n, n))
        for i0 in range(0, n, 256):                       # untrusted (only proposes E)
            G[i0:i0 + 256] = X[i0:i0 + 256] @ X.T
        np.fill_diagonal(G, -2.0)
        M = float(G.max())
        I_, J_ = np.nonzero(np.triu(G >= M - tol, 1))
        del G
        m = len(I_)
        pts = np.unique(np.concatenate([I_, J_]))         # the points that lie in a pair of E
        nE = len(pts)
        deg = np.bincount(np.concatenate([I_, J_]), minlength=n)
        log.update(d=d, n=n, tol=tol, kind="flexibility", contacts=int(m), points_in_contacts=int(nE),
                   max_contacts_at_a_point=int(deg.max()), max_cos_double=M, unknowns=int(nE + m))
        parent = list(range(n))                           # untrusted remark: is E a forest?

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        cyc = False
        for a_, b_ in zip(I_, J_):
            ra, rb = find(int(a_)), find(int(b_))
            if ra == rb:
                cyc = True
                break
            parent[ra] = rb
        log["forest"] = not cyc
        Nn, Nc = nE + m, nE * d
        if max_rows is not None and Nn > max_rows:
            log["status"] = "above size limit"
            return log
        # ---- L (exact: entries 2 x and x), rows: norm rows of pts, then contact rows; columns: the blocks of pts
        loc = -np.ones(n, int)
        loc[pts] = np.arange(nE)
        ar = np.arange(d)
        rN = np.repeat(np.arange(nE), d)
        cN = (loc[pts][:, None] * d + ar).ravel()
        vN = 2.0 * X[pts].ravel()
        rC = np.repeat(nE + np.arange(m), 2 * d)
        cC = np.concatenate([loc[I_][:, None] * d + ar, loc[J_][:, None] * d + ar], 1).ravel()
        vC = np.concatenate([X[J_], X[I_]], 1).ravel()
        L = sps.csr_matrix((np.concatenate([vN, vC]), (np.concatenate([rN, rC]), np.concatenate([cN, cC]))),
                           shape=(Nn, Nc))
        LT = L.T.tocsr()
        absL = abs(L)

        def trusted_check(Cblock):
            """TRUSTED. Cblock(j0, j1) returns the columns j0..j1-1 of C (any floats). Returns beta and ||C||_inf."""
            row_acc = np.zeros(Nn)
            crow = np.zeros(Nc)
            for j0 in range(0, Nn, BLOCK):
                Cb = Cblock(j0, min(j0 + BLOCK, Nn))
                P = np.asarray(L @ Cb)
                idx = np.arange(P.shape[1])
                Eh = -P
                Eh[j0 + idx, idx] = 1.0 - P[j0 + idx, idx]
                row_acc += np.abs(Eh).sum(1)
                crow += np.abs(Cb).sum(1)
            crow_u = up(crow)
            t1 = up(up(row_acc))
            t2 = up(gam(2 * d) * up(np.asarray(absL @ crow_u).ravel()))
            return float(np.max(up(t1 + t2))), float(np.max(crow_u))

        beta, cnorm, how = None, None, None
        Gm = (L @ LT).tocsr()                              # untrusted: L L^T in floating point
        try:
            R = approx_inverse(Gm, np.float64)
            beta, cnorm = trusted_check(lambda j0, j1: np.asarray(LT @ R[:, j0:j1]))
            how = "normal equations"
            del R
        except (np.linalg.LinAlgError, ValueError):
            pass
        log["flex_beta_normal"] = beta
        if (beta is None or not beta < 0.5) and Nn * Nc <= qr_max_entries:
            import scipy.linalg as sl                       # untrusted: C = Q R^-T from L^T = Q R (economic QR)
            Qf, Rf = sl.qr(L.T.toarray(), mode="economic", check_finite=False)
            try:
                Cfull = sl.solve_triangular(Rf, Qf.T, lower=False, check_finite=False).T   # C = Q R^-T
                del Qf, Rf
                b2, c2 = trusted_check(lambda j0, j1: Cfull[:, j0:j1])
                if beta is None or b2 < beta:
                    beta, cnorm, how = b2, c2, "QR"
            except (np.linalg.LinAlgError, ValueError):
                pass
        log.update(flex_beta=beta, flex_right_inverse=how, flex_C_norm_inf=cnorm)
        if beta is None or not beta < 1.0:
            log["status"] = "FAIL flex (row sum)"
            return log
        h = float(np.nextafter((1.0 - beta) * DEN / (2.0 * cnorm * (1 + 2.0 ** -40)), 0))
        log["flex_radius_l1"] = h                          # radius in max_i |q_i - x_i|_1
        ok = h > 0
        if dec is not None:
            ne = stored_norm_error(dec, X)
            dn = float(up(np.max(up(ne.sum(1)))))
            h_st = float(np.nextafter(h - dn, 0))
            log.update(stored_norm_dist_l1=dn, flex_radius_l1_stored=h_st)
            ok = h_st > 0
        log["status"] = "FLEX" if ok else "FAIL flex (radius)"
        return log
    finally:
        log["cpu_seconds"] = round(time.process_time() - cpu0, 2)
        log["wall_seconds"] = round(time.time() - t0, 2)


if __name__ == "__main__":
    main()
