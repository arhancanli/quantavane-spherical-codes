"""Strictness of hyperstatic contact graphs by symmetry reduction (census version 6, Phase 2).

  nice -n 10 python3 -I symstrict.py OUTDIR [--only d_n,...] [--shard I/K] [--entries open|controls] [--red KIND]

For an entry that the census certified LOCALMIN (a local minimum near the stored code; strictness open), with contacts
E on the non-rattlers V:

 1. UNTRUSTED. The census certificate is re-run (cq.certify) to obtain the zero P' of the square subsystem; P' is
    refined to 2^-400 (cq.refine, exact residuals). Permutations of V that preserve every inner product of P' (to
    1e-9) and are realised by orthogonal maps are found by backtracking (generators pi_s with numerical matrices h_s).
    An exact model of the representation is proposed: matrices M_s and a positive definite form Q over a field
    K = Q(sqrt D), from (a) simultaneous eigenspaces when every element is an involution, (b) rotation planes when
    the group is cyclic, or (c) lines that are 1-dimensional eigenspaces (eigenvalue +1 or -1) of group elements,
    with the coordinates recognised in K.
 2. TRUSTED, exact over K. Q is symmetric positive definite and every M_s satisfies M_s^T Q M_s = Q. The group is
    enumerated by closure: every product of generators that gives the same permutation of V gives the same matrix,
    so sigma -> M_sigma is a homomorphism of the permutation group Gamma into the Q-orthogonal group, and
    M_sigma p_i = p_{sigma(i)} defines the configurations with this symmetry. E is invariant under every generator.
    For one point r of each orbit, a basis F_r of the vectors fixed by the stabiliser of r is computed and checked
    exactly (M_sigma F_r = F_r for every sigma fixing r). The centraliser of the representation in the Lie algebra of
    the Q-orthogonal group is computed exactly (basis A_l).
 3. TRUSTED, ball arithmetic (Arb). The reduced system in the unknowns w (p_r = F_r w_r) and gamma: |p_r|_Q^2 = 1
    for one point r of each orbit, <p_a, p_b>_Q = gamma for one pair of each orbit of E, and one slice equation per
    A_l. A Krawczyk test (Rump 2010, Thm 13.3) proves a unique zero z* in a box. By the symmetry every point of
    P* = (M_{sigma_i} p_{r(i)}) has Q-norm 1 and every pair of E has Q-inner product gamma* exactly.
 4. TRUSTED. In an orthonormal frame (Q = L L^T, y_i = L^T p_i, enclosed in Arb): every pair of V outside E has inner
    product below gamma* (one pair per orbit of such pairs); the points span R^d; on the box around P* the Jacobian of
    the square subsystem of E' is nonsingular and the combined stress of E is enclosed and strictly positive (the
    kernel of cq.py, as in strict.py); the rattlers, placed by untrusted code, have every inner product below gamma*;
    the distance to the stored code is bounded with an exactly orthogonal rational matrix (strict.py).
Then P*, with the rattlers, satisfies the first-order criterion (paper, Theorem 2.4) and is a strict local minimum of
the largest cosine modulo rotations and rattlers, with E(P*) = E and largest cosine gamma* (status STRICT).

Red controls (--red): 'model' perturbs one entry of an exact generator matrix (the exact check must fail), 'perm'
composes one generator's permutation with a transposition inside an orbit, 'swap' replaces one contact by the closest
non-contact pair. None of them may give STRICT.
"""
import hashlib, json, math, os, sys, time, traceback, random, warnings
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, CERT)
sys.path.insert(0, os.path.join(CERT, "v5"))
warnings.filterwarnings("ignore")
import numpy as np
import mpmath
import flint
import cq
import sphcert as S
from sphcert import up
import strict as ST
import snapshot

fq = flint.fmpq
PREC = 320                    # bits of the ball arithmetic
K2 = 400                      # bits of the refined P'
GROUP_CAP = 6000              # the group is enumerated only up to this order
DCAND = [1, 5, 2, 3, 6, 15, 10, 7, 13, 30]


def _sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


SHA = dict(ST.SHA, **{"symstrict.py": _sha(os.path.abspath(__file__))})


# ============================================================================================ exact field K = Q(sqrt D)
class KM:
    """TRUSTED. A matrix A0 + A1 sqrt(D) over K = Q(sqrt D) (D squarefree, D = 1 meaning K = Q), as two FLINT rational
    matrices; every operation is exact."""
    __slots__ = ("a", "b", "D")

    def __init__(self, a, b, D):
        self.a, self.b, self.D = a, b, D

    @staticmethod
    def zero(m, n, D):
        return KM(flint.fmpq_mat(m, n), flint.fmpq_mat(m, n), D)

    @staticmethod
    def eye(n, D):
        return KM(flint.fmpq_mat(n, n, [int(i == j) for i in range(n) for j in range(n)]), flint.fmpq_mat(n, n), D)

    @staticmethod
    def from_pairs(rows, D):
        m, n = len(rows), len(rows[0])
        A = flint.fmpq_mat(m, n, [fq(x[0].numerator, x[0].denominator) for r in rows for x in r])
        B = flint.fmpq_mat(m, n, [fq(x[1].numerator, x[1].denominator) for r in rows for x in r])
        return KM(A, B, D)

    def shape(self):
        return self.a.nrows(), self.a.ncols()

    def __mul__(s, o):
        return KM(s.a * o.a + s.b * o.b * s.D, s.a * o.b + s.b * o.a, s.D)

    def __add__(s, o):
        return KM(s.a + o.a, s.b + o.b, s.D)

    def __sub__(s, o):
        return KM(s.a - o.a, s.b - o.b, s.D)

    def T(s):
        return KM(s.a.transpose(), s.b.transpose(), s.D)

    def __eq__(s, o):
        return s.a == o.a and s.b == o.b

    def key(s):
        return (str(s.a.entries()), str(s.b.entries()))

    def col(s, j):
        m = s.a.nrows()
        return KM(flint.fmpq_mat(m, 1, [s.a[i, j] for i in range(m)]), flint.fmpq_mat(m, 1, [s.b[i, j] for i in range(m)]), s.D)

    def to_arb(s):
        m, n = s.shape()
        r = flint.arb(s.D).sqrt()
        return [[flint.arb(s.a[i, j]) + flint.arb(s.b[i, j]) * r for j in range(n)] for i in range(m)]

    def to_float(s):
        m, n = s.shape()
        r = math.sqrt(s.D)
        return np.array([[float(Fraction(int(s.a[i, j].p), int(s.a[i, j].q))) +
                          r * float(Fraction(int(s.b[i, j].p), int(s.b[i, j].q))) for j in range(n)] for i in range(m)])


def hstack(mats):
    D = mats[0].D
    m = mats[0].shape()[0]
    n = sum(x.shape()[1] for x in mats)
    A, B = flint.fmpq_mat(m, n), flint.fmpq_mat(m, n)
    c = 0
    for x in mats:
        for j in range(x.shape()[1]):
            for i in range(m):
                A[i, c] = x.a[i, j]
                B[i, c] = x.b[i, j]
            c += 1
    return KM(A, B, D)


def nullspace_K(A):
    """TRUSTED (exact). A basis of {x in K^n : A x = 0}, as the columns of a KM matrix. The system is split into its
    rational and sqrt(D) parts: [[A0, D A1], [A1, A0]] (x0; x1) = 0 over Q; its rational null vectors give K-vectors
    x0 + x1 sqrt(D), of which a K-independent subset is kept (checked by the rank of the doubled matrix of the subset)."""
    m, n = A.shape()
    D = A.D
    if D == 1:
        R, rk = A.a.rref()
        piv = []
        r = 0
        for j in range(n):
            if r < rk and R[r, j] != 0:
                piv.append(j)
                r += 1
        free = [j for j in range(n) if j not in piv]
        cols = []
        for f in free:
            v = [fq(0)] * n
            v[f] = fq(1)
            for i, pj in enumerate(piv):
                v[pj] = -R[i, f]
            cols.append(v)
        if not cols:
            return KM.zero(n, 0, D)
        M = flint.fmpq_mat(n, len(cols), [cols[j][i] for i in range(n) for j in range(len(cols))])
        return KM(M, flint.fmpq_mat(n, len(cols)), D)
    big = flint.fmpq_mat(2 * m, 2 * n)
    for i in range(m):
        for j in range(n):
            big[i, j] = A.a[i, j]
            big[i, n + j] = A.b[i, j] * D
            big[m + i, j] = A.b[i, j]
            big[m + i, n + j] = A.a[i, j]
    N = nullspace_K(KM(big, flint.fmpq_mat(2 * m, 2 * n), 1)).a      # rational null vectors (2n x k)
    chosen = []
    for j in range(N.ncols()):
        trial = chosen + [j]
        # the K-span of the chosen vectors has K-dimension len(trial) iff the doubled matrix of
        # {x, sqrt(D) x : x chosen} has rational rank 2 len(trial)
        rows = []
        for t in trial:
            x0 = [N[i, t] for i in range(n)]
            x1 = [N[n + i, t] for i in range(n)]
            rows.append(x0 + x1)
            rows.append([x1[i] * D for i in range(n)] + x0)
        Mq = flint.fmpq_mat(len(rows), 2 * n, [v for r in rows for v in r])
        if Mq.rank() == 2 * len(trial):
            chosen = trial
    if not chosen:
        return KM.zero(n, 0, D)
    a = flint.fmpq_mat(n, len(chosen), [N[i, t] for i in range(n) for t in chosen])
    b = flint.fmpq_mat(n, len(chosen), [N[n + i, t] for i in range(n) for t in chosen])
    return KM(a, b, D)


def kpos_def(Q):
    """TRUSTED (exact). Q symmetric and positive definite over K: Q = Q^T and every LDL^T pivot is positive (signs of
    elements of Q(sqrt D) decided exactly)."""
    if not (Q == Q.T()):
        return False
    import kfield
    kfield.K.D = Q.D
    n = Q.shape()[0]
    A = [[kfield.K(Q.a[i, j], Q.b[i, j]) for j in range(n)] for i in range(n)]
    ok, _ = kfield.is_pos_def(A)
    return ok


# ============================================================================================ UNTRUSTED: the group
def find_auts(X, tol=1e-9, maxfound=24, seed=0, tlimit=120.0):
    """UNTRUSTED. Automorphisms of the configuration X (rows): permutations preserving every inner product to tol and
    realised by an orthogonal map; randomised backtracking on a base of d points. Returns [(perm, H)]."""
    import scipy.linalg as sl
    n, d = X.shape
    G = X @ X.T
    inv = np.stack([np.sort(G, 1)[:, -2]] + [(G ** k).sum(1) for k in (2, 3, 4)], 1)
    order = np.lexsort(inv.T[::-1])
    lab = np.zeros(n, int)
    cur = 0
    for a, b in zip(order[:-1], order[1:]):
        if np.max(np.abs(inv[a] - inv[b]) / (1 + np.abs(inv[a]))) > 1e3 * tol:
            cur += 1
        lab[b] = cur
    _, _, piv = sl.qr(X.T, pivoting=True, mode="economic")
    base = [int(i) for i in piv[:d]]
    XBinv = np.linalg.inv(X[base])
    rng = random.Random(seed)
    found, seen = [], set()
    t0 = time.time()

    def complete(img):
        H = (XBinv @ X[img]).T
        if np.max(np.abs(H.T @ H - np.eye(d))) > 1e-6:
            return None
        S_ = (X @ H.T) @ X.T
        perm = np.argmax(S_, 1)
        if np.min(S_[np.arange(n), perm]) < 1 - 1e-7 or len(set(perm.tolist())) != n:
            return None
        if np.max(np.abs(G[np.ix_(perm, perm)] - G)) > 1e3 * tol:
            return None
        return perm, H

    cand0 = [j for j in range(n) if lab[j] == lab[base[0]]]
    complete_search = False
    for _ in range(maxfound * 3):
        if len(found) >= maxfound or time.time() - t0 > tlimit:
            break
        before = len(found)
        rng.shuffle(cand0)
        ok = [False]

        def dfs(img):
            if ok[0] or time.time() - t0 > tlimit:
                return
            k = len(img)
            if k == d:
                r = complete(img)
                if r is not None and tuple(r[0].tolist()) not in seen:
                    seen.add(tuple(r[0].tolist()))
                    found.append(r)
                    ok[0] = True
                return
            b = base[k]
            cands = [j for j in range(n) if lab[j] == lab[b] and j not in img
                     and all(abs(G[img[i], j] - G[base[i], b]) < 1e3 * tol for i in range(k))]
            rng.shuffle(cands)
            for j in cands:
                dfs(img + [j])
                if ok[0]:
                    return
        for c0 in cand0:
            dfs([c0])
            if ok[0]:
                break
        if len(found) == before:
            complete_search = time.time() - t0 <= tlimit
            break
    find_auts.complete = complete_search or len(found) >= maxfound
    return found


def compose(p, q):
    """(p o q)(i) = p[q[i]]."""
    return tuple(p[i] for i in q)


def enumerate_perms(gens, cap=GROUP_CAP):
    ident = tuple(range(len(gens[0])))
    seen = {ident: []}
    frontier = [ident]
    while frontier:
        nxt = []
        for p in frontier:
            for s, g in enumerate(gens):
                q = compose(g, p)
                if q not in seen:
                    seen[q] = [s] + seen[p]       # word: q = g_{w0} o g_{w1} o ...
                    nxt.append(q)
                    if len(seen) > cap:
                        return None
        frontier = nxt
    return seen


def procrustes(P, perm):
    """UNTRUSTED. The orthogonal h minimising sum |h p_i - p_{perm(i)}|^2."""
    U, _, Vt = np.linalg.svd(P[list(perm)].T @ P)
    return U @ Vt


def recog(x, D, tol=1e-9, maxden=5000):
    """UNTRUSTED. x ~ a + b sqrt(D) with small rationals a, b; None if not found."""
    if abs(x) < 1e-11:
        return (Fraction(0), Fraction(0))
    if D == 1:
        f = Fraction(x).limit_denominator(maxden)
        return (f, Fraction(0)) if abs(float(f) - x) < tol else None
    mpmath.mp.dps = 16
    rel = mpmath.pslq([mpmath.mpf(x), 1, mpmath.sqrt(D)], tol=mpmath.mpf(tol) / 10, maxcoeff=maxden, maxsteps=400)
    if not rel or rel[0] == 0:
        return None
    a, b = Fraction(-int(rel[1]), int(rel[0])), Fraction(-int(rel[2]), int(rel[0]))
    if abs(float(a) + float(b) * math.sqrt(D) - x) > tol:
        return None
    return (a, b)


def recog_matrix(A, D):
    out = []
    for row in A:
        r = []
        for x in row:
            v = recog(float(x), D)
            if v is None:
                return None
            r.append(v)
        out.append(r)
    return out


def propose_models(gens, hs, elems_h, d, isotypic=True):
    """UNTRUSTED. Candidate exact models (B numerical frame with columns in the frame of P', M_s and Q as pairs of
    rationals over Q(sqrt D), D). Yields (how, B, Ms, Q, D)."""
    yield from _propose_direct(gens, hs, elems_h, d)
    if isotypic:
        yield from _propose_isotypic(hs, elems_h, d)


def _propose_isotypic(hs, elems_h, d):
    """UNTRUSTED. (d) The representation split into irreducible copies: the eigenspaces of a random symmetric element of
    the commutant of the generators. Copies with equal characters are mapped onto the first copy of their type by
    isometric intertwiners; each type gets an exact model on its first copy by (a)-(c) applied to the restricted
    group. The model of the whole representation is block diagonal."""
    rng = np.random.default_rng(1)
    A = np.vstack([np.kron(h.T, np.eye(d)) - np.kron(np.eye(d), h) for h in hs])
    _, sv, Vt = np.linalg.svd(A)
    null = Vt[int(np.sum(sv > 1e-8 * sv.max())):]
    if len(null) == 0:
        return
    X = sum(rng.standard_normal() * v.reshape(d, d, order="F") for v in null)
    w, U = np.linalg.eigh(X + X.T)
    clusters = []
    for i in range(d):
        if clusters and abs(w[i] - w[clusters[-1][-1]]) < 1e-7 * (1 + abs(w[i])):
            clusters[-1].append(i)
        else:
            clusters.append([i])
    if len(clusters) < 2:
        return
    copies = [U[:, c] for c in clusters]
    chars = [np.array([np.trace(W.T @ h @ W) for h in elems_h]) for W in copies]
    types = []
    for i, ch in enumerate(chars):
        for t in types:
            if copies[t[0]].shape[1] == copies[i].shape[1] and np.max(np.abs(chars[t[0]] - ch)) < 1e-6:
                t.append(i)
                break
        else:
            types.append([i])
    cols, blocksM, blocksQ = [], [], []
    Dall = 1
    for t in types:
        W1 = copies[t[0]]
        k = W1.shape[1]
        hk = [W1.T @ h @ W1 for h in hs]
        ek = [W1.T @ h @ W1 for h in elems_h]
        cand = None
        for how, Bk, Msk, Qk, Dk in _propose_direct(None, hk, ek, k):
            if Dk != 1 and Dall not in (1, Dk):
                continue
            cand = (Bk, Msk, Qk, Dk)
            break
        if cand is None:
            return
        Bk, Msk, Qk, Dk = cand
        if Dk != 1:
            Dall = Dk
        for j in t:
            Wj = copies[j]
            if j == t[0]:
                Y = np.eye(k)
            else:
                hj = [Wj.T @ h @ Wj for h in hs]
                Asy = np.vstack([np.kron(np.eye(k), a) - np.kron(b.T, np.eye(k)) for a, b in zip(hj, hk)])
                _, s2, V2 = np.linalg.svd(Asy)
                if s2[-1] > 1e-7 * max(1.0, s2[0]):
                    return
                Y = V2[-1].reshape(k, k, order="F")
                c = (Y.T @ Y)[0, 0]
                if c <= 0:
                    return
                Y = Y / math.sqrt(c)
            cols.append(Wj @ Y @ Bk)
            blocksM.append(Msk)
            blocksQ.append(Qk)
    B = np.hstack(cols)
    if B.shape != (d, d) or np.linalg.cond(B) > 1e8:
        return
    z = (Fraction(0), Fraction(0))
    Ms = []
    for s in range(len(hs)):
        M_ = [[z] * d for _ in range(d)]
        o = 0
        for Mb in blocksM:
            k = len(Mb[s])
            for i in range(k):
                for j in range(k):
                    M_[o + i][o + j] = Mb[s][i][j]
            o += k
        Ms.append(M_)
    Q = [[z] * d for _ in range(d)]
    o = 0
    for Qb in blocksQ:
        k = len(Qb)
        for i in range(k):
            for j in range(k):
                Q[o + i][o + j] = Qb[i][j]
        o += k
    yield ("isotypic", B, Ms, Q, Dall)


def _propose_direct(gens, hs, elems_h, d):
    I = np.eye(d)
    allinv = all(np.max(np.abs(h @ h - I)) < 1e-8 for h in elems_h)
    commute = all(np.max(np.abs(a @ b - b @ a)) < 1e-8 for a in hs for b in hs)
    if allinv and commute:
        # (a) elementary abelian 2-group: simultaneous eigenspaces, M_s diagonal +-1, Q = I
        cols, signs = [], []
        for eps in np.ndindex(*([2] * len(hs))):
            Pm = I.copy()
            for e, h in zip(eps, hs):
                Pm = Pm @ (I + (1 if e == 0 else -1) * h) / 2
            U, s, _ = np.linalg.svd(Pm)
            k = int(np.sum(s > 0.5))
            for j in range(k):
                cols.append(U[:, j])
                signs.append([1 if e == 0 else -1 for e in eps])
        if len(cols) == d:
            B = np.array(cols).T
            Ms = [[[(Fraction(signs[i][s] if i == j else 0), Fraction(0)) for j in range(d)] for i in range(d)]
                  for s in range(len(hs))]
            Q = [[(Fraction(int(i == j)), Fraction(0)) for j in range(d)] for i in range(d)]
            yield ("involutions", B, Ms, Q, 1)
    if len(hs) == 1 or (commute and len(elems_h) <= 60):
        # (b) cyclic: one element of maximal order generates?  Try each generator power structure through the first
        # generator only when it generates (checked by the caller through the closure).
        for gi, h in enumerate(hs):
            w, V = np.linalg.eig(h)
            ang = np.angle(w)
            m = None
            for mm in range(1, 13):
                if np.allclose(np.linalg.matrix_power(h, mm), I, atol=1e-8):
                    m = mm
                    break
            if m is None:
                continue
            cols, blocks = [], []
            for lam in (1.0, -1.0):
                _, s, Vt = np.linalg.svd(h - lam * I)
                ker = Vt[s < 1e-7].T
                for j in range(ker.shape[1]):
                    cols.append(ker[:, j])
                    blocks.append(("one", lam))
            for k in range(1, (m - 1) // 2 + 1 if m % 2 else m // 2):
                th = 2 * math.pi * k / m
                Wm = h @ h - 2 * math.cos(th) * h + I
                _, s, Vt = np.linalg.svd(Wm)
                Wb = Vt[s < 1e-7].T                     # orthonormal basis of the plane-space of angle th
                rem = Wb.copy()
                while rem.shape[1] > 0:
                    v = rem[:, 0] / np.linalg.norm(rem[:, 0])
                    hv = h @ v
                    cols += [v, hv]
                    blocks.append(("rot", th, k, m))
                    Pp = np.array([v, hv]).T
                    proj = rem - Pp @ np.linalg.lstsq(Pp, rem, rcond=None)[0]
                    U, s2, _ = np.linalg.svd(proj, full_matrices=False)
                    rem = U[:, s2 > 1e-6]
            if len(cols) != d:
                continue
            B = np.array(cols).T
            # exact blocks: cos(2 pi k / m) in Q or a quadratic field
            ckm = {}
            D = 1
            ok = True
            for bl in blocks:
                if bl[0] == "rot":
                    c = math.cos(bl[1])
                    for DD in (1, 5, 2, 3):
                        v = recog(c, DD)
                        if v is not None:
                            if v[1] != 0:
                                if D not in (1, DD):
                                    ok = False
                                D = DD
                            ckm[bl] = v
                            break
                    else:
                        ok = False
            if not ok:
                continue
            Mg = [[(Fraction(0), Fraction(0))] * d for _ in range(d)]
            Q = [[(Fraction(0), Fraction(0))] * d for _ in range(d)]
            c_ = 0
            for bl in blocks:
                if bl[0] == "one":
                    Mg[c_][c_] = (Fraction(int(bl[1])), Fraction(0))
                    Q[c_][c_] = (Fraction(1), Fraction(0))
                    c_ += 1
                else:
                    a, b = ckm[bl]
                    # basis (v, h v): h v = second vector; h (h v) = 2 cos th h v - v
                    Mg[c_][c_] = (Fraction(0), Fraction(0))
                    Mg[c_ + 1][c_] = (Fraction(1), Fraction(0))
                    Mg[c_][c_ + 1] = (Fraction(-1), Fraction(0))
                    Mg[c_ + 1][c_ + 1] = (2 * a, 2 * b)
                    Q[c_][c_] = Q[c_ + 1][c_ + 1] = (Fraction(1), Fraction(0))
                    Q[c_][c_ + 1] = Q[c_ + 1][c_] = (a, b)
                    c_ += 2
            # the other generators in this basis must be recognised too (powers of the same element)
            Binv = np.linalg.inv(B)
            Ms = []
            good = True
            for hj in hs:
                if hj is h:
                    Ms.append(Mg)
                    continue
                Mj = recog_matrix(Binv @ hj @ B, D)
                if Mj is None:
                    good = False
                    break
                Ms.append(Mj)
            if good:
                yield ("cyclic", B, Ms, Q, D)
    # (c) lines: 1-dimensional eigenspaces of group elements for the eigenvalues +1 and -1
    lines = []
    for h in elems_h:
        for lam in (1.0, -1.0):
            _, s, Vt = np.linalg.svd(h - lam * I)
            if np.sum(s < 1e-7) == 1:
                u = Vt[-1]
                u = u * (1 if u[np.argmax(np.abs(u))] > 0 else -1)
                if all(abs(abs(u @ v) - 1) > 1e-9 for v in lines):
                    lines.append(u)
    if not lines:
        return
    L = np.array(lines)
    # orbits of lines under the generators
    par = list(range(len(lines)))

    def f(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a
    for h in hs:
        img = (L @ h.T) @ L.T
        for i in range(len(lines)):
            j = int(np.argmax(np.abs(img[i])))
            if abs(abs(img[i, j]) - 1) < 1e-8:
                a, b = f(i), f(j)
                if a != b:
                    par[a] = b
    orbs = {}
    for i in range(len(lines)):
        orbs.setdefault(f(i), []).append(i)
    orbs = sorted(orbs.values(), key=len)
    import scipy.linalg as sl
    tries = []
    for o in orbs:
        if np.linalg.matrix_rank(L[o], tol=1e-8) == d:
            tries.append(o)
    # greedy unions (orbits in order of size)
    acc = []
    for o in orbs:
        acc = acc + o
        if np.linalg.matrix_rank(L[acc], tol=1e-8) == d:
            tries.append(list(acc))
            break
    for o in tries:
        _, _, piv = sl.qr(L[o].T, pivoting=True)
        B = L[[o[i] for i in piv[:d]]].T
        if np.linalg.cond(B) > 1e6:
            continue
        Binv = np.linalg.inv(B)
        Qn = B.T @ B
        for D in DCAND:
            Qr = recog_matrix(Qn, D)
            if Qr is None:
                continue
            Ms = []
            for h in hs:
                Mj = recog_matrix(Binv @ h @ B, D)
                if Mj is None:
                    break
                Ms.append(Mj)
            if len(Ms) == len(hs):
                yield ("lines", B, Ms, Qr, D)
                break


# ============================================================================================ the certificate
def sym_entry(X, dec, d, n, tol, log, red=None):
    t0 = time.time()
    keep = {}
    cl = cq.certify(X, d, n, tol, dec=dec, keep=keep)
    log["census_status_rerun"] = cl.get("status")
    if cl.get("status") != "LOCALMIN":
        log["status"] = "not LOCALMIN on re-run"
        return log
    V, E, Esq, extra = keep["V"], list(keep["E"]), keep["Esq"], keep["extra"]
    nV, m = len(V), len(E)
    log.update(nV=int(nV), contacts=int(m), extra=len(extra), rattlers=int(n - nV))
    # ---- 1. refine P' (untrusted) and its symmetry
    sysm, R = keep["sysm"], keep["R"]
    sh = K2 - keep["K"]
    Z = [[z << sh for z in row] for row in keep["Z"]]
    C = keep["C"] << sh
    ex = S.Exact(sysm, K2)
    Z, C, res = cq.refine(ex, Z, C, R, K2, iters=8)
    del keep["R"]
    P = np.array([[z / (1 << K2) for z in row] for row in Z])
    P = P / np.linalg.norm(P, axis=1, keepdims=True)
    found = find_auts(P)
    gens_all = [tuple(int(x) for x in f[0]) for f in found if any(int(f[0][i]) != i for i in range(nV))]
    log["automorphisms_found"] = len(gens_all)
    gens, order = [], 1                       # a small generating set: keep a candidate only if the group grows
    for g in gens_all:
        cl = enumerate_perms(gens + [g])
        if cl is None:
            gens = gens + [g]
            order = None
            break
        if len(cl) > order:
            gens, order = gens + [g], len(cl)
    log["generators_found"] = len(gens)
    log["automorphism_search_complete"] = bool(find_auts.complete)
    if not gens:
        log["status"] = "no symmetry" if find_auts.complete else "no symmetry found within the time limit"
        return log
    if red == "swap":                                    # red control: one contact replaced by the next pair
        G = P @ P.T
        Eset = set(E)
        best = max(((G[a, b], a, b) for a in range(nV) for b in range(a + 1, nV) if (a, b) not in Eset))
        E[0] = (best[1], best[2])
        log["red_swap"] = [list(map(int, E[0]))]
    perms = enumerate_perms(gens)
    if perms is None:
        log["status"] = "group above the enumeration cap"
        return log
    log["group_order"] = len(perms)
    hs = [procrustes(P, g) for g in gens]
    elems_h = []
    for p, word in perms.items():
        h = np.eye(d)
        for s in reversed(word):
            h = hs[s] @ h
        elems_h.append(h)
    if red == "perm":                                    # red control: a wrong permutation for one generator
        g = list(gens[0])
        moved = [i for i in range(nV) if g[i] != i]
        i_, j_ = moved[0], moved[1]
        g[i_], g[j_] = g[j_], g[i_]
        gens[0] = tuple(g)
        log["red_perm"] = [int(i_), int(j_)]
    # ---- 2. exact model (trusted checks)
    model = None
    tried = []
    for how, B, Ms_r, Q_r, D in propose_models(gens, hs, elems_h, d):
        Q = KM.from_pairs(Q_r, D)
        Ms = [KM.from_pairs(M_, D) for M_ in Ms_r]
        if red == "model":
            Ms[0].a[0, 0] = Ms[0].a[0, 0] + fq(1, 997)
        okQ = kpos_def(Q)
        okO = all(M_.T() * Q * M_ == Q for M_ in Ms)
        tried.append([how, D, bool(okQ), bool(okO)])
        if okQ and okO:
            model = (how, B, Ms, Q, D)
            break
    log["models_tried"] = tried
    if model is None:
        log["status"] = "no exact model"
        return log
    how, B, Ms, Q, D = model
    log.update(model=how, field_D=D)
    # closure: sigma -> M_sigma well defined (a homomorphism), exactly
    Mof = {tuple(range(nV)): KM.eye(d, D)}
    frontier = [tuple(range(nV))]
    while frontier:
        nxt = []
        for p in frontier:
            for s, g in enumerate(gens):
                q = compose(g, p)
                Mq = Ms[s] * Mof[p]
                if q in Mof:
                    if not (Mof[q] == Mq):
                        log["status"] = "exact model is not a representation of the permutation group"
                        return log
                else:
                    Mof[q] = Mq
                    nxt.append(q)
                    if len(Mof) > GROUP_CAP:
                        log["status"] = "group above the enumeration cap"
                        return log
        frontier = nxt
    log["group_order_exact"] = len(Mof)
    Eset = set(tuple(sorted(e)) for e in E)
    if not all(tuple(sorted((g[a], g[b]))) in Eset for g in gens for (a, b) in Eset):
        log["status"] = "E not invariant"
        return log
    # orbits, transversals, stabilisers
    rep = [-1] * nV
    trans = [None] * nV
    for r in range(nV):
        if rep[r] >= 0:
            continue
        for p in Mof:
            i = p[r]
            if rep[i] < 0:
                rep[i] = r
                trans[i] = p
    reps = sorted(set(rep))
    Fb = {}
    for r in reps:
        stab = [p for p in Mof if p[r] == r]
        rows = [Mof[p] - KM.eye(d, D) for p in stab if p != tuple(range(nV))]
        if rows:
            A = KM(flint.fmpq_mat(len(rows) * d, d, [x for M_ in rows for x in M_.a.entries()]),
                   flint.fmpq_mat(len(rows) * d, d, [x for M_ in rows for x in M_.b.entries()]), D)
            F = nullspace_K(A)
        else:
            F = KM.eye(d, D)
        if F.shape()[1] == 0:
            log["status"] = "a stabiliser fixes no nonzero vector"
            return log
        for p in stab:                                    # exact check M_sigma F = F
            if not (Mof[p] * F == F):
                log["status"] = "fixed-vector basis check failed"
                return log
        Fb[r] = F
    # edge orbits (union-find under the generators), one pair per orbit
    el = sorted(Eset)
    eid = {e: k for k, e in enumerate(el)}
    par = list(range(len(el)))

    def fnd(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a
    for g in gens:
        for k, (a, b) in enumerate(el):
            j = eid[tuple(sorted((g[a], g[b])))]
            x, y = fnd(k), fnd(j)
            if x != y:
                par[x] = y
    erep = {}
    for k, e in enumerate(el):
        erep.setdefault(fnd(k), e)
    erep = list(erep.values())
    # centraliser of the representation in the Lie algebra of the Q-orthogonal group: A^T Q + Q A = 0, A M_s = M_s A
    nA = d * d
    eqs_a, eqs_b = [], []
    for i in range(d):
        for j in range(i, d):
            ra, rb = [fq(0)] * nA, [fq(0)] * nA
            for k in range(d):                            # (A^T Q)_{ij} + (Q A)_{ij} = sum_k A_{ki} Q_{kj} + Q_{ik} A_{kj}
                ra[k * d + i] += Q.a[k, j]
                rb[k * d + i] += Q.b[k, j]
                ra[k * d + j] += Q.a[i, k]
                rb[k * d + j] += Q.b[i, k]
            eqs_a.append(ra)
            eqs_b.append(rb)
    for M_ in Ms:
        for i in range(d):
            for j in range(d):                            # (A M - M A)_{ij} = sum_k A_{ik} M_{kj} - M_{ik} A_{kj}
                ra, rb = [fq(0)] * nA, [fq(0)] * nA
                for k in range(d):
                    ra[i * d + k] += M_.a[k, j]
                    rb[i * d + k] += M_.b[k, j]
                    ra[k * d + j] -= M_.a[i, k]
                    rb[k * d + j] -= M_.b[i, k]
                eqs_a.append(ra)
                eqs_b.append(rb)
    Aeq = KM(flint.fmpq_mat(len(eqs_a), nA, [x for r in eqs_a for x in r]),
             flint.fmpq_mat(len(eqs_b), nA, [x for r in eqs_b for x in r]), D)
    Cn = nullspace_K(Aeq)
    Al = []
    for j in range(Cn.shape()[1]):
        v = Cn.col(j)
        Al.append(KM(flint.fmpq_mat(d, d, [v.a[t, 0] for t in range(nA)]),
                     flint.fmpq_mat(d, d, [v.b[t, 0] for t in range(nA)]), D))
    for A_ in Al:                                         # exact check
        if not (A_.T() * Q + Q * A_ == KM.zero(d, d, D)) or not all(A_ * M_ == M_ * A_ for M_ in Ms):
            log["status"] = "centraliser check failed"
            return log
    kdim = {r: Fb[r].shape()[1] for r in reps}
    nunk = sum(kdim.values()) + 1
    neq = len(reps) + len(erep) + len(Al)
    log.update(point_orbits=len(reps), edge_orbits=len(erep), dim_fix=nunk - 1, dim_centraliser=len(Al),
               reduced_unknowns=nunk, reduced_equations=neq)
    if neq != nunk:
        log["status"] = "reduced system not square (%d equations, %d unknowns)" % (neq, nunk)
        return log
    # ---- 3. reduced system and its Krawczyk test (Arb)
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    try:
        out = reduced_krawczyk(P, B, Ms, Q, D, Mof, Fb, reps, rep, trans, erep, Al, nV, d, log)
    finally:
        flint.ctx.prec = old
    if out is None:
        return log
    Yb, gam = out
    # ---- 4. the full framework at P* (orthonormal frame)
    return full_checks(X, dec, V, E, Yb, gam, Mof, nV, n, d, log, t0)


def kvec_arb(F):
    return F.to_arb()


def reduced_krawczyk(P, B, Ms, Q, D, Mof, Fb, reps, rep, trans, erep, Al, nV, d, log):
    """TRUSTED (Arb), except the starting point and the approximate inverse."""
    Qa = Q.to_arb()
    # coefficient matrices (exact in K, evaluated as balls)
    off = {}
    o = 0
    for r in reps:
        off[r] = o
        o += Fb[r].shape()[1]
    nw = o
    u = nw + 1
    QF = {r: (Q * Fb[r]) for r in reps}
    Nr = {r: (Fb[r].T() * Q * Fb[r]).to_arb() for r in reps}
    Cab = []
    for (a, b) in erep:
        ra, rb = rep[a], rep[b]
        Mat = (Mof[trans[a]] * Fb[ra]).T() * Q * (Mof[trans[b]] * Fb[rb])
        Cab.append((ra, rb, Mat.to_arb()))
    # starting point: model coordinates of P' (untrusted), least squares in each Fix(Stab(r))
    Pm = np.linalg.solve(B, P.T).T                         # p'_i = B x_i
    w0 = np.zeros(u)
    for r in reps:
        Ff = Fb[r].to_float()
        w0[off[r]:off[r] + Ff.shape[1]] = np.linalg.lstsq(Ff, Pm[r], rcond=None)[0]
    G = P @ P.T
    w0[-1] = float(np.mean([G[a, b] for a, b in erep]))
    # slice rows (exact reference xi_r: the start rounded to 2^-30)
    Srows = []
    for A_ in Al:
        row = [flint.arb(0)] * u
        for r in reps:
            Ff = Fb[r].to_float()
            xi = np.round((Ff @ w0[off[r]:off[r] + Ff.shape[1]]) * 2 ** 30) / 2 ** 30
            xiK = KM(flint.fmpq_mat(1, d, [fq(int(x * 2 ** 30), 2 ** 30) for x in xi]), flint.fmpq_mat(1, d), D)
            coef = (xiK * Q * A_ * Fb[r]).to_arb()[0]
            for j, cval in enumerate(coef):
                row[off[r] + j] = cval
        Srows.append(row)

    def fval(z):
        """f at a ball vector z (list of arb)."""
        out = []
        for r in reps:
            w = z[off[r]:off[r] + len(Nr[r])]
            acc = flint.arb(0)
            for i in range(len(w)):
                for j in range(len(w)):
                    acc += w[i] * Nr[r][i][j] * w[j]
            out.append(acc - 1)
        for ra, rb, Mt in Cab:
            wa = z[off[ra]:off[ra] + len(Mt)]
            wb = z[off[rb]:off[rb] + len(Mt[0])]
            acc = flint.arb(0)
            for i in range(len(wa)):
                for j in range(len(wb)):
                    acc += wa[i] * Mt[i][j] * wb[j]
            out.append(acc - z[-1])
        for row in Srows:
            acc = flint.arb(0)
            for j in range(nw):
                acc += row[j] * z[j]
            out.append(acc)
        return out

    def jac(z):
        J = [[flint.arb(0)] * u for _ in range(u)]
        k = 0
        for r in reps:
            w = z[off[r]:off[r] + len(Nr[r])]
            for i in range(len(w)):
                acc = flint.arb(0)
                for j in range(len(w)):
                    acc += (Nr[r][i][j] + Nr[r][j][i]) * w[j]
                J[k][off[r] + i] = acc
            k += 1
        for ra, rb, Mt in Cab:
            wa = z[off[ra]:off[ra] + len(Mt)]
            wb = z[off[rb]:off[rb] + len(Mt[0])]
            row = [flint.arb(0)] * u
            for i in range(len(wa)):
                acc = flint.arb(0)
                for j in range(len(wb)):
                    acc += Mt[i][j] * wb[j]
                row[off[ra] + i] = row[off[ra] + i] + acc
            for j in range(len(wb)):
                acc = flint.arb(0)
                for i in range(len(wa)):
                    acc += wa[i] * Mt[i][j]
                row[off[rb] + j] = row[off[rb] + j] + acc
            row[u - 1] = flint.arb(-1)
            J[k] = row
            k += 1
        for row in Srows:
            J[k] = row[:nw] + [flint.arb(0)]
            k += 1
        return J

    def to_np(J):
        return np.array([[float(x.mid()) for x in row] for row in J])

    # Newton (chord) steps at high precision, untrusted
    z = [flint.arb(float(x)) for x in w0]
    Cn = np.linalg.inv(to_np(jac(z)))
    for it in range(12):
        F = fval(z)
        Fm = np.array([float(x.mid()) for x in F])
        if it > 0 and np.max(np.abs(Fm)) < 2.0 ** -(PREC - 40):
            break
        # step: z - C F, F at full precision (C applied to the ball midpoints in arb)
        Ca = [[flint.arb(float(c)) for c in row] for row in Cn]
        step = [sum((Ca[i][j] * F[j].mid() for j in range(u)), flint.arb(0)) for i in range(u)]
        z = [flint.arb((z[i] - step[i]).mid()) for i in range(u)]
        if it == 0:
            Cn = np.linalg.inv(to_np(jac(z)))
    zt = [flint.arb(x.mid()) for x in z]                    # exact binary midpoints
    Ft = fval(zt)
    log["reduced_residual"] = float(max(abs(x).upper() for x in Ft))
    Cn = np.linalg.inv(to_np(jac(zt)))
    Ca = flint.arb_mat([[flint.arb(float(c)) for c in row] for row in Cn])   # exact floats
    CF = Ca * flint.arb_mat([[x] for x in Ft])
    ybar = max(float(abs(CF[i, 0]).upper()) for i in range(u))
    ok = False
    attempts = []
    for fac in (4.0, 64.0, 2.0 ** 20):
        r = max(ybar * fac, 2.0 ** -(PREC - 30))
        box = [zt[i] + flint.arb(0, r) for i in range(u)]
        Jb = flint.arb_mat(jac(box))
        Imat = flint.arb_mat([[flint.arb(int(i == j)) for j in range(u)] for i in range(u)])
        Km = Imat - Ca * Jb
        Xr = flint.arb_mat([[flint.arb(0, r)] for _ in range(u)])
        Kv = -CF + Km * Xr
        mx = max(float(abs(Kv[i, 0]).upper()) for i in range(u))
        attempts.append({"r": r, "max_over_r": mx / r})
        if mx < r:
            ok = True
            break
    log["krawczyk_attempts"] = attempts
    if not ok:
        log["status"] = "FAIL existence (reduced Krawczyk)"
        return None
    zs = [zt[i] + Kv[i, 0] for i in range(u)]              # the zero lies in z~ + K
    gam = zs[-1]
    log["gamma_mid"] = gam.mid().str(60, radius=False)
    log["gamma_rad"] = float(gam.rad())
    log["reduced_box_radius"] = float(max(x.rad() for x in zs))
    # P* in the model, then an orthonormal frame y = L^T p (Q = L L^T, Cholesky in Arb)
    pr = {}
    for r in reps:
        F = Fb[r].to_arb()
        w = zs[off[r]:off[r] + len(F[0])]
        pr[r] = [sum((F[i][j] * w[j] for j in range(len(w))), flint.arb(0)) for i in range(d)]
    Lc = [[flint.arb(0)] * d for _ in range(d)]
    for i in range(d):
        for j in range(i + 1):
            s_ = Qa[i][j]
            for k in range(j):
                s_ = s_ - Lc[i][k] * Lc[j][k]
            if i == j:
                if not s_ > 0:
                    log["status"] = "Q not positive in ball Cholesky"
                    return None
                Lc[i][i] = s_.sqrt()
            else:
                Lc[i][j] = s_ / Lc[j][j]
    Yb = []
    for i in range(nV):
        Mi = Mof[trans[i]].to_arb()
        p = [sum((Mi[a][b] * pr[rep[i]][b] for b in range(d)), flint.arb(0)) for a in range(d)]
        y = [sum((Lc[k][a] * p[k] for k in range(a, d)), flint.arb(0)) for a in range(d)]   # (L^T p)_a
        Yb.append(y)
    return Yb, gam


def full_checks(X, dec, V, E, Yb, gam, Mof, nV, n, d, log, t0):
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    Xf = np.empty((nV, d))
    rad = np.empty((nV, d))
    for i in range(nV):
        for k in range(d):
            f = float(Yb[i][k].mid())
            Xf[i, k] = f
            rad[i, k] = math.nextafter(float(abs(Yb[i][k] - flint.arb(f)).upper()), math.inf)
    rad = up(rad)
    log["enclosure_radius"] = float(rad.max())
    # every pair of V outside E below gamma*: one pair per orbit of pairs (inner products are constant on orbits)
    gens_p = list(Mof.keys())
    Eset = set(tuple(sorted(e)) for e in E)
    seen = set()
    worst = None
    for a in range(nV):
        for b in range(a + 1, nV):
            if (a, b) in Eset or (a, b) in seen:
                continue
            for p in gens_p:
                seen.add(tuple(sorted((p[a], p[b]))))
            ip = sum((Yb[a][k] * Yb[b][k] for k in range(d)), flint.arb(0))
            diff = gam - ip
            if not diff > 0:
                flint.ctx.prec = old
                log["status"] = "a pair outside E not provably below gamma*"
                log["bad_pair"] = [a, b]
                return log
            lo = float(diff.lower())
            worst = lo if worst is None else min(worst, lo)
    log["gap_lower"] = worst
    flint.ctx.prec = old
    # stress on the box (kernel of cq.py), as strict.py stage 3
    iso = nV * (d - 1) - d * (d - 1) // 2 + 1
    Elist = sorted(Eset)
    m = len(Elist)
    full = S.System(Xf, Elist, d)
    sub = cq.choose_square_subset_qr(full, Xf, d, nV, iso)
    Esq = [Elist[e] for e in sub]
    ss = set(sub)
    extra = [Elist[e] for e in range(m) if e not in ss]
    sysm = S.System(Xf, Esq, d)
    N = sysm.nu
    Jc = sysm.jac(Xf)
    try:
        R = cq.approx_inverse(Jc, np.float64)
    except (np.linalg.LinAlgError, ValueError):
        log["status"] = "square Jacobian numerically singular at P*"
        return log
    beta_f = cq.choose_beta(Esq, extra, Xf, R, N, nV, m, d, log)
    if beta_f is None:
        log["status"] = "no positive stress found at P*"
        return log
    Jr = sysm.jac_radius(rad)
    b, brad = cq.extra_rhs(extra, beta_f, Xf, rad, N, d)
    yt, ybound = cq.transposed_enclosure(sysm, Jc, Jr, R, b, brad, log, "stress")
    del R
    if yt is None:
        log["status"] = "stress enclosure failed at P*"
        return log
    om_lo = np.nextafter(yt[nV:nV + len(Esq)] - ybound[nV:nV + len(Esq)], -np.inf)
    wmin = float(min(om_lo.min(), beta_f.min()))
    log.update(stress_min_lower=wmin, stress_min_rel=wmin * m)
    if not (om_lo.min() > 0 and beta_f.min() > 0):
        log["status"] = "stress not provably positive at P*"
        return log
    ok_span, span_rs, _ = S.span_check(Xf, rad, d)
    log["span_rowsum"] = span_rs
    if not ok_span:
        log["status"] = "span failed at P*"
        return log
    # rattlers in the frame of P*
    rat = np.setdiff1d(np.arange(n), V)
    U_, _, Vt_ = np.linalg.svd(X[V].T @ Xf)
    O = U_ @ Vt_
    Ycfg = np.empty((n, d))
    Ycfg[V] = Xf
    Ycfg[rat] = X[rat] @ O
    c_lo = Fraction(math.nextafter(float(gam.lower()), -math.inf))
    log["c_star_lower"] = float(c_lo)
    if len(rat):
        Yr, rinfo = S.place_rattlers(Ycfg, V, Xf, float(c_lo), d)
        log["rattler_placement"] = rinfo
        ok_r, rbad = ST.rattler_gaps(Yr, V, rat, Xf, rad, c_lo, d)
        log["rattler_gaps_ok"] = ok_r
        if not ok_r:
            log["status"] = "rattler gaps failed"
            return log
    # distance from the stored code (exactly orthogonal rational rotation, ball arithmetic)
    dist = ST.distance_to_stored(dec, V, Yb, O, d)
    log.update(dist)
    if not dist.get("rotation_exactly_orthogonal") or not math.isfinite(dist.get("dist_ls_upper", float("nan"))):
        log["status"] = "distance not certified"
        return log
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    log["gamma_ball"] = gam.str(50, radius=True)
    flint.ctx.prec = old
    log["status"] = "STRICT"
    log["seconds"] = round(time.time() - t0, 2)
    return log


def main():
    out = sys.argv[1]
    args = sys.argv[2:]
    sh, nsh = 0, 1
    only = None
    red = None
    which = "open"
    if "--shard" in args:
        sh, nsh = (int(x) for x in args[args.index("--shard") + 1].split("/"))
    if "--only" in args:
        only = set(args[args.index("--only") + 1].split(","))
    if "--red" in args:
        red = args[args.index("--red") + 1]
    if "--entries" in args:
        which = args[args.index("--entries") + 1]
    os.makedirs(out, exist_ok=True)
    s5 = json.load(open(os.path.join(CERT, "v5", "out", "summary_v5.json")))
    want = {"open": ("local",), "controls": ("strict",)}[which]
    keys = sorted(tuple(map(int, k.split(","))) for k, v in s5["class"].items()
                  if v[0] in want and (which == "open" or v[1] == "strict_hyper"))
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    todo = []
    for d, n in keys:
        if only and f"{d}_{n}" not in only:
            continue
        r = json.load(open(os.path.join(CERT, "census", "results_v4", f"d{d}_n{n}.json")))
        a = [x for x in r["attempts"] if x.get("status") == "LOCALMIN"][0]
        todo.append((a["unknowns"], d, n, a["tol"]))
    todo.sort()
    for k, (size, d, n, tol) in enumerate(todo):
        if k % nsh != sh:
            continue
        js = os.path.join(out, f"d{d}_n{n}.json")
        if os.path.exists(js) and "--force" not in args:
            continue
        log = {"d": d, "n": n, "tol": tol, "unknowns": size, "code_sha256": SHA, "red_control": red}
        t0, c0 = time.time(), time.process_time()
        try:
            pts = snapshot.read_points(SNAP, idx[(d, n)])
            dec = [[t.strip() for t in p] for p in pts]
            X = np.array([[float(x) for x in p] for p in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            sym_entry(X, dec, d, n, tol, log, red=red)
        except Exception:
            log["status"] = "error"
            log["error"] = traceback.format_exc()[-2000:]
        log["wall"] = round(time.time() - t0, 2)
        log["cpu_seconds"] = round(time.process_time() - c0, 2)
        with open(js + ".part", "w") as fo:
            json.dump(log, fo, default=str)
        os.replace(js + ".part", js)
        print(f"d{d}_n{n} {size} {log['status']} {log['wall']}s", flush=True)


if __name__ == "__main__":
    main()
