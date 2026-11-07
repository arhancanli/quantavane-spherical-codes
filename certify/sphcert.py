"""Computer-assisted certificate that a spherical code lies next to an exact strict local optimum.

  python3 -I certify/sphcert.py CODE.txt D N [--tol 1e-10] [--bits 100] [--json OUT.json] [--duty 1]

CODE.txt holds N points of R^D, one per line, coordinates separated by commas (the format of the codes/ folder; other
lines are ignored). PASS means that the following was proved, for an exact configuration p* of N unit vectors whose
distance to the code is printed (the paper, Sections 2 and 3, gives the statements and proofs):

  (a) existence: the points of V (the code's points minus its rattlers) are within the printed radius of a solution of
      |p_i|^2 = 1 (i in V), <p_i, p_j> = c (ij in E, the contacts), sum_i <A_l x_i, p_i> = 0 (l = 1..D(D-1)/2: a slice
      that removes rotations; x is the code, A_l the elementary skew-symmetric matrices), with c in the printed interval;
      the solution is unique in that box and the Jacobian there is non-singular (Rump, Acta Numerica 2010, Theorem 13.3);
  (b) span: the points of V at p* span R^D;
  (c) stress: the unique equilibrium stress of the contacts (normalised to total 1) is strictly positive on every
      contact;
  (d) gaps: every pair of points outside E, rattlers included, has inner product strictly below c at p*.
By the classical criterion recalled in the paper (Section 2), p* is then a strict local minimum of the largest inner
product, modulo rotations and motions of the rattlers, and its largest inner product is c.

The program has two kinds of functions. UNTRUSTED ones (choice of contacts and rattlers, Newton refinement, approximate
inverse, rattler placement) only propose candidates; nothing they return is assumed to be correct. TRUSTED ones
(marked so) check the candidates, with exact integer arithmetic or with floating-point arithmetic and rigorous error
bounds (Lemma R and Lemma S below). Only IEEE double precision with rounding to nearest is assumed; no dense
matrix-matrix product of a BLAS library is used in the trusted part (only scipy's sparse products, elementwise numpy
operations and numpy sums), so the error bounds hold whatever the order of summation.

Lemma R (sums of nonnegative terms). Let s be the computed value of a sum of k <= 2^20 terms, each a nonnegative
float or a product of two nonnegative floats, in any order, with gradual underflow. Then the exact sum is at most
up(s) = fl(fl(s * (1 + 2^-29)) + 2^-996).  Proof: with u = 2^-53 and eta = 2^-1074, every product is rounded to at
least p(1 - u) - eta/2 and every addition of nonnegative numbers to at least (a + b)(1 - u); so s >= (1 - k u) S - k eta/2
for the exact sum S, and the two operations in up() give at least (s (1 + 2^-29)(1 - u) - eta/2 + 2^-996)(1 - u),
which is at least (s + k eta/2) / (1 - k u) >= S because (1 - u)^2 (1 + 2^-29) >= 1/(1 - 2^-33) and k <= 2^20.
Lemma S (signed dot products; Rump 2010, Theorem 2.1, with underflow added). The computed value of sum_i a_i b_i
(k <= 2^20 terms, any order) differs from the exact one by at most gamma_k sum_i |a_i b_i| + k eta, gamma_k = k u /
(1 - k u) <= 2 k u; hence by at most up(2 k u * up(fl(sum_i |a_i| |b_i|))).
"""
import argparse, json, math, os, sys, time
from fractions import Fraction

import numpy as np
import scipy.sparse as sps

U = 2.0 ** -53            # unit roundoff, binary64, rounding to nearest
ETA = 2.0 ** -1074        # smallest positive subnormal
UP_F = 1.0 + 2.0 ** -29
UP_T = 2.0 ** -996
KMAX = 2 ** 20
BLOCK = 128               # rows or columns per block in the large products (memory BLOCK x size)


class Pace:
    """Keeps the process a light background job on a loaded machine: after every `slice_` seconds of CPU time it sleeps
    so that its CPU time is at most the fraction `duty` of the wall time (duty 1: no sleeping). Only timing changes."""

    def __init__(self, duty=1.0, slice_=0.04):
        self.duty, self.slice = duty, slice_
        self.cpu = time.process_time()

    def tick(self):
        if self.duty >= 1.0:
            return
        used = time.process_time() - self.cpu
        if used >= self.slice:
            time.sleep(used * (1.0 / self.duty - 1.0))
            self.cpu = time.process_time()


PACE = Pace()


# ============================================================================================ TRUSTED: arithmetic
def up(x):
    """TRUSTED. Upper bound of the exact value of a nonnegative computed sum (Lemma R)."""
    return x * UP_F + UP_T


def gam(k):
    """TRUSTED. A float >= gamma_k = k u / (1 - k u), for 1 <= k <= 2^20 (exact: 2 k u)."""
    assert 1 <= k <= KMAX
    return 2.0 * k * U


def frac_up(q):
    """TRUSTED. The least float >= the rational q (or the next one)."""
    f = float(q)
    if Fraction(f) < q:
        f = math.nextafter(f, math.inf)
    return f


def frac_down(q):
    """TRUSTED. A float <= the rational q."""
    f = float(q)
    if Fraction(f) > q:
        f = math.nextafter(f, -math.inf)
    return f


# ============================================================================================ reading
def read_code(path, d, n):
    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
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
    return X


# ============================================================================================ UNTRUSTED: structure
def contact_structure(X, d, tol):
    """UNTRUSTED. Contacts: pairs within tol of the largest cosine (double precision). Rattlers: points with fewer than d
    contacts among the remaining points, removed repeatedly. Returns V (sorted indices), E (pairs of positions in V),
    the largest cosine and the Gram matrix."""
    Xn = X / np.linalg.norm(X, axis=1, keepdims=True)
    G = Xn @ Xn.T
    np.fill_diagonal(G, -2.0)
    M = G.max()
    A = G >= M - tol
    alive = np.ones(len(X), bool)
    while True:
        deg = (A & alive[None, :] & alive[:, None]).sum(1)
        drop = alive & (deg < d)
        if not drop.any():
            break
        alive &= ~drop
    V = np.nonzero(alive)[0]
    pos = -np.ones(len(X), int)
    pos[V] = np.arange(len(V))
    I, J = np.nonzero(np.triu(A & alive[None, :] & alive[:, None], 1))
    E = [(int(pos[i]), int(pos[j])) for i, j in zip(I, J)]
    return V, E, float(M), G


class System:
    """The square system F(z) = 0, z = (p_1, ..., p_n, c) with p_a in R^d (the points of V in order), equations
      S_a = |p_a|^2 - 1                     a = 0..n-1
      C_e = <p_a, p_b> - c                  e = (a, b) in E
      G_l = sum_a <A_l xs_a, p_a>           l = (s, t), s < t  (A_l = E_st - E_ts; xs = the code's points of V)
    so G_l = sum_a (xs_a[t] p_a[s] - xs_a[s] p_a[t]). Unknown (a, k) is column a*d + k, c is column n*d."""

    def __init__(self, xs, E, d):
        self.d, self.xs, self.E = d, xs, E
        self.n, self.m, self.g = len(xs), len(E), d * (d - 1) // 2
        self.nu = self.n * d + 1
        self.neq = self.n + self.m + self.g
        self.ST = [(s, t) for s in range(d) for t in range(s + 1, d)]
        self.Ea = np.array([e[0] for e in E], int)
        self.Eb = np.array([e[1] for e in E], int)

    def _pattern(self):
        n, d, m = self.n, self.d, self.m
        ar = np.arange(d)
        rS = np.repeat(np.arange(n), d)
        cS = np.arange(n * d)
        rC = np.repeat(n + np.arange(m), 2 * d)
        cC = np.concatenate([self.Ea[:, None] * d + ar, self.Eb[:, None] * d + ar], 1).ravel()
        return rS, cS, rC, cC

    def jac(self, P):
        """TRUSTED (exact). The Jacobian at points P (n x d floats): every entry is 2 P, P, -1 or a coordinate of xs, so
        it is exact in floating point."""
        n, d, m = self.n, self.d, self.m
        rS, cS, rC, cC = self._pattern()
        vS = 2.0 * P.ravel()
        vC = np.concatenate([P[self.Eb], P[self.Ea]], 1).ravel()
        rCc = n + np.arange(m)
        cCc = np.full(m, n * d)
        vCc = -np.ones(m)
        rG, cG, vG = [], [], []
        for l, (s, t) in enumerate(self.ST):
            row = n + m + l
            rG += [np.full(n, row), np.full(n, row)]
            cG += [np.arange(n) * d + s, np.arange(n) * d + t]
            vG += [self.xs[:, t], -self.xs[:, s]]
        rows = np.concatenate([rS, rC, rCc] + rG)
        cols = np.concatenate([cS, cC, cCc] + cG)
        vals = np.concatenate([vS, vC, vCc] + vG)
        return sps.csr_matrix((vals, (rows, cols)), shape=(self.neq, self.nu))

    def jac_radius(self, rad):
        """TRUSTED. Entrywise upper bounds of |J(z) - J(Xf)| when every coordinate (a, k) of z is within rad[a, k] of Xf:
        2 rad on the S rows, rad on the C rows, 0 on the c column and the G rows (constant)."""
        n, d, m = self.n, self.d, self.m
        rS, cS, rC, cC = self._pattern()
        vS = up(2.0 * rad.ravel())
        vC = np.concatenate([rad[self.Eb], rad[self.Ea]], 1).ravel()
        rows = np.concatenate([rS, rC])
        cols = np.concatenate([cS, cC])
        vals = np.concatenate([vS, vC])
        return sps.csr_matrix((vals, (rows, cols)), shape=(self.neq, self.nu))


class Exact:
    """TRUSTED. Exact residual F(z~) for z~ = (Z / 2^K, C / 2^K) with integer Z (n x d lists) and integer C."""

    def __init__(self, system, K):
        self.s, self.K = system, K
        ratios = [[float(v).as_integer_ratio() for v in row] for row in system.xs]
        self.Lx = max(den.bit_length() - 1 for row in ratios for _, den in row)
        self.Xs = [[num * ((1 << self.Lx) // den) for num, den in row] for row in ratios]  # xs = Xs / 2^Lx exactly

    def residual(self, Z, C):
        """Floats F_f (correctly rounded F) and Ferr >= |F - F_f| componentwise; also the exact numerators."""
        s, K = self.s, self.K
        out = []
        one = 1 << (2 * K)
        for a in range(s.n):
            za = Z[a]
            out.append((sum(z * z for z in za) - one, 2 * K))
        CK = C << K
        for e, (a, b) in enumerate(s.E):
            za, zb = Z[a], Z[b]
            out.append((sum(x * y for x, y in zip(za, zb)) - CK, 2 * K))
            if e % 256 == 0:
                PACE.tick()
        for st, (ss, t) in enumerate(s.ST):
            acc = 0
            for a in range(s.n):
                acc += self.Xs[a][t] * Z[a][ss] - self.Xs[a][ss] * Z[a][t]
            out.append((acc, self.Lx + K))
            PACE.tick()
        Ff = np.array([num / (1 << e) for num, e in out])  # int / int: correctly rounded
        Ferr = np.abs(Ff) * 2.0 ** -52 + ETA
        return Ff, Ferr


# ============================================================================================ UNTRUSTED: refinement
def paced_inverse(A, nb=96, chunk=192, rchunk=48, overwrite=False):
    """UNTRUSTED. Approximate inverse by a blocked LU factorisation with partial pivoting (right-looking, as LAPACK's
    getrf) and blocked triangular solves, in steps short enough for Pace. Any approximate inverse would do. With
    overwrite=True the float64 array A is factorised in place (it is destroyed), which saves one n x n array."""
    import scipy.linalg as sl
    if not (overwrite and isinstance(A, np.ndarray) and A.dtype == np.float64):
        A = np.array(A, dtype=float)
    n = A.shape[0]
    perm = np.arange(n)
    for j0 in range(0, n, nb):
        j1 = min(j0 + nb, n)
        lu, piv = sl.lu_factor(A[j0:, j0:j1], check_finite=False)
        A[j0:, j0:j1] = lu
        for i, p in enumerate(piv):
            if p != i:
                a, b = j0 + i, j0 + int(p)
                A[[a, b], :j0] = A[[b, a], :j0]
                A[[a, b], j1:] = A[[b, a], j1:]
                perm[[a, b]] = perm[[b, a]]
        PACE.tick()
        if j1 < n:
            A[j0:j1, j1:] = sl.solve_triangular(A[j0:j1, j0:j1], A[j0:j1, j1:], lower=True, unit_diagonal=True,
                                                check_finite=False)
            for r0 in range(j1, n, chunk):
                A[r0:r0 + chunk, j1:] -= A[r0:r0 + chunk, j0:j1] @ A[j0:j1, j1:]
                PACE.tick()
    if not np.all(np.isfinite(A)) or np.min(np.abs(np.diag(A))) == 0.0:
        raise np.linalg.LinAlgError("singular")
    R = np.empty((n, n))
    for c0 in range(0, n, rchunk):
        cols = np.arange(c0, min(c0 + rchunk, n))
        B = (perm[:, None] == cols[None, :]).astype(float)  # columns of P = I[perm]
        Y = sl.solve_triangular(A, B, lower=True, unit_diagonal=True, check_finite=False)
        PACE.tick()
        R[:, c0:c0 + len(cols)] = sl.solve_triangular(A, Y, lower=False, check_finite=False)
        PACE.tick()
    return R


def refine(ex, Z, C, R, K, iters=12, log=None):
    """UNTRUSTED. Newton's chord method with exact residuals: z <- z - R F(z), rounded to the grid 2^-K."""
    s = ex.s
    hist = []
    for it in range(iters):
        Ff, _ = ex.residual(Z, C)
        res = float(np.abs(Ff).max())
        hist.append(res)
        if res < 2.0 ** -(K - 6) or (len(hist) > 2 and res > 0.5 * hist[-2]):
            break
        step = -(R @ Ff)
        dz = np.rint(step * 2.0 ** K)
        for a in range(s.n):
            row = Z[a]
            base = a * s.d
            for k in range(s.d):
                row[k] += int(dz[base + k])
        C += int(dz[s.n * s.d])
    if log is not None:
        log["newton_residuals"] = [float(f"{h:.3e}") for h in hist]
    return Z, C, hist[-1]


def place_rattlers(Xpub, V, Xf, c, d, near=3e-6, steps=(1e-6, 1e-7, 1e-8, 1e-9)):
    """UNTRUSTED. Positions for the rattlers (points outside V). A rattler stays where the code has it unless one of its
    pairs is within `near` of c. Then, one rattler at a time (the others held where they are), it moves along a great
    circle by a small angle t in a unit tangent direction v. First choice: a linear programme for v that lowers every
    pair within 1e-11 of c to first order, <v, x_j> <= -s, while the looser pairs within `near` rise by at most a third
    of their gap, t <v, x_j> <= gap / 3 (t from `steps`, largest first). Otherwise, when c > 0, a v orthogonal to the
    tight pairs, which lowers each of their inner products g to g cos t < g, with t limited by the looser pairs. The gap
    check decides whether the result is good."""
    from scipy.optimize import linprog
    N = len(Xpub)
    inV = np.zeros(N, bool)
    inV[V] = True
    rat = np.nonzero(~inV)[0]
    Y = Xpub.copy()
    Y[V] = Xf
    info = {"rattlers": int(len(rat)), "moved_first_order": 0, "moved_second_order": 0, "stuck": 0}
    Yn = Y / np.linalg.norm(Y, axis=1, keepdims=True)
    stuck = []
    for r in rat:
        PACE.tick()
        g = Yn @ Yn[r]
        g[r] = -2.0
        js = np.nonzero(g > c - near)[0]
        if not len(js):
            continue
        x = Yn[r]
        A = Yn[js] - np.outer(Yn[js] @ x, x)  # tangent parts of the neighbours
        gap = c - g[js]
        tight = gap <= 1e-11
        v, step, kind = None, 0.0, None
        for t in steps:
            A_ub = np.hstack([A, np.where(tight, 1.0, 0.0)[:, None]])  # tight: <a, v> + s <= 0; loose: <a, v> <= gap/3t
            b_ub = np.where(tight, 0.0, gap / (3 * t))
            res = linprog(np.r_[np.zeros(d), -1.0], A_ub=A_ub, b_ub=b_ub, A_eq=np.r_[x, 0.0][None, :], b_eq=[0.0],
                          bounds=[(-1.0, 1.0)] * d + [(None, 1.0)], method="highs")
            if res.status == 0 and res.x[-1] > 1e-6:
                v, step, kind = res.x[:d], t, "first"
                info["moved_first_order"] += 1
                break
        if v is None and c > 0 and tight.any():
            Mt = np.vstack([A[tight], x[None, :]])
            _, sv, Vt = np.linalg.svd(Mt)
            if len(sv) < d or sv[-1] < 1e-9 * sv[0]:
                w = Vt[-1] - (Vt[-1] @ x) * x
                w /= np.linalg.norm(w)
                best = (0.0, 1)
                for sign in (1, -1):  # the step must not close the looser pairs: t <w, x_j> < gap_j / 3
                    rate = sign * (A[~tight] @ w)
                    lim = min([steps[0]] + [gap[~tight][i] / (3 * rate[i]) for i in range(len(rate)) if rate[i] > 0])
                    best = max(best, (lim, sign))
                if best[0] > 0:
                    v, step, kind = best[1] * w, best[0], "second"
                    info["moved_second_order"] += 1
        if v is None:
            info["stuck"] += 1
            stuck.append(r)
            continue
        v = v - (v @ x) * x
        if kind == "second":  # the unit direction w: along the great circle
            Yn[r] = math.cos(step) * x + math.sin(step) * v
        else:  # the linear programme's v: x + t v, normalised (changes t <v, x_j> to first order)
            y = x + step * v
            Yn[r] = y / np.linalg.norm(y)
    if stuck:  # rattlers that block one another: move them together (one linear programme for all of them)
        mov = list(stuck)
        col = {r: i for i, r in enumerate(mov)}
        rows, rhs, tight_rows = [], [], []
        for r in mov:
            g = Yn @ Yn[r]
            g[r] = -2.0
            for j in np.nonzero(g > c - near)[0]:
                if j in col and j < r:
                    continue  # pair already listed from j's side
                row = np.zeros(len(mov) * d + 1)
                row[col[r] * d:(col[r] + 1) * d] = Yn[j] - (Yn[j] @ Yn[r]) * Yn[r]
                if j in col:
                    row[col[j] * d:(col[j] + 1) * d] = Yn[r] - (Yn[r] @ Yn[j]) * Yn[j]
                gap = c - g[j]
                if gap <= 1e-11:
                    row[-1] = 1.0
                    rhs.append(0.0)
                else:
                    rhs.append(gap / (3 * steps[0]))
                rows.append(row)
        A_eq = np.zeros((len(mov), len(mov) * d + 1))
        for r, i in col.items():
            A_eq[i, i * d:(i + 1) * d] = Yn[r]
        res = linprog(np.r_[np.zeros(len(mov) * d), -1.0], A_ub=np.array(rows), b_ub=np.array(rhs), A_eq=A_eq,
                      b_eq=np.zeros(len(mov)), bounds=[(-1.0, 1.0)] * (len(mov) * d) + [(None, 1.0)], method="highs")
        info["joint_margin"] = float(res.x[-1]) if res.status == 0 else None
        if res.status == 0 and res.x[-1] > 1e-6:
            v = res.x[:-1].reshape(len(mov), d)
            for r, i in col.items():
                y = Yn[r] + steps[0] * v[i]
                Yn[r] = y / np.linalg.norm(y)
            info["moved_jointly"] = len(mov)
            info["stuck"] = 0
    Y[rat] = Yn[rat]
    return Y, info


# ============================================================================================ TRUSTED: the checks
def krawczyk(sysm, Ff, Ferr, Jc, R, dfx, log):
    """TRUSTED. Rump 2010, Theorem 13.3 with X = [-r, r]^N: checks -R F(z~) + (I - R J(z~ + X)) X in int(X) through the
    componentwise bound Ybar + rs r < r, where Ybar >= |R F(z~)| and rs >= the row sums of |I - R J| over the box.
    Returns (ok, r, err) with err >= |z* - z~| componentwise for the unique zero z* in z~ + X."""
    N = sysm.nu
    k = N + 1
    g = gam(k)
    absJc = abs(Jc)
    Ybar = np.empty(N)
    for i0 in range(0, N, BLOCK):
        Rb = R[i0:i0 + BLOCK]
        RF = (Rb * Ff[None, :]).sum(1)
        A = (np.abs(Rb) * np.abs(Ff)[None, :]).sum(1)
        B = (np.abs(Rb) * Ferr[None, :]).sum(1)
        Ybar[i0:i0 + BLOCK] = up(np.abs(RF) + up(g * up(A)) + up(B))
        PACE.tick()
    r = float(up(2.0 * Ybar.max()))
    for attempt in range(4):
        rad = up(dfx + r)
        Jr = sysm.jac_radius(rad)
        rs = np.empty(N)
        for i0 in range(0, N, BLOCK):
            Rb = R[i0:i0 + BLOCK]
            aR = np.abs(Rb)
            P = np.asarray(Rb @ Jc)
            Eh = -P
            idx = np.arange(Eh.shape[0])
            Eh[idx, i0 + idx] = 1.0 - P[idx, i0 + idx]
            A1 = np.asarray(aR @ absJc)
            A2 = np.asarray(aR @ Jr)
            Ebar = up(up(np.abs(Eh)) + up(g * up(A1)) + up(A2))
            rs[i0:i0 + BLOCK] = up(Ebar.sum(1))
            PACE.tick()
        err = up(Ybar + up(rs * r))
        log.setdefault("krawczyk_attempts", []).append({"r": r, "max_rowsum": float(rs.max()),
                                                        "max_err_over_r": float((err / r).max())})
        if np.all(err < r):
            return True, r, err, Jr, rs
        if rs.max() >= 1.0:
            return False, r, err, Jr, rs
        r = float(up(4.0 * Ybar.max() / (1.0 - float(rs.max()))))
    return False, r, err, Jr, rs


def stress_enclosure(sysm, Jc, Jr, R, log):
    """TRUSTED (the approximate solution is untrusted). Encloses y* = J(z*)^{-T} b, b = -e_c, for every J in Jc +- Jr:
    with T = R^T, |y* - y~| <= w + cs ||y* - y~||_inf, ||y* - y~||_inf <= max(w) / (1 - beta), where w >= |R|^T |b - J^T y~|
    and cs >= the column sums of |I - J R| (beta = max cs < 1). Returns (y~, bound) or (None, None)."""
    N = sysm.nu
    k = N + 1
    g = gam(k)
    b = np.zeros(N)
    b[N - 1] = -1.0
    yt = -R[N - 1, :].copy()
    for _ in range(3):  # untrusted refinement
        rho = b - Jc.T @ yt
        yt = yt + R.T @ rho
    absJc = abs(Jc)
    t = Jc.T @ yt
    rh = b - t
    S1 = absJc.T @ np.abs(yt)
    S2 = Jr.T @ np.abs(yt)
    Resbar = up(up(np.abs(rh)) + up(g * up(S1)) + up(S2))
    w = np.empty(N)
    cs = np.empty(N)
    for j0 in range(0, N, BLOCK):
        Rc = R[:, j0:j0 + BLOCK]
        aR = np.abs(Rc)
        w[j0:j0 + BLOCK] = up((aR * Resbar[:, None]).sum(0))
        Q = np.asarray(Jc @ Rc)
        Eh = -Q
        idx = np.arange(Eh.shape[1])
        Eh[j0 + idx, idx] = 1.0 - Q[j0 + idx, idx]
        A1 = np.asarray(absJc @ aR)
        A2 = np.asarray(Jr @ aR)
        Ebar = up(up(np.abs(Eh)) + up(g * up(A1)) + up(A2))
        cs[j0:j0 + BLOCK] = up(Ebar.sum(0))
        PACE.tick()
    beta = float(cs.max())
    log["stress_beta"] = beta
    if not beta < 1.0:
        return None, None
    den = (1.0 - beta) * (1.0 - 2.0 ** -50)
    dinf = up(up(w.max()) / den)
    bound = up(w + up(cs * dinf))
    return yt, bound


def span_check(Xf, rad, d):
    """TRUSTED (the choice of points and the inverse are untrusted). Picks d points by pivoted QR and proves that every
    matrix within rad of their rows is non-singular: max row sum of an upper bound of |I - Rb B*| < 1."""
    import scipy.linalg as sl
    _, _, piv = sl.qr(Xf.T, pivoting=True, mode="economic")
    idx = np.sort(piv[:d])
    B = Xf[idx]
    Rb = np.linalg.inv(B)
    P = (Rb[:, :, None] * B[None, :, :]).sum(1)
    Eh = -P
    Eh[np.arange(d), np.arange(d)] = 1.0 - P[np.arange(d), np.arange(d)]
    A1 = (np.abs(Rb)[:, :, None] * np.abs(B)[None, :, :]).sum(1)
    A2 = (np.abs(Rb)[:, :, None] * rad[idx][None, :, :]).sum(1)
    Ebar = up(up(np.abs(Eh)) + up(gam(d) * up(A1)) + up(A2))
    rs = up(Ebar.sum(1))
    return bool(rs.max() < 1.0), float(rs.max()), [int(i) for i in idx]


def gap_check(sysm, V, E_glob, Y, Xf, err_pts, Z, C, K, err_max, c_lo, c_err, d):
    """TRUSTED. Every pair i < j of the N points not in E_glob has <p*_i, p*_j> < c*. Points of V: p* within
    dfx + err of Xf in each coordinate; other points: p*_r = Y_r / |Y_r|, within |1 - |Y_r|^2| of Y_r (Euclidean).
    Floating-point filter (Lemma S) first; pairs it cannot decide among V are checked exactly from Z."""
    N = len(Y)
    inV = np.zeros(N, bool)
    inV[V] = True
    A = Y.copy()
    A[V] = Xf
    eps = np.zeros(N)
    eps[V] = err_pts
    for r in np.nonzero(~inV)[0]:
        PACE.tick()
        q = sum(Fraction(float(x)) ** 2 for x in A[r])
        eps[r] = frac_up(abs(1 - q))
    nrm2 = up((A * A).sum(1))
    amax = float(up(math.sqrt(float(nrm2.max()))))
    emax = float(eps.max())
    errG = float(up(up(gam(d) * up(amax * amax)) + up(2.0 * amax * emax) + up(emax * emax)))
    thr = float(up(errG))
    contacts = set(E_glob)
    undecided, min_gap, worst = [], math.inf, None
    for i0 in range(0, N, BLOCK):
        Gb = (A[i0:i0 + BLOCK, None, :] * A[None, :, :]).sum(-1)
        margin = c_lo - Gb
        PACE.tick()
        for ii in range(Gb.shape[0]):
            i = i0 + ii
            if ii % 32 == 0:
                PACE.tick()
            js = np.nonzero(margin[ii, i + 1:] <= thr)[0] + i + 1
            for j in js:
                if (i, int(j)) not in contacts:
                    undecided.append((i, int(j)))
            mrow = margin[ii, i + 1:]
            if len(mrow):
                cand = mrow.copy()
                for j in js:
                    if (i, int(j)) in contacts:
                        cand[j - i - 1] = np.inf
                jm = int(np.argmin(cand))
                if cand[jm] < min_gap:
                    min_gap, worst = float(cand[jm]), (i, i + 1 + jm)
    pos = {int(v): a for a, v in enumerate(V)}
    exact_ok, failed = 0, []
    two = Fraction(1, 1 << K)
    cstar_lo = Fraction(C, 1 << K) - Fraction(c_err)
    e = Fraction(err_max)
    for i, j in undecided:
        PACE.tick()
        if i in pos and j in pos:  # both in V: p* within e of z~ = Z / 2^K in every coordinate
            za, zb = Z[pos[i]], Z[pos[j]]
            s = Fraction(sum(x * y for x, y in zip(za, zb)), 1 << (2 * K))
            l1 = Fraction(sum(abs(x) for x in za) + sum(abs(y) for y in zb), 1 << K)
            if s + e * l1 + d * e * e < cstar_lo:
                exact_ok += 1
                continue
        elif cstar_lo > 0:  # a rattler p*_r = Y_r / |Y_r|: decide <Y_r, q> < c* |Y_r| exactly (q the other point)
            if i in pos:
                i, j = j, i
            yr = [Fraction(float(x)) for x in A[i]]
            nr2 = sum(x * x for x in yr)
            if j in pos:
                zj = [Fraction(z, 1 << K) for z in Z[pos[j]]]
                lhs = sum(x * y for x, y in zip(yr, zj)) + e * sum(abs(x) for x in yr)
                rhs2 = cstar_lo * cstar_lo * nr2
            else:
                ys = [Fraction(float(x)) for x in A[j]]
                lhs = sum(x * y for x, y in zip(yr, ys))
                rhs2 = cstar_lo * cstar_lo * nr2 * sum(x * x for x in ys)
            if lhs <= 0 or lhs * lhs < rhs2:
                exact_ok += 1
                continue
        failed.append((i, j))
    return not failed, {"min_gap_float": min_gap - thr, "closest_noncontact": worst, "undecided_by_float": len(undecided),
                        "decided_exactly": exact_ok, "failed": failed[:10], "n_failed": len(failed)}


# ============================================================================================ driver
def certify(path, d, n, tol=1e-10, K=100, log=None):
    log = {} if log is None else log
    cpu0 = time.process_time()
    try:
        return _certify(path, d, n, tol, K, log)
    finally:
        log["cpu_seconds"] = round(time.process_time() - cpu0, 2)


def _certify(path, d, n, tol, K, log):
    t0 = time.time()
    Xpub = read_code(path, d, n)
    V, E, M, _ = contact_structure(Xpub, d, tol)
    nV, m = len(V), len(E)
    iso = nV * (d - 1) - d * (d - 1) // 2 + 1
    log.update(d=d, n=n, file=os.path.basename(path), tol=tol, bits=K, n_nonrattlers=int(nV), contacts=int(m),
               isostatic_count=int(iso), rattlers=int(n - nV), max_cos_double=M)
    if m != iso:
        log["status"] = "not square: contacts - isostatic count = %d" % (m - iso)
        return log
    sysm = System(Xpub[V], E, d)
    X0 = Xpub[V] / np.linalg.norm(Xpub[V], axis=1, keepdims=True)
    c0 = float(np.mean([X0[a] @ X0[b] for a, b in E]))
    Z = [[int(v) for v in row] for row in (X0 * 2.0 ** K)]
    C = int(c0 * 2.0 ** K)
    ex = Exact(sysm, K)
    J0 = sysm.jac(X0).toarray()
    t1 = time.time()
    try:
        R = paced_inverse(J0, overwrite=True)
    except (np.linalg.LinAlgError, ValueError):
        log["status"] = "Jacobian numerically singular"
        return log
    del J0
    log["inverse_norm_inf"] = float(np.abs(R).sum(1).max())  # ||R||_inf, R ~ J^-1 (untrusted, for the log)
    Z, C, res = refine(ex, Z, C, R, K, log=log)
    t2 = time.time()
    # ---- trusted part
    Ff, Ferr = ex.residual(Z, C)
    Xf = np.array([[z / (1 << K) for z in row] for row in Z])
    dfx = np.abs(Xf) * 2.0 ** -52 + ETA
    Jc = sysm.jac(Xf)
    ok_k, r, err, Jr, rs = krawczyk(sysm, Ff, Ferr, Jc, R, dfx, log)
    log["krawczyk"] = bool(ok_k)
    log["residual_max"] = float(np.abs(Ff).max())
    if not ok_k:
        log["status"] = "FAIL existence (Krawczyk)"
        return log
    err_pts = up(dfx + err[:nV * d].reshape(nV, d))
    c_err = float(err[-1])
    cq = Fraction(C, 1 << K)
    c_lo, c_hi = frac_down(cq - Fraction(c_err)), frac_up(cq + Fraction(c_err))
    log["c_mid"] = str(C * 10 ** 40 // (1 << K))  # 40 decimals of c~ (truncated)
    log["c_err"] = c_err
    log["coord_err_max"] = float(err[:-1].max())
    t3 = time.time()
    yt, ybound = stress_enclosure(sysm, Jc, Jr, R, log)
    if yt is None:
        log["status"] = "FAIL stress enclosure"
        return log
    om = yt[nV:nV + m]
    ob = ybound[nV:nV + m]
    ok_s = bool(np.all(om > ob))
    log.update(stress_min=float(om.min()), stress_max=float(om.max()), stress_err_max=float(ob.max()),
               stress_positive=ok_s, stress_sum=float(om.sum()), slice_multipliers_max=float(np.abs(yt[nV + m:]).max()),
               stress_quantiles_over_mean=[float(x) for x in np.quantile(om * m, [0, 0.01, 0.1, 0.5, 0.9, 0.99, 1])])
    t4 = time.time()
    ok_span, span_rs, span_idx = span_check(Xf, err_pts, d)
    log.update(span=ok_span, span_rowsum=span_rs)
    Y, rinfo = place_rattlers(Xpub, V, Xf, float(cq), d)
    log["rattler_placement"] = rinfo
    E_glob = sorted((int(V[a]), int(V[b])) for a, b in E)
    ok_g, ginfo = gap_check(sysm, V, E_glob, Y, Xf, np.sqrt(d) * err_pts.max(axis=1) * UP_F, Z, C, K,
                            float(err.max() + dfx.max()) * UP_F, c_lo, c_err, d)
    log["gaps"] = ginfo
    log["gaps_ok"] = ok_g
    # distances to the code as published (floats read from the file; the decimal strings differ by < 2^-53 |x|)
    dV = float(up(np.abs(Xf - Xpub[V]).max() + err_pts.max() + 2.0 ** -52))
    rat = np.setdiff1d(np.arange(n), V)
    dR = 0.0
    if len(rat):  # q_r = Y_r / |Y_r| is within |1 - |Y_r|^2| of Y_r
        e_r = max(frac_up(abs(1 - sum(Fraction(float(x)) ** 2 for x in Y[r]))) for r in rat)
        dR = float(up(up(np.abs(Y[rat] - Xpub[rat]).max() + e_r) + 2.0 ** -52))
    log.update(dist_nonrattlers=dV, dist_rattlers=dR, c_lo=repr(c_lo), c_hi=repr(c_hi))
    lo = (cq - Fraction(c_err)) * 10 ** 25
    hi = (cq + Fraction(c_err)) * 10 ** 25
    log["c_interval_25"] = [f"0.{(lo.numerator // lo.denominator):025d}", f"0.{(-(-hi.numerator // hi.denominator)):025d}"]
    ok = ok_k and ok_s and ok_span and ok_g
    log["status"] = "PASS" if ok else "FAIL " + ",".join(
        x for x, f in (("stress", ok_s), ("span", ok_span), ("gaps", ok_g)) if not f)
    log["seconds"] = {"setup_inverse": round(t1 - t0, 2), "inverse_refine": round(t2 - t1, 2),
                      "krawczyk": round(t3 - t2, 2), "stress": round(t4 - t3, 2), "total": round(time.time() - t0, 2)}
    return log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("code")
    ap.add_argument("d", type=int)
    ap.add_argument("n", type=int)
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--bits", type=int, default=100)
    ap.add_argument("--json")
    ap.add_argument("--duty", type=float, default=1.0, help="CPU share for a light background run (see Pace)")
    a = ap.parse_args()
    global PACE
    PACE = Pace(a.duty)
    log = certify(a.code, a.d, a.n, a.tol, a.bits)
    line = (f"{log['file']} d={a.d} n={a.n} {log['status']} nonrattlers={log['n_nonrattlers']} contacts={log['contacts']}"
            f" rattlers={log['rattlers']}")
    if "c_interval_25" in log:
        line += f" c in [{log['c_interval_25'][0]}, {log['c_interval_25'][1]}] stress_min={log['stress_min']:.3e}"
        line += f" min_gap={log['gaps']['min_gap_float']:.2e} secs={log['seconds']['total']}"
    print(line, flush=True)
    if a.json:
        with open(a.json, "w") as f:
            json.dump(log, f, indent=1, default=str)


if __name__ == "__main__":
    main()
