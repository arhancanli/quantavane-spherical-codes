"""Independent re-check, in exact integer and rational arithmetic, of the census certificates (paper version 4).
Written separately from sphcert.py, cq.py and strict.py (it imports none of them) and with no floating-point error
analysis: every quantity that enters a decision is an exact integer or rational. Floating point is used only to
propose candidates (Newton steps, approximate inverses, the subset E', the coefficients of the combined stress).

  python3 -I certify/census/exact_big.py CODE.txt D N TOL [--hyper | --flex | --strict] [--json OUT]

Default (isostatic graph at TOL):
 1. Contacts (pairs within TOL of the largest cosine) and rattlers (fewer than d contacts, removed repeatedly) are
    proposed in double precision.
 2. Unknowns z = (p_1..p_n, c); equations |p_i|^2 = 1, <p_a, p_b> = c (ab in E), and the slice
    sum_i (xi_i[t] p_i[s] - xi_i[s] p_i[t]) = 0 (s < t), with xi the code normalised in double precision.
 3. Existence (Rump 2010, Thm 13.3, with the exact midpoint z~ on the grid 2^-K and an approximate inverse C' rounded to
    a fixed-point grid): |C' F(z~)| + r (|I - C' J(z~)| + |C'| Jrad(r)) 1 < r componentwise, checked exactly.
 4. Stress: y* = J(z*)^-T b enclosed exactly, b = -e_c (isostatic) and positive on E.
 5. Gaps: every pair of points of V outside E has <p*_a, p*_b> < c*, checked exactly from z~ and the enclosure;
    and every rattler at its stored position (normalised exactly) has every inner product below c* (exactly).
--hyper: E is hyperstatic; E' is chosen (untrusted), 3 is applied to E', and with coefficients lambda_f > 0 (untrusted
    least squares) the combined stress J^-T (-e_c - sum_f lambda_f grad F_f(z*)) is enclosed exactly (with the
    enclosure of z* in the right-hand side) and checked positive; tau' (the deviation of the other contacts) is
    enclosed exactly; the sharpness constant kappa, the radius rho and the displacement 4 sqrt(n d) tau'/kappa are
    recomputed exactly, and the local minimum is confirmed when the displacement is below rho; gaps as in 5, with the
    margin of the hyperstatic certificate.
--flex: E = all pairs within TOL (no rattlers removed); with a right inverse C' of the row matrix L (rounded to a
    fixed-point grid), max row sum of |I - L C'| < 1 and h = (1 - beta) / (2 ||C'||_inf) are computed exactly.
--strict: the inner products of the exact configuration P' of the hyperstatic certificate are recognised as rationals
    (continued fractions; untrusted); then exactly: unit diagonal, contacts equal to c*, other pairs below c*, the
    matrix has rank d and is positive semidefinite (G = G_VB G_BB^-1 G_BV, G_BB positive definite), a strictly
    positive equilibrium stress exists (exact rational solution), and the rigidity matrix in the coordinates
    w_a = sum_b u_ab p_b (b in B) has rank n d - d(d-1)/2 over F_p (a lower bound for its rank over Q). This is the
    classical criterion at an exact configuration, in exact arithmetic only (rational Gram matrices).
"""
import json, math, os, sys, time, warnings
from fractions import Fraction

sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/flintlib"))
warnings.filterwarnings("ignore", message=".*encountered in matmul")
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spl

K = 100
P61 = (1 << 61) - 1


def abs_sum(a, axis=None):
    """Exact sum of |a| for an int64 array with |a| < 2^62, as Python ints (split into 31-bit halves)."""
    b = np.abs(a)
    hi = (b >> 31).sum(axis=axis)
    lo = (b & ((1 << 31) - 1)).sum(axis=axis)
    if axis is None:
        return (int(hi) << 31) + int(lo)
    return [(int(h) << 31) + int(l) for h, l in zip(hi, lo)]


def read(path, d, n):
    rows, dec = [], []
    for line in open(path):
        p = [t.strip() for t in line.strip().split(",")]
        if len(p) == d:
            try:
                rows.append([float(x) for x in p])
                dec.append(p)
            except ValueError:
                pass
    X = np.array(rows)
    assert X.shape == (n, d)
    return X / np.linalg.norm(X, axis=1, keepdims=True), dec


def pairs_within(X, tol):
    G = X @ X.T
    np.fill_diagonal(G, -2.0)
    M = G.max()
    return G >= M - tol


def contacts(X, d, tol):
    A = pairs_within(X, tol)
    alive = np.ones(len(X), bool)
    while True:
        deg = (A & alive[:, None] & alive[None, :]).sum(1)
        drop = alive & (deg < d)
        if not drop.any():
            break
        alive &= ~drop
    V = np.nonzero(alive)[0]
    pos = {int(v): k for k, v in enumerate(V)}
    E = [(pos[int(i)], pos[int(j)]) for i, j in zip(*np.nonzero(np.triu(A, 1))) if alive[i] and alive[j]]
    return V, E


def isqrt_up(q):
    """A rational >= sqrt(q) for a rational q >= 0."""
    a, b = q.numerator, q.denominator
    s = math.isqrt(a * b << 120)
    return Fraction(s + 1, b << 60)


def jac_float(P, E, xi, d):
    """Jacobian at P (n x d floats): rows norms, contacts, slice; columns p (n d), c."""
    n = len(P)
    g = d * (d - 1) // 2
    rows, cols, vals = [], [], []
    for a in range(n):
        for k in range(d):
            rows.append(a); cols.append(a * d + k); vals.append(2 * P[a, k])
    for e, (a, b) in enumerate(E):
        for k in range(d):
            rows += [n + e, n + e]; cols += [a * d + k, b * d + k]; vals += [P[b, k], P[a, k]]
        rows.append(n + e); cols.append(n * d); vals.append(-1.0)
    l = 0
    for s in range(d):
        for t in range(s + 1, d):
            for a in range(n):
                rows += [n + len(E) + l, n + len(E) + l]; cols += [a * d + s, a * d + t]
                vals += [xi[a, t], -xi[a, s]]
            l += 1
    return sps.csr_matrix((vals, (rows, cols)), shape=(n + len(E) + g, n * d + 1))


# ============================================================================================ existence (exact)
def existence(xi, E, d, log):
    """Steps 2-3 for the square system of E on the points xi (the code's non-rattlers). Returns a dict with z~ (Z, Cz),
    the enclosure eps, C' (Cint / 2^tC), the exact Jacobian J(z~) (by rows and columns over 2^L) and the row sums."""
    nv, m = len(xi), len(E)
    N = nv * d + 1
    fr = [[Fraction(float(v)) for v in row] for row in xi]
    Lx = max(f.denominator.bit_length() - 1 for row in fr for f in row)
    XI = [[int(f * (1 << Lx)) for f in row] for row in fr]
    P = xi.copy()
    c = float(np.mean([P[a] @ P[b] for a, b in E]))

    def Ffloat(P, c):
        f = np.empty(N)
        f[:nv] = (P * P).sum(1) - 1
        f[nv:nv + m] = np.array([P[a] @ P[b] for a, b in E]) - c
        l = 0
        for s in range(d):
            for t in range(s + 1, d):
                f[nv + m + l] = (xi[:, t] * P[:, s] - xi[:, s] * P[:, t]).sum()
                l += 1
        return f

    for _ in range(6):
        J = jac_float(P, E, xi, d).tocsc()
        step = spl.splu(J).solve(-Ffloat(P, c))
        P = P + step[:-1].reshape(nv, d)
        c = c + step[-1]
    Z = [[int(round(v * 2.0 ** K)) for v in row] for row in P]
    Cz = int(round(c * 2.0 ** K))

    def Fexact(Z, Cz):
        f = []
        one = 1 << (2 * K)
        for a in range(nv):
            f.append(Fraction(sum(z * z for z in Z[a]) - one, one))
        for a, b in E:
            f.append(Fraction(sum(x * y for x, y in zip(Z[a], Z[b])) - (Cz << K), one))
        for s in range(d):
            for t in range(s + 1, d):
                acc = sum(XI[a][t] * Z[a][s] - XI[a][s] * Z[a][t] for a in range(nv))
                f.append(Fraction(acc, 1 << (K + Lx)))
        return f

    Pf = np.array([[z / 2.0 ** K for z in row] for row in Z])
    Cinv = np.linalg.inv(jac_float(Pf, E, xi, d).toarray())
    for _ in range(4):
        f = np.array([float(v) for v in Fexact(Z, Cz)])
        if np.abs(f).max() < 2.0 ** -(K - 4):
            break
        dz = -Cinv @ f
        for a in range(nv):
            for k in range(d):
                Z[a][k] += int(round(dz[a * d + k] * 2.0 ** K))
        Cz += int(round(dz[-1] * 2.0 ** K))
    Fz = Fexact(Z, Cz)
    log["residual_max"] = float(max(abs(v) for v in Fz))
    tC = 62 - int(math.ceil(math.log2(np.abs(Cinv).max() + 1e-300)))
    Cint = np.rint(Cinv * 2.0 ** tC).astype(np.int64)
    del Cinv
    L = max(K, Lx)
    rowsJ, colsJ, valsJ = [], [], []
    for a in range(nv):
        for k in range(d):
            rowsJ.append(a); colsJ.append(a * d + k); valsJ.append(2 * Z[a][k] << (L - K))
    for e, (a, b) in enumerate(E):
        for k in range(d):
            rowsJ += [nv + e, nv + e]; colsJ += [a * d + k, b * d + k]
            valsJ += [Z[b][k] << (L - K), Z[a][k] << (L - K)]
        rowsJ.append(nv + e); colsJ.append(N - 1); valsJ.append(-(1 << L))
    l = 0
    for s in range(d):
        for t in range(s + 1, d):
            for a in range(nv):
                rowsJ += [nv + m + l, nv + m + l]; colsJ += [a * d + s, a * d + t]
                valsJ += [XI[a][t] << (L - Lx), -XI[a][s] << (L - Lx)]
            l += 1
    bycol = [[] for _ in range(N)]
    byrow = [[] for _ in range(N)]
    for i, j, v in zip(rowsJ, colsJ, valsJ):
        if v:
            bycol[j].append((i, v))
            byrow[i].append((j, v))
    S_C_rows = nv + m
    denF = 1 << (2 * K + Lx + K)
    Fint = [int(v * denF) for v in Fz]
    assert all(Fraction(fi, denF) == v for fi, v in zip(Fint, Fz))
    CF = []
    for i0 in range(0, N, 256):
        blk = Cint[i0:i0 + 256].astype(object)
        CF += list(blk.dot(np.array(Fint, dtype=object)))
    absCF = [Fraction(abs(v), denF << tC) for v in CF]
    scale = 1 << (tC + L)
    rowI, rowRad, rowCE = [], [], []
    for i0 in range(0, N, 256):
        blk = Cint[i0:i0 + 256].astype(object)
        B = blk.shape[0]
        Pm = np.zeros((B, N), dtype=object)
        for j in range(N):
            acc = None
            for k, v in bycol[j]:
                term = blk[:, k] * v
                acc = term if acc is None else acc + term
            if acc is not None:
                Pm[:, j] = acc
        for ii in range(B):
            Pm[ii, i0 + ii] -= scale
        rowI += [sum(abs(x) for x in Pm[ii]) for ii in range(B)]
        rowRad += abs_sum(Cint[i0:i0 + B, :S_C_rows], axis=1)
        rowCE += abs_sum(Cint[i0:i0 + B, nv:], axis=1)       # columns of the contact and slice equations
        del Pm
    r = 2 * max(absCF) or Fraction(1, 1 << (2 * K))      # (a zero residual: any small radius)
    ok, lhs, rs = False, None, None
    for attempt in range(3):
        rs = [Fraction(rowI[i], scale) + 2 * d * r * Fraction(rowRad[i], 1 << tC) for i in range(N)]
        lhs = [absCF[i] + r * rs[i] for i in range(N)]
        if all(x < r for x in lhs):
            ok = True
            break
        r *= 4
    log["existence"] = ok
    return dict(ok=ok, Z=Z, Cz=Cz, eps=lhs, r=r, rs=rs, Cint=Cint, tC=tC, L=L, byrow=byrow, bycol=bycol, N=N, nv=nv,
                m=m, scale=scale, rowCE=rowCE, Lx=Lx, XI=XI)


def transposed(ex, b, brad, epsmax, E, d, log):
    """Exact enclosure of y* = J(z*)^-T b' for every b' within brad of b (b, brad: lists of Fractions), with J(z*) for
    every z* within epsmax of z~ in the point coordinates. Returns (y~, dev) with |y* - y~| <= dev, or (None, None)."""
    N, nv, Cint, tC, L, byrow = ex["N"], ex["nv"], ex["Cint"], ex["tC"], ex["L"], ex["byrow"]
    scale = ex["scale"]
    degs = [0] * nv
    for a, bb in E:
        degs[a] += 1
        degs[bb] += 1
    cJ = 2 + max(degs)
    colI, colRad = [0] * N, [0] * N
    for j0 in range(0, N, 256):
        blk = Cint[:, j0:j0 + 256].astype(object)
        B = blk.shape[1]
        Qm = np.zeros((N, B), dtype=object)
        for i in range(N):
            acc = None
            for k, v in byrow[i]:
                term = blk[k, :] * v
                acc = term if acc is None else acc + term
            if acc is not None:
                Qm[i, :] = acc
        for jj in range(B):
            Qm[j0 + jj, jj] -= scale
        cr = abs_sum(Cint[:, j0:j0 + B], axis=0)
        for jj in range(B):
            colI[j0 + jj] = sum(abs(x) for x in Qm[:, jj])
            colRad[j0 + jj] = cr[jj]
        del Qm
    beta = max(Fraction(colI[j], scale) + cJ * epsmax * Fraction(colRad[j], 1 << tC) for j in range(N))
    log["stress_beta"] = float(beta)
    if not beta < 1:
        return None, None
    # y~ = C'^T b (exact)
    yt = [Fraction(0)] * N
    nzb = [(k, bk) for k, bk in enumerate(b) if bk != 0]
    for i in range(N):
        yt[i] = sum((Fraction(int(Cint[k, i]), 1 << tC) * bk for k, bk in nzb), Fraction(0))
    res = list(b)
    for i in range(N):
        if yt[i] == 0:
            continue
        for j, v in byrow[i]:
            res[j] -= Fraction(v, 1 << L) * yt[i]
    rad = list(brad)
    for a in range(nv):
        for k in range(d):
            rad[a * d + k] += 2 * epsmax * abs(yt[a])
    for e, (a, bb) in enumerate(E):
        w = epsmax * abs(yt[nv + e])
        for k in range(d):
            rad[a * d + k] += w
            rad[bb * d + k] += w
    absres = [abs(res[j]) + rad[j] for j in range(N)]
    num = max(absres) * Fraction(max(abs_sum(Cint, axis=0)), 1 << tC)
    return yt, num / (1 - beta)


def gaps_exact(ex, V, E_all, X, dec, n, d, c_lo, log):
    """Every pair of V outside E_all below c_lo at every z* within eps of z~ (exactly), and every rattler at its stored
    position, normalised exactly, below c_lo with every other point (exactly, by squares)."""
    Z, eps = ex["Z"], ex["eps"]
    e = max(eps[:-1])
    Es = set(E_all)
    nv = len(V)
    bad = 0
    two = 1 << (2 * K)
    for a in range(nv):
        za = Z[a]
        la = Fraction(sum(abs(x) for x in za), 1 << K)
        for b in range(a + 1, nv):
            if (a, b) in Es:
                continue
            zb = Z[b]
            s = Fraction(sum(x * y for x, y in zip(za, zb)), two)
            if s + e * (la + Fraction(sum(abs(y) for y in zb), 1 << K)) + d * e * e >= c_lo:
                bad += 1
    log["gaps_V_failed"] = bad
    rat = [i for i in range(n) if i not in set(int(v) for v in V)]
    rbad = 0
    pos = {int(v): k for k, v in enumerate(V)}
    for r in rat:
        q = [Fraction(t) for t in dec[r]]
        nq = sum(x * x for x in q)
        for j in range(n):
            if j == r or (j in rat and j < r):
                continue
            if j in pos:
                zj = [Fraction(z, 1 << K) for z in Z[pos[j]]]
                lhs = sum(x * y for x, y in zip(q, zj)) + e * sum(abs(x) for x in q)
                rhs2 = c_lo * c_lo * nq
            else:
                qj = [Fraction(t) for t in dec[j]]
                lhs = sum(x * y for x, y in zip(q, qj))
                rhs2 = c_lo * c_lo * nq * sum(x * x for x in qj)
            ok = (lhs <= 0 or lhs * lhs < rhs2) if c_lo > 0 else (lhs < 0 and lhs * lhs > rhs2)
            if not ok:
                rbad += 1
    log["rattlers"] = len(rat)
    log["gaps_rattlers_failed_at_stored_positions"] = rbad
    return bad == 0


# ============================================================================================ the modes
def run_iso(X, dec, d, n, tol, log):
    V, E = contacts(X, d, tol)
    nv, m = len(V), len(E)
    g = d * (d - 1) // 2
    log.update(nV=nv, m=m, iso=nv * (d - 1) - g + 1)
    assert m == log["iso"], "not isostatic"
    ex = existence(X[V], E, d, log)
    log["unknowns"] = ex["N"]
    if not ex["ok"]:
        log["status"] = "FAIL existence"
        return log
    eps, N = ex["eps"], ex["N"]
    cmid = Fraction(ex["Cz"], 1 << K)
    log["c_lo_e30"] = math.floor((cmid - eps[-1]) * 10 ** 30)
    log["c_hi_e30"] = math.ceil((cmid + eps[-1]) * 10 ** 30)
    b = [Fraction(0)] * N
    b[N - 1] = Fraction(-1)
    yt, dev = transposed(ex, b, [Fraction(0)] * N, max(eps[:-1]), E, d, log)
    if yt is None:
        log["status"] = "FAIL stress beta"
        return log
    om = yt[nv:nv + m]
    log["stress_min_lower"] = float(min(om) - dev)
    log["stress_positive"] = bool(min(om) > dev)
    ok_g = gaps_exact(ex, V, E, X, dec, n, d, cmid - eps[-1], log)
    log["status"] = ("PASS existence+stress" if log["stress_positive"] else "FAIL stress") + \
        (", gaps of V" if ok_g else ", FAIL gaps of V")
    return log


def choose_subset(X, E, d, iso):
    """UNTRUSTED. E' by QR with column pivoting on the contact rows projected off the norm and slice rows."""
    import scipy.linalg as sl
    J = jac_float(X, E, X, d).tocsr()
    nv, m = len(X), len(E)
    Bm = sps.vstack([J[:nv], J[nv + m:]]).toarray()
    Qb, _ = np.linalg.qr(Bm.T)
    Cr = J[nv:nv + m].toarray()
    Cr = Cr - (Cr @ Qb) @ Qb.T
    _, piv = sl.qr(Cr.T, mode="r", pivoting=True)
    return sorted(int(e) for e in piv[:iso])


def run_hyper(X, dec, d, n, tol, log):
    V, E = contacts(X, d, tol)
    nv, m = len(V), len(E)
    g = d * (d - 1) // 2
    iso = nv * (d - 1) - g + 1
    log.update(nV=nv, m=m, iso=iso)
    assert m > iso, "not hyperstatic"
    xi = X[V]
    sub = choose_subset(xi, E, d, iso)
    Esq = [E[e] for e in sub]
    extra = [E[e] for e in range(m) if e not in set(sub)]
    ex = existence(xi, Esq, d, log)
    N = ex["N"]
    log["unknowns"] = N
    if not ex["ok"]:
        log["status"] = "FAIL existence"
        return log
    Z, eps = ex["Z"], ex["eps"]
    e = max(eps[:-1])
    cmid = Fraction(ex["Cz"], 1 << K)
    # coefficients lambda_f (untrusted): least squares towards the uniform stress, from float solves
    Pf = np.array([[z / 2.0 ** K for z in row] for row in Z])
    Jd = jac_float(Pf, Esq, xi, d).toarray()
    Ji = np.linalg.inv(Jd)
    y0 = -Ji[N - 1]
    Jf = np.zeros((len(extra), N))
    for f, (a, b) in enumerate(extra):
        Jf[f, a * d:(a + 1) * d] = Pf[b]
        Jf[f, b * d:(b + 1) * d] = Pf[a]
        Jf[f, N - 1] = -1.0
    Y = Jf @ Ji
    om0, Ysub = y0[nv:nv + iso], Y[:, nv:nv + iso]
    u = 1.0 / m
    lam = np.linalg.solve(Ysub @ Ysub.T + np.eye(len(extra)), Ysub @ (om0 - u) + u)
    log["coefficients"] = "least squares"
    if lam.min() <= 0 or (om0 - Ysub.T @ lam).min() <= 0:
        from scipy.optimize import linprog             # untrusted: maximise the smallest entry of the combined stress
        kx = len(extra)
        A_ub = np.vstack([np.hstack([Ysub.T, np.ones((iso, 1))]), np.hstack([-np.eye(kx), np.ones((kx, 1))])])
        lp = linprog(np.r_[np.zeros(kx), -1.0], A_ub=A_ub, b_ub=np.r_[om0, np.zeros(kx)],
                     bounds=[(0, None)] * kx + [(None, 1.0)], method="highs")
        if lp.status != 0 or lp.x[-1] <= 0:
            log["status"] = "no positive coefficients (least squares and linear program)"
            return log
        lam = lp.x[:kx]
        log["coefficients"] = "linear program"
    lamq = [Fraction(float(x)) for x in lam]
    # b = -e_c - sum_f lambda_f grad F_f(z~), brad = sum_f lambda_f e (the gradients at z* within eps)
    b = [Fraction(0)] * N
    brad = [Fraction(0)] * N
    b[N - 1] = Fraction(-1) + sum(lamq)
    for f, (a, bb) in enumerate(extra):
        for k in range(d):
            b[a * d + k] -= lamq[f] * Fraction(Z[bb][k], 1 << K)
            b[bb * d + k] -= lamq[f] * Fraction(Z[a][k], 1 << K)
            brad[a * d + k] += lamq[f] * e
            brad[bb * d + k] += lamq[f] * e
    yt, dev = transposed(ex, b, brad, e, Esq, d, log)
    if yt is None:
        log["status"] = "FAIL stress beta"
        return log
    om = yt[nv:nv + iso]
    wmin = min(min(om) - dev, min(lamq))
    log["stress_min_lower"] = float(wmin)
    if not (min(om) > dev and min(lamq) > 0):
        log["status"] = "FAIL combined stress"
        return log
    # tau' exactly
    two = 1 << (2 * K)
    tau = Fraction(0)
    for a, bb in extra:
        s = Fraction(sum(x * y for x, y in zip(Z[a], Z[bb])), two)
        l1 = Fraction(sum(abs(x) for x in Z[a]) + sum(abs(y) for y in Z[bb]), 1 << K)
        tau = max(tau, abs(s - cmid) + e * l1 + d * e * e + eps[-1])
    log["tau_extra"] = float(tau)
    # kappa, rho and the displacement, exactly (Lemma on the inverse and the quantitative theorem of the paper)
    rs = ex["rs"]
    betaK = max(rs)
    Cnorm = Fraction(max(abs_sum(ex["Cint"], axis=1)), 1 << ex["tC"])
    mu = Cnorm / (1 - betaK)
    mu_p = max(Fraction(ex["rowCE"][i], 1 << ex["tC"]) + rs[i] * mu for i in range(nv * d))
    dxi = max(abs(Fraction(z, 1 << K) - Fraction(float(xi[a, k]))) for a in range(nv) for k, z in enumerate(Z[a])) + e
    theta = 2 * nv * dxi * mu_p
    if not theta < 1:
        log["status"] = "FAIL theta"
        return log
    sigma = (1 - theta) / mu_p
    kappa = sigma * wmin / (1 - wmin)
    Lam = 1 + abs(cmid) + eps[-1] + tau
    sq_d = isqrt_up(Fraction(d))
    rho = min(kappa / (4 * Lam * sq_d), Fraction(1, 2))
    disp = 4 * isqrt_up(Fraction(nv * d)) * tau / kappa
    log.update(kappa=float(kappa), rho=float(rho), localmin_displacement=float(disp))
    ok_g = gaps_exact(ex, V, E, X, dec, n, d, cmid - eps[-1] - 2 * tau - 4 * sq_d * tau / kappa, log)
    log["status"] = ("LOCALMIN confirmed" if disp < rho else "FAIL displacement") + \
        (", gaps of V" if ok_g else ", FAIL gaps of V")
    return log


def run_flex(X, dec, d, n, tol, log):
    A = pairs_within(X, tol)
    I, J = np.nonzero(np.triu(A, 1))
    m = len(I)
    pts = sorted(set(I.tolist()) | set(J.tolist()))
    loc = {p: k for k, p in enumerate(pts)}
    nE = len(pts)
    Nn, Nc = nE + m, nE * d
    log.update(contacts=m, points_in_contacts=nE, rows=Nn)
    rows, cols, vals = [], [], []
    for p in pts:
        for k in range(d):
            rows.append(loc[p]); cols.append(loc[p] * d + k); vals.append(2.0 * X[p, k])
    for e, (a, b) in enumerate(zip(I, J)):
        for k in range(d):
            rows += [nE + e, nE + e]; cols += [loc[a] * d + k, loc[b] * d + k]; vals += [X[b, k], X[a, k]]
    L = sps.csr_matrix((vals, (rows, cols)), shape=(Nn, Nc))
    import scipy.linalg as sl
    Qf, Rf = sl.qr(L.T.toarray(), mode="economic")
    C = sl.solve_triangular(Rf, Qf.T, lower=False).T            # untrusted right inverse
    tC = 60 - int(math.ceil(math.log2(np.abs(C).max() + 1e-300)))
    Ci = np.rint(C * 2.0 ** tC).astype(np.int64)
    # exact: L has float entries (rationals over 2^Lq)
    fr = [Fraction(float(v)) for v in L.data]
    Lq = max(f.denominator.bit_length() - 1 for f in fr)
    Lint = sps.csr_matrix(([int(f * (1 << Lq)) for f in fr], L.indices, L.indptr), shape=L.shape, dtype=object)
    scale = 1 << (Lq + tC)
    rowsum = []
    for i in range(Nn):
        acc = np.zeros(Nn, dtype=object)
        for jj in range(Lint.indptr[i], Lint.indptr[i + 1]):
            acc = acc + Ci[Lint.indices[jj]].astype(object) * Lint.data[jj]
        acc[i] -= scale
        rowsum.append(sum(abs(x) for x in acc))
    beta = Fraction(max(rowsum), scale)
    cnorm = Fraction(max(abs_sum(Ci, axis=1)), 1 << tC)
    log["flex_beta"] = float(beta)
    if not beta < 1:
        log["status"] = "FAIL flex"
        return log
    h = (1 - beta) / (2 * cnorm)
    log["flex_radius_l1"] = float(h)
    log["status"] = "FLEX confirmed"
    return log


def run_strict(X, dec, d, n, tol, log):
    """Rational Gram matrices only."""
    V, E = contacts(X, d, tol)
    nv, m = len(V), len(E)
    g = d * (d - 1) // 2
    iso = nv * (d - 1) - g + 1
    log.update(nV=nv, m=m, iso=iso)
    xi = X[V]
    sub = choose_subset(xi, E, d, iso) if m > iso else list(range(m))
    Esq = [E[e] for e in sub]
    ex = existence(xi, Esq, d, log)
    if not ex["ok"]:
        log["status"] = "FAIL existence"
        return log
    Z = ex["Z"]
    two = 1 << (2 * K)
    # untrusted: rationals of small height for every inner product of P'
    G = [[None] * nv for _ in range(nv)]
    for a in range(nv):
        G[a][a] = Fraction(1)
        for b in range(a + 1, nv):
            v = Fraction(sum(x * y for x, y in zip(Z[a], Z[b])), two)
            q = v.limit_denominator(10 ** 6)
            if abs(q - v) > Fraction(1, 10 ** 24):
                log["status"] = "an inner product not recognised as a rational"
                return log
            G[a][b] = G[b][a] = q
    c = G[E[0][0]][E[0][1]]
    Es = set(E)
    # exact: contacts equal c, the other pairs below c
    for a in range(nv):
        for b in range(a + 1, nv):
            if ((a, b) in Es) != (G[a][b] == c) or G[a][b] > c:
                log["status"] = "contacts or gaps fail in the rational Gram matrix"
                return log
    # exact: B (untrusted choice), G_BB positive definite (LDL^T pivots), G = G_VB G_BB^-1 G_BV
    import scipy.linalg as sl
    _, _, piv = sl.qr(xi.T, pivoting=True, mode="economic")
    B = sorted(int(i) for i in piv[:d])
    M = [[G[i][j] for j in B] for i in B]
    piv_ok = True
    for k in range(d):
        if M[k][k] <= 0:
            piv_ok = False
            break
        for i in range(k + 1, d):
            f = M[i][k] / M[k][k]
            for j in range(k, d):
                M[i][j] -= f * M[k][j]
    if not piv_ok:
        log["status"] = "G_BB not positive definite"
        return log
    # W = G_BB^-1 G_BV by exact elimination
    Aug = [[G[i][j] for j in B] + [G[i][v] for v in range(nv)] for i in B]
    for k in range(d):
        p = Aug[k][k]
        Aug[k] = [x / p for x in Aug[k]]
        for i in range(d):
            if i != k and Aug[i][k] != 0:
                f = Aug[i][k]
                Aug[i] = [x - f * y for x, y in zip(Aug[i], Aug[k])]
    W = [row[d:] for row in Aug]
    for a in range(nv):
        for b in range(a, nv):
            if sum(G[a][B[k]] * W[k][b] for k in range(d)) != G[a][b]:
                log["status"] = "the rational Gram matrix does not have rank d"
                return log
    # exact: a strictly positive stress at the Gram level: sum_b w_ab G_bk = lam_a G_ak for k in B, sum w = 1.
    # untrusted: the uniform stress first, else an LP point; exact: rref with the free unknowns fixed
    nbr = {a: [] for a in range(nv)}
    for e_, (a, b) in enumerate(E):
        nbr[a].append((b, e_))
        nbr[b].append((a, e_))
    def stress_ok(w):
        for a in range(nv):
            s = [sum(w[e_] * G[b][k] for b, e_ in nbr[a]) for k in B]
            lam = s[0] / G[a][B[0]] if G[a][B[0]] != 0 else None
            if lam is None:
                kk = next(i for i, k in enumerate(B) if G[a][k] != 0)
                lam = s[kk] / G[a][B[kk]]
            if any(s[i] != lam * G[a][k] for i, k in enumerate(B)):
                return False
        return True
    w = [Fraction(1, m)] * m
    kind = "uniform" if stress_ok(w) else None
    if kind is None:
        w = gram_stress(G, W, B, E, nv, d, m, log)
        if w is None or not stress_ok(w):
            log["status"] = "no exact positive stress found"
            return log
        kind = "rational (exact solution near a linear program's)"
    # exact: bar-rigidity, rank over F_q of the rigidity matrix in the coordinates u (w_a = sum_b u_ab p_b)
    Q31 = (1 << 31) - 1

    def modp(q):
        return q.numerator % Q31 * pow(q.denominator, -1, Q31) % Q31
    gm = {}
    for a in range(nv):
        for kk, b in enumerate(B):
            gm[(a, kk)] = modp(G[a][b])
    rows = np.zeros((nv + m, nv * d), dtype=np.int64)
    for a in range(nv):
        for kk in range(d):
            rows[a, a * d + kk] = gm[(a, kk)]
    for e_, (a, a2) in enumerate(E):
        for kk in range(d):
            rows[nv + e_, a2 * d + kk] = (rows[nv + e_, a2 * d + kk] + gm[(a, kk)]) % Q31
            rows[nv + e_, a * d + kk] = (rows[nv + e_, a * d + kk] + gm[(a2, kk)]) % Q31
    rk = rank_mod_p(rows, nv * d)
    del rows
    log.update(stress=kind, rank_mod_p=rk, rank_needed=nv * d - g)
    log["status"] = "STRICT confirmed (exact)" if rk == nv * d - g else "FAIL rank"
    return log


def gram_stress(G, W, B, E, nv, d, m, log):
    """A strictly positive rational stress at the Gram level. Untrusted: a linear program (floating point) for the
    stress of largest smallest entry, with the equilibrium written in the coordinates W (p_v = sum_k W[k][v] p_{B_k}):
    sum_b w_ab W[.][b] - lam_a W[.][a] = 0, sum w = 1. Exact: the reduced row echelon form of these equations over Q
    (python-flint), the free unknowns fixed at the rounded LP values and the pivot unknowns solved exactly."""
    import flint
    from scipy.optimize import linprog
    ncol = m + nv
    Wf = np.array([[float(W[k][v]) for v in range(nv)] for k in range(d)])
    Aeq = np.zeros((nv * d + 1, ncol))
    for e_, (a, b) in enumerate(E):
        Aeq[a * d:(a + 1) * d, e_] += Wf[:, b]
        Aeq[b * d:(b + 1) * d, e_] += Wf[:, a]
    for a in range(nv):
        Aeq[a * d:(a + 1) * d, m + a] = -Wf[:, a]
    Aeq[nv * d, :m] = 1.0
    lp = linprog(np.r_[np.zeros(ncol), -1.0], A_ub=np.hstack([-np.eye(m), np.zeros((m, nv)), np.ones((m, 1))]),
                 b_ub=np.zeros(m), A_eq=np.hstack([Aeq, np.zeros((nv * d + 1, 1))]), b_eq=np.r_[np.zeros(nv * d), 1.0],
                 bounds=[(None, None)] * ncol + [(None, 1.0)], method="highs")
    if lp.status != 0 or lp.x[-1] <= 1e-12:
        return None
    z0 = [Fraction(float(x)).limit_denominator(10 ** 12) for x in lp.x[:ncol]]
    fq = lambda q: flint.fmpq(q.numerator, q.denominator)
    Ab = flint.fmpq_mat(nv * d + 1, ncol + 1)
    for e_, (a, b) in enumerate(E):
        for k in range(d):
            Ab[a * d + k, e_] = Ab[a * d + k, e_] + fq(W[k][b])
            Ab[b * d + k, e_] = Ab[b * d + k, e_] + fq(W[k][a])
        Ab[nv * d, e_] = 1
    for a in range(nv):
        for k in range(d):
            Ab[a * d + k, m + a] = -fq(W[k][a])
    Ab[nv * d, ncol] = 1
    R, rk = Ab.rref()
    piv = []
    for i in range(rk):
        j = next(j for j in range(ncol + 1) if R[i, j] != 0)
        if j == ncol:
            return None
        piv.append(j)
    ps = set(piv)
    z = [fq(q) for q in z0]
    for i, j in enumerate(piv):
        s = R[i, ncol]
        for jj in range(ncol):
            if jj not in ps and R[i, jj] != 0:
                s = s - R[i, jj] * z[jj]
        z[j] = s
    w = [Fraction(int(x.p), int(x.q)) for x in z[:m]]
    if min(w) <= 0:
        return None
    log["stress_min_exact_rel"] = float(min(w)) * m
    return w


def rank_mod_p(rows, ncols, p=None):
    """The rank over F_q, q = 2^31 - 1, of an integer matrix given by rows of residues (Gaussian elimination in int64;
    products of two residues stay below 2^62). A lower bound for the rank over Q of a rational matrix whose
    denominators are prime to q."""
    q = (1 << 31) - 1
    A = np.asarray(rows, dtype=np.int64)
    A %= q
    rank = 0
    nr = A.shape[0]
    for c in range(ncols):
        if rank == nr:
            break
        nz = np.nonzero(A[rank:, c])[0]
        if len(nz) == 0:
            continue
        pr = rank + int(nz[0])
        if pr != rank:
            A[[rank, pr]] = A[[pr, rank]]
        inv = pow(int(A[rank, c]), q - 2, q)
        A[rank] = (A[rank] * inv) % q
        col = A[:, c].copy()
        col[rank] = 0
        idx = np.nonzero(col)[0]
        piv = A[rank][None, :]
        for k0 in range(0, len(idx), 512):      # in chunks, to bound the temporaries
            ii = idx[k0:k0 + 512]
            A[ii] = (A[ii] - (col[ii, None] * piv) % q) % q
        rank += 1
    return rank


def main():
    path, d, n, tol = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
    out = sys.argv[sys.argv.index("--json") + 1] if "--json" in sys.argv else None
    t0 = time.time()
    X, dec = read(path, d, n)
    import hashlib
    log = {"d": d, "n": n, "tol": tol,
           "code_sha256": {"exact_big.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()}}
    mode = "iso"
    for f in ("--hyper", "--flex", "--strict"):
        if f in sys.argv:
            mode = f[2:]
    log["mode"] = mode
    {"iso": run_iso, "hyper": run_hyper, "flex": run_flex, "strict": run_strict}[mode](X, dec, d, n, tol, log)
    log["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(log))
    if out:
        json.dump(log, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
