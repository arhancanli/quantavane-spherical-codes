"""A second, separately written check of the existence and stress parts of a certificate, in exact rational arithmetic.

  python3 -I certify/exact_check.py CODE.txt D N [--tol 1e-10] [--bits 100]

It takes the same candidates as sphcert.py (contacts, the refined point z~ = Z / 2^K, the floating-point approximate
inverse R and the approximate stress y~; how they were found does not matter) and re-checks, with Python integers
only, (1) Rump's Theorem 13.3 for the box z~ + [-r, r]^N:
      |R F(z~)| + (|I - R J(z~)| + |R| Jrad) r 1 < r     componentwise,
where Jrad bounds |J(z) - J(z~)| for |z - z~| <= r (2r on the rows |p_a|^2 - 1, r on the rows <p_a, p_b> - c), and
(2) the enclosure of the stress y* = J(z*)^-T (-e_c):
      beta = max column sum of (|I - J(z~) R| + Jrad |R|) < 1,  w = |R|^T (|b - J(z~)^T y~| + Jrad^T |y~|),
      |y* - y~| <= w + colsum * max(w) / (1 - beta),  and y~_e > that bound on every contact e.
No floating-point operation enters a decision: floats are converted to exact rationals first. It is slow (pure
Python, about N * nnz(J) big-integer products) and meant for the smallest codes, to cross-check sphcert.py.
"""
import os, sys, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import sphcert as S


def ints_of(vals):
    """Exact integers M and a shift L with vals = M / 2^L (vals: floats)."""
    rat = [float(v).as_integer_ratio() for v in vals]
    L = max(den.bit_length() - 1 for _, den in rat)
    return [num << (L - (den.bit_length() - 1)) for num, den in rat], L


def main():
    path, d, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    tol = float(sys.argv[sys.argv.index("--tol") + 1]) if "--tol" in sys.argv else 1e-10
    K = int(sys.argv[sys.argv.index("--bits") + 1]) if "--bits" in sys.argv else 100
    t0, cpu0 = time.time(), time.process_time()
    # ---- candidates (untrusted), exactly as sphcert.certify makes them
    Xpub = S.read_code(path, d, n)
    V, E, M, _ = S.contact_structure(Xpub, d, tol)
    nV, m = len(V), len(E)
    assert m == nV * (d - 1) - d * (d - 1) // 2 + 1, "not square"
    sysm = S.System(Xpub[V], E, d)
    X0 = Xpub[V] / np.linalg.norm(Xpub[V], axis=1, keepdims=True)
    c0 = float(np.mean([X0[a] @ X0[b] for a, b in E]))
    Z = [[int(v) for v in row] for row in (X0 * 2.0 ** K)]
    C = int(c0 * 2.0 ** K)
    ex = S.Exact(sysm, K)
    R = S.paced_inverse(sysm.jac(X0).toarray())
    Z, C, _ = S.refine(ex, Z, C, R, K)
    Xf = np.array([[z / (1 << K) for z in row] for row in Z])
    Jf = sysm.jac(Xf)
    N = sysm.nu
    b = np.zeros(N)
    b[N - 1] = -1.0
    yt = -R[N - 1, :].copy()
    for _ in range(3):
        yt = yt + R.T @ (b - Jf.T @ yt)
    # ---- exact data. Common scale 2^LJ for J(z~): S rows 2 Z / 2^K, C rows Z / 2^K and -1, G rows Xs / 2^Lx.
    LJ = max(K, ex.Lx)
    J = [dict() for _ in range(N)]  # J[row][col] = integer over 2^LJ
    for a in range(nV):
        for k in range(d):
            J[a][a * d + k] = (2 * Z[a][k]) << (LJ - K)
    for e, (a, bb) in enumerate(E):
        row = J[nV + e]
        for k in range(d):
            row[a * d + k] = Z[bb][k] << (LJ - K)
            row[bb * d + k] = Z[a][k] << (LJ - K)
        row[nV * d] = -(1 << LJ)
    for l, (s, t) in enumerate(sysm.ST):
        row = J[nV + m + l]
        for a in range(nV):
            row[a * d + s] = ex.Xs[a][t] << (LJ - ex.Lx)
            row[a * d + t] = -(ex.Xs[a][s] << (LJ - ex.Lx))
    cols = [[] for _ in range(N)]
    for i, row in enumerate(J):
        for j, v in row.items():
            cols[j].append((i, v))
    Rint, LR = ints_of(R.ravel())
    Rr = [Rint[i * N:(i + 1) * N] for i in range(N)]
    Fnum = []  # exact F as integers over 2^LF
    Ff, _ = ex.residual(Z, C)
    # recompute F exactly with integers (independently of the floats Ff)
    LF = max(2 * K, ex.Lx + K)
    for a in range(nV):
        Fnum.append((sum(z * z for z in Z[a]) - (1 << (2 * K))) << (LF - 2 * K))
    for a, bb in E:
        Fnum.append((sum(x * y for x, y in zip(Z[a], Z[bb])) - (C << K)) << (LF - 2 * K))
    for s, t in sysm.ST:
        acc = sum(ex.Xs[a][t] * Z[a][s] - ex.Xs[a][s] * Z[a][t] for a in range(nV))
        Fnum.append(acc << (LF - ex.Lx - K))
    # radius: twice the largest |R F| (any choice is fine; the test decides)
    RF = []
    for i in range(N):  # |R F| over 2^(LR+LF)
        S.PACE.tick()
        RF.append(abs(sum(Rr[i][j] * Fnum[j] for j in range(N) if Fnum[j])))
    rmax = Fraction(max(RF), 1 << (LR + LF))
    r = 2 * rmax if rmax > 0 else Fraction(1, 1 << (2 * K))
    # Jrad as integers over 2^LJ times r: S entries 2r, C entries r (positions as J, excluding the c column and G rows)
    ok = True
    worst = Fraction(0)
    one = 1 << (LR + LJ)
    for i in range(N):
        S.PACE.tick()
        Ri = Rr[i]
        rowsum_E = 0  # |I - R J| row sum, over 2^(LR+LJ)
        rowsum_rad = 0  # |R| Jrad row sum divided by r, over 2^LR
        for j in range(N):
            acc = 0
            radacc = 0
            for k, v in cols[j]:
                acc += Ri[k] * v
                if k < nV:
                    radacc += abs(Ri[k]) * 2
                elif k < nV + m and j != nV * d:
                    radacc += abs(Ri[k])
            if i == j:
                acc = one - acc
            rowsum_E += abs(acc)
            rowsum_rad += radacc
        lhs = Fraction(RF[i], 1 << (LR + LF)) + (Fraction(rowsum_E, one) + r * Fraction(rowsum_rad, 1 << LR)) * r
        worst = max(worst, lhs / r)
        if not lhs < r:
            ok = False
    t1 = time.time()
    print(f"existence: {'PASS' if ok else 'FAIL'}  r = {float(r):.3e}  max (lhs / r) = {float(worst):.6f}  "
          f"({t1 - t0:.0f} s, CPU {time.process_time() - cpu0:.1f} s)", flush=True)
    if not ok:
        return
    # ---- stress
    Yint, LY = ints_of(yt)
    # residual b - J(z~)^T y~ over 2^(LJ+LY); b = -e_c
    res = []
    for j in range(N):
        S.PACE.tick()
        acc = sum(v * Yint[i] for i, v in cols[j])
        bj = -(1 << (LJ + LY)) if j == N - 1 else 0
        res.append(bj - acc)
    # |res*| <= |res| + Jrad^T |y~|
    resb = []
    for j in range(N):
        radacc = 0
        if j < nV * d:
            for i, _ in cols[j]:
                if i < nV:
                    radacc += 2 * abs(Yint[i])
                elif i < nV + m:
                    radacc += abs(Yint[i])
        resb.append(Fraction(abs(res[j]), 1 << (LJ + LY)) + r * Fraction(radacc, 1 << LY))
    # w = |R|^T resb; column sums of |I - J R| + Jrad |R|
    w = []
    cs = []
    for c in range(N):
        S.PACE.tick()
        w.append(sum(Fraction(abs(Rr[k][c]), 1 << LR) * resb[k] for k in range(N) if resb[k]))
    oneJ = 1 << (LJ + LR)
    for c in range(N):
        S.PACE.tick()
        tot = 0
        totrad = 0
        Rc = [Rr[k][c] for k in range(N)]
        for i in range(N):
            acc = sum(v * Rc[k] for k, v in J[i].items())
            if i == c:
                acc = oneJ - acc
            tot += abs(acc)
            if i < nV:
                totrad += 2 * sum(abs(Rc[k]) for k in J[i])
            elif i < nV + m:
                totrad += sum(abs(Rc[k]) for k in J[i] if k != nV * d)
        cs.append(Fraction(tot, oneJ) + r * Fraction(totrad, 1 << LR))
    beta = max(cs)
    if not beta < 1:
        print(f"stress: FAIL beta = {float(beta)}")
        return
    dinf = max(w) / (1 - beta)
    smin = None
    ok2 = True
    for e in range(m):
        i = nV + e
        bound = w[i] + cs[i] * dinf
        y = Fraction(Yint[i], 1 << LY)
        if not y > bound:
            ok2 = False
        q = float(y - bound)
        smin = q if smin is None else min(smin, q)
    print(f"stress: {'PASS' if ok2 else 'FAIL'}  beta = {float(beta):.3e}  min (y~ - bound) over contacts = {smin:.3e}  "
          f"({time.time() - t1:.0f} s; CPU in all {time.process_time() - cpu0:.1f} s)", flush=True)


if __name__ == "__main__":
    main()
