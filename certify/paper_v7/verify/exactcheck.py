"""Exact re-verification of the census certificates (census version 7): the decision of record.

  nice -n 10 python3 -I exactcheck.py OUTDIR [--shard I/K] [--only d_n,...] [--status PASS,LOCALMIN,NONOPT]

Every decision is made in exact integer and rational arithmetic (FLINT's fmpz matrices and Python integers), with Arb
balls (python-flint) only where a square root enters: the normalisation of the stored code, sqrt(d) and sqrt(nd). No
floating-point bound of the census kernel (sphcert.py, cq.py; Lemma 3.19, Appendix A of the paper) is used.

UNTRUSTED (taken from a re-run of the census certificate, cq.certify with keep): the proposal (V, E, E'), the slice
reference xi (the stored code rounded and normalised in double precision), the grid point z~ = (Z, C) / 2^100, the
approximate inverse C of the square Jacobian, the coefficients lambda_f, the rattler positions, the d points of the
span test. They are only data: whatever is reported as proved is proved whichever data were given.

TRUSTED, exact, for the square system F of the paper (Section 3.1) with E' and the slice reference xi:
  (K)  C is replaced by C' = diag(2^-s_i) Cint with integer Cint (each row rounded to 60 bits; any matrix will do).
       J0 = DF(z~) is exact (scaled to integers), C'J0 is formed exactly (fmpz_mat), F(z~) is exact, and for a dyadic
       radius r: |C'F(z~)| + r (|I - C'J0| + r |C'| Delta(1)) 1 < r 1, compared exactly (Theorem 3.3 with delta = 0).
  (S)  The stress (row gamma of DF(z*)^-1, or b^T DF(z*)^-1 for a hyperstatic graph) is enclosed through
       J^-1 = C' + (I - C'J) C' + (I - C'J)^2 J^-1, with |J^-1| <= ||C'||_inf / (1 - beta_K) entrywise and
       |C'(J - J0)| <= |C'| Delta(r 1); every product is exact. Lower bounds > 0 on the contacts (PASS, LOCALMIN), or
       an upper bound < 0 on one (NONOPT).
  (P)  d points of z~: B^-1 exact (fmpq_mat), and ||B^-1||_inf d r < 1, so every matrix within r of B is nonsingular.
  (G)  The Gram matrix of z~ exactly (fmpz_mat): every pair of V outside E below the bound required, with the
       perturbation r (||p_i||_1 + ||p_j||_1) + d r^2; rattlers decided exactly through squares.
  (Q)  mu_p, the slice term, zeta (the rigidity constant), kappa, Lambda, rho, chi, the hyperstatic inequality
       4 sqrt(nd) chi < kappa rho and the tube 4 sqrt(nd) chi / kappa, all exact (rational bounds of the roots).
  (D)  The least-squares and largest coordinate distance to the exactly normalised stored code, in Arb.
The output records the certified quantities (as decimal strings rounded outward and as floats) and the SHA-256 of this
program.
"""
import hashlib, json, math, os, sys, time, traceback
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, CERT)
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import flint

SHA = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
CBITS = 60            # bits kept in each row of the approximate inverse
BLOCK = 384           # rows of C per exact block product
APREC = 256           # bits of the Arb balls (roots only)


# ============================================================================================ exact helpers
def F2(num, e):
    """The rational num / 2^e."""
    return Fraction(num, 1 << e) if e >= 0 else Fraction(num * (1 << -e))


def dy(x):
    """A double as an exact Fraction."""
    return Fraction(float(x))


def up_float(q):
    """A float >= the rational q (for reporting only)."""
    f = float(q)
    return f if Fraction(f) >= q else math.nextafter(f, math.inf)


def dn_float(q):
    f = float(q)
    return f if Fraction(f) <= q else math.nextafter(f, -math.inf)


def sqrt_up(q, bits=APREC):
    """A rational >= sqrt(q) (q >= 0 rational), from an Arb ball."""
    old = flint.ctx.prec
    flint.ctx.prec = bits
    a = flint.arb(flint.fmpq(q.numerator, q.denominator)).sqrt()
    u = a.upper()
    flint.ctx.prec = old
    m, e = u.mid().man_exp()
    r = Fraction(int(m)) * (Fraction(2) ** int(e))
    assert r * r >= q
    return r


def sqrt_dn(q, bits=APREC):
    old = flint.ctx.prec
    flint.ctx.prec = bits
    a = flint.arb(flint.fmpq(q.numerator, q.denominator)).sqrt()
    lo = a.lower()
    flint.ctx.prec = old
    m, e = lo.mid().man_exp()
    r = Fraction(int(m)) * (Fraction(2) ** int(e))
    if r < 0:
        r = Fraction(0)
    assert r * r <= q
    return r


def dyadic_exp(arr):
    """The smallest L with every double of arr an integer multiple of 2^-L, and the integers."""
    L = 0
    ratios = []
    for x in np.asarray(arr, float).ravel():
        n_, d_ = float(x).as_integer_ratio()
        ratios.append((n_, d_))
        L = max(L, d_.bit_length() - 1)
    ints = [n_ * ((1 << L) // d_) for n_, d_ in ratios]
    return L, ints


def scale_rows(R):
    """Integer rows Cint (|entries| <= 2^60) and shifts s with C' = diag(2^-s) Cint (the approximate inverse used)."""
    U = R.shape[0]
    Cint = np.zeros(R.shape, dtype=np.int64)
    s = np.zeros(U, dtype=np.int64)
    for i in range(U):
        mx = float(np.max(np.abs(R[i])))
        if mx == 0.0 or not math.isfinite(mx):
            raise ValueError("approximate inverse has a zero or non-finite row")
        e = math.frexp(mx)[1]                        # mx < 2^e
        s[i] = CBITS - e
        Cint[i] = np.rint(np.ldexp(R[i], int(s[i]))).astype(np.int64)
    return Cint, s


# ============================================================================================ the square system, exactly
class Square:
    """The square system of the paper for points 0..n-1, contact pairs Esq and slice reference xi (doubles), with the
    centre given exactly as integers Zc / 2^Lc (n x d) and gammac / 2^Lc."""

    def __init__(self, n, d, Esq, xi, Zc, gc, Lc):
        self.n, self.d, self.Esq = n, d, list(Esq)
        self.m = len(self.Esq)
        self.ST = [(s, t) for s in range(d) for t in range(s + 1, d)]
        self.g = len(self.ST)
        self.U = n * d + 1
        assert self.n + self.m + self.g == self.U, "not square"
        self.Lx, xin = dyadic_exp(xi)
        self.Xi = [xin[a * d:(a + 1) * d] for a in range(n)]
        self.Zc, self.gc, self.Lc = Zc, gc, Lc
        self.L = max(self.Lc, self.Lx)
        # the Jacobian at the centre as integers / 2^L, column structure for sparse products
        sc = 1 << (self.L - self.Lc)
        sx = 1 << (self.L - self.Lx)
        cols = [[] for _ in range(self.U)]           # column j: list of (row, value)
        for a in range(n):
            for k in range(d):
                cols[a * d + k].append((a, 2 * Zc[a][k] * sc))
        for e, (a, b) in enumerate(self.Esq):
            r_ = n + e
            for k in range(d):
                cols[a * d + k].append((r_, Zc[b][k] * sc))
                cols[b * d + k].append((r_, Zc[a][k] * sc))
            cols[n * d].append((r_, -(1 << self.L)))
        for l, (s_, t_) in enumerate(self.ST):
            r_ = n + self.m + l
            for a in range(n):
                cols[a * d + s_].append((r_, self.Xi[a][t_] * sx))
                cols[a * d + t_].append((r_, -self.Xi[a][s_] * sx))
        self.cols = cols
        self.deg = [0] * n
        self.inc = [[] for _ in range(n)]            # contact rows at each point
        for e, (a, b) in enumerate(self.Esq):
            self.inc[a].append(n + e)
            self.inc[b].append(n + e)

    def J0_fmpz(self):
        M = flint.fmpz_mat(self.U, self.U)
        for j, lst in enumerate(self.cols):
            for i, v in lst:
                M[i, j] = v
        return M

    def residual(self):
        """F(centre) exactly, as integers / 2^LF."""
        n, d, Lc = self.n, self.d, self.Lc
        LF = max(2 * Lc, self.Lx + Lc)
        out = []
        one = 1 << (2 * Lc)
        for a in range(n):
            out.append((sum(z * z for z in self.Zc[a]) - one) << (LF - 2 * Lc))
        gsh = self.gc << Lc
        for (a, b) in self.Esq:
            out.append((sum(x * y for x, y in zip(self.Zc[a], self.Zc[b])) - gsh) << (LF - 2 * Lc))
        for (s_, t_) in self.ST:
            acc = 0
            for a in range(n):
                acc += self.Xi[a][t_] * self.Zc[a][s_] - self.Xi[a][s_] * self.Zc[a][t_]
            out.append(acc << (LF - self.Lx - Lc))
        return out, LF

    def vecJ0(self, w, ew):
        """The row vector (w / 2^ew) J0, exactly, as integers / 2^(ew + L)."""
        return [sum(w[i] * v for i, v in lst) for lst in self.cols], ew + self.L


class Inverse:
    """The approximate inverse C' = diag(2^-s) Cint, with exact block products."""

    def __init__(self, R):
        self.Cint, self.s = scale_rows(R)
        self.U = R.shape[0]
        self.S = int(self.s.max())
        self.blocks = []
        for i0 in range(0, self.U, BLOCK):
            rows = self.Cint[i0:i0 + BLOCK].tolist()
            self.blocks.append((i0, flint.fmpz_mat(rows)))
        self.rowabs = [int(np.abs(self.Cint[i]).astype(object).sum()) for i in range(self.U)]

    def row(self, i):
        return [int(x) for x in self.Cint[i]], int(self.s[i])

    def vecC(self, w, ew, absolute=False):
        """(w / 2^ew) C' (or |C'|), exactly, as integers / 2^(ew + S): sum over blocks of (w_j 2^(S - s_j)) Cint_j."""
        acc = None
        for i0, B in self.blocks:
            b = B.nrows()
            ww = [int(w[i0 + k]) << (self.S - int(self.s[i0 + k])) for k in range(b)]
            if not any(ww):
                continue
            Bm = B if not absolute else _abs_cache(self, i0, B)
            v = flint.fmpz_mat(1, b, ww) * Bm
            acc = v if acc is None else acc + v
        if acc is None:
            return [0] * self.U, ew + self.S
        return [int(x) for x in acc.entries()], ew + self.S


_ABS = {}


def _abs_cache(inv, i0, B):
    key = (id(inv), i0)
    if key not in _ABS:
        rows = np.abs(inv.Cint[i0:i0 + B.nrows()]).tolist()
        _ABS[key] = flint.fmpz_mat(rows)
    return _ABS[key]


def square_test(sq, inv, mode, rad=None, bspec=None, want_constants=False, contact_cols=None):
    """The exact tests on the square system sq with the approximate inverse inv.
    mode 'krawczyk': find a dyadic r with (K); the box is z~ + [-r, r]^U. mode 'box': rad (n x d doubles) is a given
    componentwise radius of the box around the centre (no existence test; nonsingularity on the box).
    bspec: None for b = -e_gamma; else (extra pairs, lambda doubles) for b = -e_gamma - sum lambda_f grad F_f(p).
    Returns a dict of exact rationals and flags."""
    n, d, m, U = sq.n, sq.d, sq.m, sq.U
    out = {}
    L = sq.L
    # ---- exact block products: A_i = row sums of |I - C'J0|, and (for the stress) nothing else is needed
    J0 = sq.J0_fmpz()
    A = [None] * U
    for i0, B in inv.blocks:
        P = B * J0
        ent = P.entries()
        b = B.nrows()
        for k in range(b):
            i = i0 + k
            si = int(inv.s[i])
            diag = 1 << (si + L)
            rowv = ent[k * U:(k + 1) * U]
            tot = 0
            for j in range(U):
                x = int(rowv[j])
                if j == i:
                    x -= diag
                tot += x if x >= 0 else -x
            A[i] = F2(tot, si + L)
        del P, ent
    del J0
    # ---- |C'| Delta(1) 1: 2d sum over the norm and contact rows (k < n + m) of |C'_ik|  (box mode: Delta(rad) 1)
    if mode == "box":
        Lr, radi = dyadic_exp(rad)
        radI = [radi[a * d:(a + 1) * d] for a in range(n)]
        drow = []                                     # (Delta(rad) 1)_k, integers / 2^Lr
        for a in range(n):
            drow.append(2 * sum(radI[a]))
        for (a, b) in sq.Esq:
            drow.append(sum(radI[a]) + sum(radI[b]))
        drow += [0] * sq.g
        # |C'| Delta(rad) 1 = row-wise sum_k |C'_ik| drow_k: a matrix-vector product with |Cint|
        Bv = []
        col = flint.fmpz_mat(U, 1, drow + [0] * (U - len(drow)))
        for i0, B in inv.blocks:
            v = _abs_cache(inv, i0, B) * col
            for k in range(B.nrows()):
                i = i0 + k
                Bv.append(F2(int(v[k, 0]), int(inv.s[i]) + Lr))
        beta = [A[i] + Bv[i] for i in range(U)]
        betaK = max(beta)
        out["beta_K"] = betaK
        if not betaK < 1:
            out["ok_nonsingular"] = False
            return out
        out["ok_nonsingular"] = True
        r = None
    else:
        Bq = []
        for i in range(U):
            row = inv.Cint[i]
            s_nc = int(np.abs(row[:n + m]).astype(object).sum())
            Bq.append(F2(2 * d * s_nc, int(inv.s[i])))
        Fint, LF = sq.residual()
        col = flint.fmpz_mat(U, 1, Fint)
        CF = []
        for i0, B in inv.blocks:
            v = B * col
            for k in range(B.nrows()):
                i = i0 + k
                CF.append(abs(F2(int(v[k, 0]), int(inv.s[i]) + LF)))
        ymax = max(CF)
        # r: a power of two at least 2 max |C'F|, tried up to three times larger
        ok = False
        e0 = (2 * ymax).denominator.bit_length() - (2 * ymax).numerator.bit_length() if ymax > 0 else 400
        for e in range(e0, e0 - 6, -1):
            r = Fraction(1, 1 << e) if e >= 0 else Fraction(1 << -e)
            if r < 2 * ymax:
                continue
            if all(CF[i] + r * A[i] + r * r * Bq[i] < r for i in range(U)):
                ok = True
                break
        out["ok_K"] = ok
        if not ok:
            return out
        beta = [A[i] + r * Bq[i] for i in range(U)]
        betaK = max(beta)
        eps = [CF[i] + r * beta[i] for i in range(U)]
        out.update(r=r, beta_K=betaK, eps=eps)
        # epsilon-inflation: the largest r' (within a factor 2) for which (K) still holds
        q = min((1 - A[i]) / Bq[i] for i in range(U) if Bq[i] > 0)
        ru = None
        rr = q / 2
        for _ in range(60):
            if rr <= r:
                break
            if all(CF[i] + rr * A[i] + rr * rr * Bq[i] < rr for i in range(U)):
                ru = rr
                break
            rr /= 2
        out["r_unique"] = ru if ru is not None else r
    # ---- uncertainty of J over the box, row-wise: D(v)_j = sum_k |v_k| Delta_kj
    def Dvec(v):
        """For a row vector v (rationals over the equations), the vector |v| Delta over the unknowns."""
        res = [Fraction(0)] * U
        for a in range(n):
            if mode == "box":
                ra = [F2(x, Lr) for x in radI[a]]
                base = 2 * abs(v[a])
                for k in range(d):
                    # Delta_{norm a, (a,k)} = 2 rad_ak; Delta_{contact e, (a,k)} = rad of the other end b, coordinate k
                    t = base * ra[k]
                    for e_ in sq.inc[a]:
                        a2, b2 = sq.Esq[e_ - n]
                        other = b2 if a2 == a else a2
                        t += abs(v[e_]) * F2(radI[other][k], Lr)
                    res[a * d + k] = t
            else:
                t = 2 * abs(v[a]) + sum(abs(v[e_]) for e_ in sq.inc[a])
                for k in range(d):
                    res[a * d + k] = r * t
        return res
    # ---- the right-hand side b (as exact rationals over the equations), and its radius over the box
    b = [Fraction(0)] * U
    b[U - 1] = Fraction(-1)
    brad_sum = Fraction(0)
    lam = []
    if bspec is not None:
        extra, lamf = bspec
        for f, (a, bb) in enumerate(extra):
            lf = dy(lamf[f])
            lam.append(lf)
            for k in range(d):
                b[a * d + k] -= lf * F2(sq.Zc[bb][k], sq.Lc)
                b[bb * d + k] -= lf * F2(sq.Zc[a][k], sq.Lc)
            b[U - 1] += lf
            if mode == "box":
                brad_sum += lf * (sum(F2(x, Lr) for x in radI[a]) + sum(F2(x, Lr) for x in radI[bb]))
            else:
                brad_sum += lf * 2 * d * r
    # NOTE: b is indexed by the unknowns (it is the right-hand side of J^T y = b); y = J^-T b, y^T = b^T J^-1.
    # ---- y^T = b^T J^-1 = b^T C' + b^T (I - C'J) C' + b^T (I - C'J)^2 J^-1, all over the equations index of y
    den_b = 1
    for x in b:
        den_b = max(den_b, x.denominator)
    eb = (den_b - 1).bit_length()
    assert all(((x * (1 << eb)).denominator == 1) for x in b)
    bint = [int(x * (1 << eb)) for x in b]
    bC, ebC = inv.vecC(bint, eb)                      # b^T C' over the equations: integers / 2^ebC
    bCq = [F2(x, ebC) for x in bC]
    # v = b^T (I - C'J0) = b^T - (b^T C') J0  (a vector over the unknowns)
    bCJ, ebCJ = sq.vecJ0(bC, ebC)
    v = [b[j] - F2(bCJ[j], ebCJ) for j in range(U)]
    Db = Dvec(bCq)                                    # |b^T C'| Delta: bound of |b^T C' (J - J0)|
    # first order: v C' +- Db |C'|
    den_v = max(x.denominator for x in v)
    ev_ = (den_v - 1).bit_length()
    vint = [int(x * (1 << ev_)) for x in v]
    assert all(F2(vint[j], ev_) == v[j] for j in range(U))
    vC, evC = inv.vecC(vint, ev_)
    den_D = max(x.denominator for x in Db) if Db else 1
    eD = (den_D - 1).bit_length()
    Dint = [int(x * (1 << eD)) + (0 if (x * (1 << eD)).denominator == 1 else 1) for x in Db]
    DC, eDC = inv.vecC(Dint, eD, absolute=True)
    # second order and the radius of b
    normC = max(F2(inv.rowabs[i], int(inv.s[i])) for i in range(U))
    Jinv = normC / (1 - betaK)                        # >= every |entry| of J^-1 on the box
    t2 = sum((abs(v[k]) + Db[k]) * beta[k] for k in range(U)) * Jinv
    t3 = brad_sum * Jinv
    lo, hi = [], []
    for e in range(m):
        j = n + e                                     # the contact equation e
        mid = bCq[j] + F2(vC[j], evC)
        rad_ = F2(DC[j], eDC) + t2 + t3
        lo.append(mid - rad_)
        hi.append(mid + rad_)
    out.update(stress_lo=lo, stress_hi=hi, lam=lam, normC=normC, Jinv=Jinv)
    if want_constants:
        # mu_p >= max over the unknown rows i < nd of sum_{j in contact and slice equations} |J^-1_ij|
        ind = [0] * n + [1] * (m + sq.g)
        col = flint.fmpz_mat(U, 1, ind)
        nuE = []
        for i0, B in inv.blocks:
            vv = _abs_cache(inv, i0, B) * col
            for k in range(B.nrows()):
                nuE.append(F2(int(vv[k, 0]), int(inv.s[i0 + k])))
        mu_p = max(nuE[i] + beta[i] * Jinv for i in range(n * d))
        out["mu_p"] = mu_p
    return out


# ============================================================================================ census certificates
def span_exact(Zc, Lc, d, rmax):
    """(P): d points chosen by pivoted QR (untrusted); B^-1 exact; ||B^-1||_inf d rmax < 1."""
    import scipy.linalg as sl
    Xf = np.array([[float(F2(z, Lc)) for z in row] for row in Zc])
    _, _, piv = sl.qr(Xf.T, pivoting=True, mode="economic")
    idx = sorted(int(i) for i in piv[:d])
    B = flint.fmpq_mat(d, d, [flint.fmpq(Zc[i][k], 1 << Lc) for i in idx for k in range(d)])
    try:
        Bi = B.inv()
    except Exception:
        return False, None, idx
    nrm = max(sum(abs(Fraction(int(Bi[i, k].p), int(Bi[i, k].q))) for k in range(d)) for i in range(d))
    val = nrm * d * rmax
    return val < 1, val, idx


def gaps_exact(Zc, Lc, d, Eset, rmax, c_lo):
    """(G) on V: every pair outside Eset has <p_i, p_j> < c_lo for all p within rmax (coordinatewise) of Zc / 2^Lc.
    Returns (ok, the smallest proven gap)."""
    n = len(Zc)
    M = flint.fmpz_mat(Zc)
    G = M * M.transpose()
    ent = G.entries()
    l1 = [sum(abs(z) for z in row) for row in Zc]
    worst = None
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in Eset:
                continue
            s = F2(int(ent[i * n + j]), 2 * Lc)
            bound = s + rmax * F2(l1[i] + l1[j], Lc) + d * rmax * rmax
            gap = c_lo - bound
            if worst is None or gap < worst:
                worst = gap
    if worst is None:
        return True, None
    return worst > 0, worst


def rattlers_exact(Y, V, Zc, Lc, rmax, c_lo, d):
    """(G) for the rattlers: q_r = Y_r / |Y_r| below c_lo with every q_j (j in V within rmax of Zc_j / 2^Lc, or another
    rattler), decided exactly through squares. Returns (ok, list of failures)."""
    N = len(Y)
    pos = {int(v): a for a, v in enumerate(V)}
    rat = [r for r in range(N) if r not in pos]
    bad = []
    for ii, r in enumerate(rat):
        yr = [dy(x) for x in Y[r]]
        nr2 = sum(x * x for x in yr)
        for j in range(N):
            if j == r or (j in rat and rat.index(j) < ii):
                continue
            if j in pos:
                zj = [F2(z, Lc) for z in Zc[pos[j]]]
                lhs = sum(x * y for x, y in zip(yr, zj)) + rmax * sum(abs(x) for x in yr)
                # |q_j| <= 1 + d rmax ... q_j is a unit vector within rmax of zj: <Y_r, q_j> <= lhs; need < c_lo |Y_r|
                rhs2 = c_lo * c_lo * nr2
            else:
                ys = [dy(x) for x in Y[j]]
                lhs = sum(x * y for x, y in zip(yr, ys))
                rhs2 = c_lo * c_lo * nr2 * sum(x * x for x in ys)
            if c_lo > 0:
                ok = lhs <= 0 or lhs * lhs < rhs2
            else:
                ok = lhs < 0 and lhs * lhs > rhs2
            if not ok:
                bad.append((r, j))
    return not bad, bad


def distance_arb(Zc, Lc, dec_rows, rmax, Yrat=None, dec_rat=None):
    """Upper bounds of the least-squares and largest coordinate distance between the points (within rmax of
    Zc / 2^Lc) and the exactly normalised stored points dec_rows; rattlers Y_r / |Y_r| against theirs."""
    old = flint.ctx.prec
    flint.ctx.prec = APREC
    tot = flint.arb(0)
    mx = flint.arb(0)
    rm = flint.arb(flint.fmpq(rmax.numerator, rmax.denominator))
    for a, row in enumerate(dec_rows):
        q = [Fraction(t) for t in row]
        nr = flint.arb(flint.fmpq(*_fq(sum(x * x for x in q)))).sqrt()
        for k in range(len(q)):
            xh = flint.arb(flint.fmpq(*_fq(q[k]))) / nr
            diff = abs(flint.arb(flint.fmpq(Zc[a][k], 1 << Lc)) - xh) + rm
            tot += diff * diff
            if diff.upper() > mx.upper():
                mx = diff
    dr = flint.arb(0)
    if Yrat is not None:
        for y, row in zip(Yrat, dec_rat):
            q = [Fraction(t) for t in row]
            nr = flint.arb(flint.fmpq(*_fq(sum(x * x for x in q)))).sqrt()
            yy = [flint.arb(flint.fmpq(*_fq(dy(t)))) for t in y]
            ny = sum((t * t for t in yy), flint.arb(0)).sqrt()
            for k in range(len(q)):
                diff = abs(yy[k] / ny - flint.arb(flint.fmpq(*_fq(q[k]))) / nr)
                if diff.upper() > dr.upper():
                    dr = diff
    D = tot.sqrt()
    res = (_arb_up(D), _arb_up(mx), _arb_up(dr))
    flint.ctx.prec = old
    return res


def _fq(q):
    return int(q.numerator), int(q.denominator)


def _arb_up(a):
    u = a.upper()
    m, e = u.mid().man_exp()
    return Fraction(int(m)) * (Fraction(2) ** int(e))


def check_entry(X, dec, d, n, tol, census_status, log):
    """Re-verifies the census certificate of an entry exactly (statuses PASS, NONOPT, LOCALMIN)."""
    import cq
    t0 = time.time()
    keep = {}
    cl = cq.certify(X, d, n, tol, dec=dec, keep=keep)
    log["census_status_rerun"] = cl.get("status")
    if cl.get("status") != census_status:
        log["status"] = "census certificate not reproduced on re-run"
        return log
    V, E, Esq, extra = keep["V"], keep["E"], keep["Esq"], keep["extra"]
    nV = len(V)
    hyper = len(extra) > 0
    Z, C, K = keep["Z"], keep["C"], keep["K"]
    xs = keep["xs"]
    R = np.asarray(keep["R"], dtype=np.float64)
    log.update(nV=int(nV), contacts=int(len(E)), extra=len(extra), unknowns=int(nV * d + 1), tol=tol,
               E_sha256=hashlib.sha256(json.dumps(sorted([int(V[a]), int(V[b])] for a, b in E)).encode()).hexdigest())
    sq = Square(nV, d, Esq, xs, Z, C, K)
    inv = Inverse(R)
    del R, keep["R"]
    bspec = (extra, keep["beta_f"]) if hyper else None
    res = square_test(sq, inv, "krawczyk", bspec=bspec, want_constants=True)
    _ABS.clear()
    if not res.get("ok_K"):
        log["status"] = "FAIL (K)"
        return log
    r = res["r"]
    log.update(r=float(r), beta_K=up_float(res["beta_K"]), r_unique=dn_float(res["r_unique"]))
    eps_g = res["eps"][-1]
    gam_c = F2(C, K)
    c_lo, c_hi = gam_c - eps_g, gam_c + eps_g
    lo30 = (c_lo * 10 ** 30)
    hi30 = (c_hi * 10 ** 30)
    log["c_interval_30"] = [str(lo30.numerator // lo30.denominator), str(-(-hi30.numerator // hi30.denominator))]
    lo, hi = res["stress_lo"], res["stress_hi"]
    lam = res["lam"]
    om_min = min(lo + lam) if lam else min(lo)
    stress_pos = min(lo) > 0 and (not lam or min(lam) > 0)
    stress_neg = (not hyper) and min(hi) < 0           # the stress of total 1 then has entries of both signs
    log.update(stress_min_lower=dn_float(om_min), stress_min_rel=dn_float(om_min * len(E)),
               stress_positive=bool(stress_pos), stress_provably_negative_entry=bool(stress_neg))
    if stress_neg:
        log["stress_min_upper"] = up_float(min(hi))
    # (P)
    ok_P, valP, _ = span_exact(Z, K, d, r)
    log.update(span=bool(ok_P), span_value=up_float(valP) if valP is not None else None)
    # chi (hyperstatic): every extra contact within chi of c'
    chi = Fraction(0)
    if hyper:
        for (a, bb) in extra:
            s = F2(sum(x * y for x, y in zip(Z[a], Z[bb])), 2 * K)
            l1 = F2(sum(abs(x) for x in Z[a]) + sum(abs(y) for y in Z[bb]), K)
            chi = max(chi, abs(s - gam_c) + r * l1 + d * r * r + eps_g)
        log["chi"] = up_float(chi)
    # (Q): constants
    mu_p = res["mu_p"]
    ups = max(abs(dy(xs[a][k]) - F2(Z[a][k], K)) for a in range(nV) for k in range(d)) + r
    theta = 2 * nV * ups * mu_p
    zeta = kappa = rho = None
    sqd_up = sqrt_up(Fraction(d))
    sqnd_up = sqrt_up(Fraction(nV * d))
    if stress_pos and theta < 1:
        zeta = (1 - theta) / mu_p
        kappa = zeta * om_min / (1 - om_min)
        Lam = 1 + max(abs(c_lo), abs(c_hi)) + chi
        rho = min(kappa / (4 * Lam * sqd_up), Fraction(1, 2))
        log.update(zeta=dn_float(zeta), kappa=dn_float(kappa), kappa_half=dn_float(kappa / (2 * sqd_up)),
                   rho=dn_float(rho), Lambda=up_float(Lam))
    log.update(mu_p=up_float(mu_p), theta=up_float(theta))
    # the hyperstatic inequality and the tube
    tube = None
    if hyper and kappa is not None:
        tube = 4 * sqnd_up * chi / kappa
        log["tube"] = up_float(tube)
        log["hyper_inequality"] = bool(tube < rho)
    # (G) on V: margin for a hyperstatic graph 2 chi + 4 sqrt(d) chi / kappa (as the census)
    Eset = set(tuple(sorted(e)) for e in E)
    margin = Fraction(0)
    if hyper and kappa is not None:
        margin = 2 * chi + 4 * sqd_up * chi / kappa
    c_req = c_lo - margin
    okG, gap = gaps_exact(Z, K, d, Eset, r, c_req)
    log.update(gaps=bool(okG), gap_lower=dn_float(gap) if gap is not None else None)
    Y = keep["Y"]
    okR, bad = rattlers_exact(Y, V, Z, K, r, c_req, d)
    log.update(rattler_gaps=bool(okR), rattler_failures=[list(map(int, b_)) for b_ in bad[:5]])
    # (D)
    rat = [i for i in range(n) if i not in set(int(v) for v in V)]
    D, Dc, Dr = distance_arb(Z, K, [dec[v] for v in V], r, Yrat=[Y[i] for i in rat], dec_rat=[dec[i] for i in rat])
    log.update(dist_ls=up_float(D), dist_coord=up_float(Dc), dist_rattlers=up_float(Dr))
    if rho is not None:
        log["code_inside_rho"] = bool(D < rho)
    if hyper and tube is not None:
        log["dist_ls_plus_tube"] = up_float(D + tube)
    # verdict
    base = ok_P and okG and okR
    if not hyper:
        if base and stress_pos:
            log["status"] = "PASS"
        elif base and stress_neg:
            log["status"] = "NONOPT"
        else:
            log["status"] = "FAIL " + ",".join(x for x, f in (("stress", stress_pos), ("span", ok_P), ("gaps", okG),
                                                             ("rattlers", okR)) if not f)
    else:
        if base and stress_pos and kappa is not None and tube is not None and tube < rho:
            log["status"] = "LOCALMIN"
        else:
            log["status"] = "FAIL " + ",".join(x for x, f in (("stress", stress_pos), ("span", ok_P), ("gaps", okG),
                                                             ("rattlers", okR), ("kappa", kappa is not None),
                                                             ("tube", tube is not None and rho is not None and tube < rho))
                                               if not f)
    log["seconds"] = round(time.time() - t0, 2)
    return log


def main():
    import snapshot
    out = sys.argv[1]
    args = sys.argv[2:]
    sh, nsh = 0, 1
    only = None
    statuses = ("PASS", "LOCALMIN", "NONOPT")
    if "--shard" in args:
        sh, nsh = (int(x) for x in args[args.index("--shard") + 1].split("/"))
    if "--only" in args:
        only = set(args[args.index("--only") + 1].split(","))
    if "--status" in args:
        statuses = tuple(args[args.index("--status") + 1].split(","))
    os.makedirs(out, exist_ok=True)
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    import glob
    todo = []
    for f in sorted(glob.glob(os.path.join(CERT, "census", "results_v4", "*.json"))):
        j = json.load(open(f))
        for st in statuses:
            att = [a for a in j["attempts"] if a.get("status") == st]
            if not att:
                continue
            if st == "NONOPT" and j["status"] in ("PASS", "LOCALMIN"):
                continue
            if st in ("PASS", "LOCALMIN") and j["status"] != st:
                continue
            a = att[0]
            key = f"{j['d']}_{j['n']}"
            if only and key not in only:
                continue
            todo.append((a["unknowns"], j["d"], j["n"], a["tol"], st))
    todo = sorted(set(todo))
    for k, (U, d, n, tol, st) in enumerate(todo):
        if k % nsh != sh:
            continue
        js = os.path.join(out, f"d{d}_n{n}_{st}_{tol:.0e}.json")
        if os.path.exists(js):
            continue
        log = {"d": d, "n": n, "census_status": st, "program_sha256": SHA}
        t0, c0 = time.time(), time.process_time()
        try:
            pts = snapshot.read_points(SNAP, idx[(d, n)])
            dec = [[t.strip() for t in p] for p in pts]
            X = np.array([[float(x) for x in p] for p in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            check_entry(X, dec, d, n, tol, st, log)
        except Exception:
            log["status"] = "error"
            log["error"] = traceback.format_exc()[-2000:]
        log["wall"] = round(time.time() - t0, 2)
        log["cpu_seconds"] = round(time.process_time() - c0, 2)
        _ABS.clear()
        with open(js + ".part", "w") as fo:
            json.dump(log, fo, default=str)
        os.replace(js + ".part", js)
        print(f"d{d}_n{n} {st} {U} {log['status']} {log['wall']}s", flush=True)


if __name__ == "__main__":
    main()
