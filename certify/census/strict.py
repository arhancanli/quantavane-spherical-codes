"""Strictness for hyperstatic contact graphs (paper, Section 3.5): an exact equal-contact configuration from an exact
Gram matrix, then the first-order criterion at that configuration.

  python3 -I certify/census/strict.py SNAPDIR RESULTS_DIR OUTDIR [--shard I/K] [--only d_n,...]

For an entry that the census certified LOCALMIN at the tolerance tau, with contacts E on the non-rattlers V:
 1. UNTRUSTED. The hyperstatic certificate is re-run (cq.certify) to obtain the zero P' of the square system of E' on
    the grid 2^-100, and P' is refined by chord steps with exact residuals to the grid 2^-K2. The common cosine is
    recognised as a root of an integer polynomial p: the table's minimal polynomial when it lists one, otherwise one
    found by integer relation detection. Every distinct inner product of P' is recognised as an element of
    Q[x]/(p) by PSLQ.
 2. TRUSTED, exact. p is irreducible over Q (FLINT), and c* is the unique real root of p in an isolating ball. The
    recognised symmetric matrix G over Q(c*) (unit diagonal) has G_ab = c* exactly on E, G_ab < c* on the other pairs
    of V (signs by ball arithmetic at c*), and G = G_{V,B} G_BB^{-1} G_{B,V} exactly for a set B of d indices with
    G_BB positive definite. So G is the Gram matrix of a configuration P* in (S^{d-1})^V with E(P*) = E, M(P*) = c*.
 3. TRUSTED, ball arithmetic and the kernel of cq.py. P* (the Cholesky factor of G) is enclosed coordinatewise. At
    every configuration in the enclosure: the points span R^d (span_check); the Jacobian of the square system of E'
    is nonsingular, and the combined stress of E (coefficients beta_f on E \\ E', chosen by untrusted code) is
    enclosed and strictly positive (transposed_enclosure). So (P*, E) is bar-rigid and carries a strictly positive
    equilibrium stress.
 4. TRUSTED. The rattlers, placed by untrusted code in the frame of P*, have every inner product below c*.
Then P*, with the rattlers, satisfies Theorem 2.5 and is a strict local minimum of the largest cosine modulo rotations
and rattlers, with largest cosine exactly c*. Its distance from the stored code (exactly normalised), modulo rotations,
is bounded with an exactly orthogonal rational matrix (a Cayley transform) in ball arithmetic.
"""
import hashlib, json, math, os, sys, time, traceback, warnings
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/flintlib"))
warnings.filterwarnings("ignore", message=".*encountered in matmul")
import numpy as np
import cq
import sphcert as S
from sphcert import up, frac_up, frac_down
import flint
import mpmath

SHA = dict(cq.CODE_SHA256, **{"strict.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()})
K2 = 400                      # bits of the refined P'
DPS = 110                     # decimal digits used in the integer relations
MAX_VALUES = 4000             # recognition is not attempted for more distinct inner products than this
MAX_DEGREE = 12               # nor for a cosine of larger degree
MAX_WORK = 6 * 10 ** 6        # nor when |V|^2 times the degree exceeds this (memory of the exact rank test)
PREC = 480                    # bits of the ball arithmetic


# ============================================================================================ UNTRUSTED: recognition
def parse_poly(s):
    """The table's polynomial, e.g. '7*x^3+x^2-3*x-1', as integer coefficients (constant first)."""
    s = s.replace(" ", "").replace("-", "+-")
    coef = {}
    for term in s.split("+"):
        if not term:
            continue
        if "x" in term:
            c, _, e = term.partition("x")
            c = c.rstrip("*")
            c = -1 if c == "-" else (1 if c == "" else int(c))
            e = int(e[1:]) if e.startswith("^") else 1
        else:
            c, e = int(term), 0
        coef[e] = coef.get(e, 0) + c
    return [coef.get(k, 0) for k in range(max(coef) + 1)]


def find_poly(c, maxdeg=10, maxcoeff=10 ** 8):
    """UNTRUSTED. An integer polynomial with root c (mpmath findpoly), lowest degree first; a relation is kept only if
    its coefficients use at most DPS - 25 digits in all (shorter relations than the precision, so not chance ones)."""
    for deg in range(1, maxdeg + 1):
        p = mpmath.findpoly(c, deg, maxcoeff=maxcoeff, maxsteps=50000)
        if p:
            size = sum(math.log10(abs(a) + 1) for a in p)
            return [int(a) for a in reversed(p)] if size < DPS - 25 else None
    return None


def recognise(v, cpows, maxcoeff=10 ** 12):
    """UNTRUSTED. Rationals a_0..a_{D-1} with v = sum a_k c^k, by PSLQ on (v, 1, c, ..., c^{D-1})."""
    if abs(v) < mpmath.mpf(10) ** (-(DPS - 8)):
        return [Fraction(0)] * len(cpows)
    rel = mpmath.pslq([v] + cpows, maxcoeff=maxcoeff, maxsteps=200000)
    if not rel or rel[0] == 0:
        return None
    return [Fraction(-a, rel[0]) for a in rel[1:]]


# ============================================================================================ TRUSTED: the field
class Field:
    """Q[x]/(p) for an irreducible p, elements as fmpq_poly of degree < D reduced mod p, with a real embedding x -> c*
    given by an isolating ball."""

    def __init__(self, coeffs, c_float):
        self.p = flint.fmpq_poly(coeffs)
        fac = self.p.factor()[1]
        if len(fac) != 1 or fac[0][1] != 1:
            raise ValueError("polynomial not irreducible")
        self.D = self.p.degree()
        old = flint.ctx.prec
        flint.ctx.prec = PREC
        roots = self.p.numer().complex_roots()
        flint.ctx.prec = old
        # FLINT isolates every root in a ball and returns the real roots with imaginary part exactly zero
        real = [r for r, mlt in roots if r.imag == 0 and abs(float(r.real) - c_float) < 1e-9]
        if len(real) != 1:
            raise ValueError("no unique real root of p near the certified value")
        self.c = real[0].real
        self.x = flint.fmpq_poly([0, 1]) % self.p

    def elem(self, coeffs):
        return flint.fmpq_poly([flint.fmpq(q.numerator, q.denominator) for q in coeffs]) % self.p

    def value(self, a):
        """A ball that contains a(c*)."""
        old = flint.ctx.prec
        flint.ctx.prec = PREC
        acc = flint.arb(0)
        for co in reversed(a.coeffs()):
            acc = acc * self.c + flint.arb(co)
        flint.ctx.prec = old
        return acc

    def mul(self, a, b):
        return (a * b) % self.p

    def inv(self, a):
        g, s, _ = a.xgcd(self.p)
        if g.degree() != 0:
            raise ZeroDivisionError
        return (s / g) % self.p


def coeff_mats(F, A, shape):
    """The matrix A (dict (i, j) -> element, or a list of lists) as D rational matrices (power basis)."""
    mats = [flint.fmpq_mat(shape[0], shape[1]) for _ in range(F.D)]
    for i in range(shape[0]):
        for j in range(shape[1]):
            co = A[i][j].coeffs()
            for k, q in enumerate(co):
                mats[k][i, j] = q
    return mats


def field_matmul(F, A, B):
    """TRUSTED. A B over Q[x]/(p) for coefficient-matrix lists A (D mats m x k) and B (D mats k x n): sum over r, s of
    A_r B_s x^(r+s), with x^t reduced by p."""
    D = F.D
    red = []
    for t in range(2 * D - 1):
        red.append([(flint.fmpq_poly([0] * t + [1]) % F.p).coeffs()])
    m, n = A[0].nrows(), B[0].ncols()
    out = [flint.fmpq_mat(m, n) for _ in range(D)]
    for r in range(D):
        for s in range(D):
            P = A[r] * B[s]
            co = red[r + s][0]
            for k, q in enumerate(co):
                if q != 0:
                    out[k] = out[k] + P * q
    return out


# ============================================================================================ the stage
def gram_ints(Z, pairs=None):
    n = len(Z)
    G = {}
    for a in range(n):
        za = Z[a]
        for b in range(a + 1, n):
            G[(a, b)] = sum(x * y for x, y in zip(za, Z[b]))
    return G


def strict_entry(X, dec, d, n, tol, poly_str, log):
    t0 = time.time()
    keep = {}
    cl = cq.certify(X, d, n, tol, dec=dec, keep=keep)
    log["census_status_rerun"] = cl.get("status")
    if cl.get("status") != "LOCALMIN":
        log["status"] = "not LOCALMIN on re-run"
        return log
    V, E, Esq, extra = keep["V"], keep["E"], keep["Esq"], keep["extra"]
    nV, m = len(V), len(E)
    log.update(nV=int(nV), contacts=int(m), extra=len(extra), rattlers=int(n - nV))
    # ---- 1. refine P' (untrusted)
    sysm, R = keep["sysm"], keep["R"]
    sh = K2 - keep["K"]
    Z = [[z << sh for z in row] for row in keep["Z"]]
    C = keep["C"] << sh
    ex = S.Exact(sysm, K2)
    Z, C, res = cq.refine(ex, Z, C, R, K2, iters=8)
    log["refined_residual"] = float(res)
    del keep["R"]
    mpmath.mp.dps = DPS + 10
    cmp = mpmath.mpf(C) / mpmath.mpf(2) ** K2
    # ---- distinct inner products of P' and their recognition (untrusted)
    Gi = gram_ints(Z)
    shift = 2 * K2 - 200
    clusters = {}
    for key, g in Gi.items():
        clusters.setdefault((g + (1 << (shift - 1))) >> shift, []).append(key)
    log["distinct_inner_products"] = len(clusters)
    if len(clusters) > MAX_VALUES:
        log["status"] = "too many distinct inner products"
        return log
    # ---- the polynomial
    if poly_str:
        coeffs, src = parse_poly(poly_str), "table"
    else:
        coeffs, src = find_poly(cmp), "integer relation"
    if coeffs is None:
        log["status"] = "no polynomial for the cosine"
        return log
    log.update(poly=coeffs, poly_source=src)
    try:
        F = Field(coeffs, float(cmp))
    except ValueError as e:
        # a reducible relation: keep the factor that vanishes at c
        if src == "integer relation":
            fac = flint.fmpq_poly(coeffs).factor()[1]
            best = min(fac, key=lambda fm: abs(float(mpmath.polyval([mpmath.mpf(int(q.p)) / int(q.q) for q in reversed(fm[0].coeffs())], cmp))))
            coeffs = [Fraction(int(q.p), int(q.q)) for q in best[0].coeffs()]
            den = 1
            for q in coeffs:
                den = den * q.denominator // math.gcd(den, q.denominator)
            coeffs = [int(q * den) for q in coeffs]
            log["poly"] = coeffs
            try:
                F = Field(coeffs, float(cmp))
            except ValueError as e2:
                log["status"] = "polynomial: " + str(e2)
                return log
        else:
            log["status"] = "polynomial: " + str(e)
            return log
    log["field_degree"] = F.D
    if F.D > MAX_DEGREE:
        log["status"] = "cosine of degree above %d" % MAX_DEGREE
        return log
    if nV * nV * F.D > MAX_WORK:
        log["status"] = "exact rank test above the work limit"
        return log
    mpmath.mp.dps = DPS
    cpows = [mpmath.mpf(1)]
    for _ in range(F.D - 1):
        cpows.append(cpows[-1] * cmp)
    elem_of = {}
    for ck, members in clusters.items():
        a, b = members[0]
        v = mpmath.mpf(Gi[(a, b)]) / mpmath.mpf(2) ** (2 * K2)
        co = recognise(v, cpows)
        if co is None:
            log["status"] = "an inner product not recognised in Q(c)"
            return log
        rec = sum(mpmath.mpf(q.numerator) / q.denominator * cp for q, cp in zip(co, cpows))
        if abs(rec - v) > mpmath.mpf(10) ** (-(DPS - 8)):
            log["status"] = "an inner product not recognised in Q(c)"
            return log
        elem_of[ck] = F.elem(co)
    # ---- 2. exact verification (trusted)
    one = F.elem([Fraction(1)])
    Gx = [[None] * nV for _ in range(nV)]
    for ck, members in clusters.items():
        el = elem_of[ck]
        for a, b in members:
            Gx[a][b] = Gx[b][a] = el
    for a in range(nV):
        Gx[a][a] = one
    Eset = set(E)
    for (a, b) in E:
        if Gx[a][b] != F.x:
            log["status"] = "a contact is not exactly c* in the recognised Gram matrix"
            return log
    vals = {}
    worst = None
    for ck, el in elem_of.items():
        diff = F.value(F.x - el) if el != F.x else None
        vals[ck] = diff
    for ck, members in clusters.items():
        if elem_of[ck] == F.x:
            if any(tuple(mb) not in Eset for mb in members):
                log["status"] = "a pair outside E has inner product c*"
                return log
            continue
        diff = vals[ck]
        if not diff > 0:
            log["status"] = "a pair outside E not provably below c*"
            return log
        lo = float(diff.lower())
        worst = lo if worst is None else min(worst, lo)
    log["gap_lower"] = worst
    # B: d well-conditioned points (untrusted choice)
    import scipy.linalg as sl
    Xf0 = np.array([[z / (1 << K2) for z in row] for row in Z])
    _, _, piv = sl.qr(Xf0.T, pivoting=True, mode="economic")
    B = sorted(int(i) for i in piv[:d])
    GBB = [[Gx[i][j] for j in B] for i in B]
    # G_BB positive definite: Cholesky in ball arithmetic with positive pivots
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    Gv = {}
    for ck, el in elem_of.items():
        Gv[ck] = F.value(el)
    def gval(i, j):
        if i == j:
            return flint.arb(1)
        return F.value(Gx[i][j]) if (i, j) not in Gv_cache else Gv_cache[(i, j)]
    Gv_cache = {}
    for ck, members in clusters.items():
        for a, b in members:
            Gv_cache[(a, b)] = Gv[ck]
            Gv_cache[(b, a)] = Gv[ck]
    Lc = [[flint.arb(0)] * d for _ in range(d)]
    for i in range(d):
        for j in range(i + 1):
            s_ = gval(B[i], B[j])
            for k in range(j):
                s_ = s_ - Lc[i][k] * Lc[j][k]
            if i == j:
                if not s_ > 0:
                    flint.ctx.prec = old
                    log["status"] = "G_BB not provably positive definite"
                    return log
                Lc[i][i] = s_.sqrt()
            else:
                Lc[i][j] = s_ / Lc[j][j]
    flint.ctx.prec = old
    # G = G_{V,B} G_BB^{-1} G_{B,V}: W = G_BB^{-1} G_{B,V} over the field, then compare with G
    Ainv = [[F.elem([Fraction(int(i == j))]) for j in range(d)] for i in range(d)]
    M_ = [row[:] for row in GBB]
    for c_ in range(d):
        pv = next((r for r in range(c_, d) if M_[r][c_] != 0), None)
        if pv is None:
            log["status"] = "G_BB singular"
            return log
        M_[c_], M_[pv] = M_[pv], M_[c_]
        Ainv[c_], Ainv[pv] = Ainv[pv], Ainv[c_]
        iv = F.inv(M_[c_][c_])
        M_[c_] = [F.mul(iv, e) for e in M_[c_]]
        Ainv[c_] = [F.mul(iv, e) for e in Ainv[c_]]
        for r in range(d):
            if r != c_ and M_[r][c_] != 0:
                f = M_[r][c_]
                M_[r] = [(e - F.mul(f, g)) % F.p for e, g in zip(M_[r], M_[c_])]
                Ainv[r] = [(e - F.mul(f, g)) % F.p for e, g in zip(Ainv[r], Ainv[c_])]
    GBV = [[Gx[b][j] for j in range(nV)] for b in B]
    GVB = [[Gx[i][b] for b in B] for i in range(nV)]
    W = field_matmul(F, coeff_mats(F, Ainv, (d, d)), coeff_mats(F, GBV, (d, nV)))
    P = field_matmul(F, coeff_mats(F, GVB, (nV, d)), W)
    del W
    zero = flint.fmpq(0)
    for k in range(F.D):                      # every entry of G_{V,B} G_BB^-1 G_{B,V} equals the entry of G
        ent = P[k].entries()
        for i in range(nV):
            row = Gx[i]
            base = i * nV
            for j in range(nV):
                co = row[j].coeffs()
                if ent[base + j] != (co[k] if k < len(co) else zero):
                    log["status"] = "the recognised Gram matrix does not have rank d"
                    return log
        P[k] = None
    log["exact_gram"] = True
    # ---- 3. enclosure of P* (Cholesky gauge) and the interval tests at P*
    flint.ctx.prec = PREC
    Yb = []
    for i in range(nV):
        y = [flint.arb(0)] * d
        for k in range(d):
            s_ = gval(i, B[k])
            for j in range(k):
                s_ = s_ - Lc[k][j] * y[j]
            y[k] = s_ / Lc[k][k]
        Yb.append(y)
    flint.ctx.prec = old
    Xf = np.empty((nV, d))
    rad = np.empty((nV, d))
    flint.ctx.prec = PREC
    for i in range(nV):
        for k in range(d):
            f = float(Yb[i][k])
            Xf[i, k] = f
            rad[i, k] = math.nextafter(float(abs(Yb[i][k] - flint.arb(f)).upper()), math.inf)
    flint.ctx.prec = old
    rad = up(rad)
    log["enclosure_radius"] = float(rad.max())
    iso = nV * (d - 1) - d * (d - 1) // 2 + 1
    full = S.System(Xf, E, d)
    sub = cq.choose_square_subset_qr(full, Xf, d, nV, iso)
    Esq = [E[e] for e in sub]
    ss = set(sub)
    extra = [E[e] for e in range(m) if e not in ss]
    sysm = S.System(Xf, Esq, d)
    N = sysm.nu
    Jc = sysm.jac(Xf)
    try:
        R = cq.approx_inverse(Jc, np.float64 if N <= 12000 else np.float32)
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
    # ---- 4. rattlers in the frame of P*: the stored rattlers rotated by the Procrustes map (untrusted), then placed
    rat = np.setdiff1d(np.arange(n), V)
    U_, _, Vt_ = np.linalg.svd(X[V].T @ Xf)
    O = U_ @ Vt_                                           # X_V O ~ Xf
    Ycfg = np.empty((n, d))
    Ycfg[V] = Xf
    Ycfg[rat] = X[rat] @ O
    c_lo = Fraction(math.nextafter(float(F.c.lower()), -math.inf))     # a rational lower bound of c*
    log["c_star_lower"] = float(c_lo)
    if len(rat):
        Yr, rinfo = S.place_rattlers(Ycfg, V, Xf, float(c_lo), d)
        log["rattler_placement"] = rinfo
        ok_r, rbad = rattler_gaps(Yr, V, rat, Xf, rad, c_lo, d)
        log["rattler_gaps_ok"] = ok_r
        if not ok_r:
            log["status"] = "rattler gaps failed"
            log["rattler_gap_failures"] = rbad[:5]
            return log
        Ycfg = Yr
    # ---- the distance from the stored code, modulo rotations (exactly orthogonal rational Q)
    try:
        log.update(distance_to_stored(dec, V, Yb, O, d))
    except Exception:
        log["distance_error"] = traceback.format_exc()[-300:]
    log["status"] = "STRICT"
    log["c_star"] = mpmath.nstr(cmp, 40)
    log["seconds"] = round(time.time() - t0, 2)
    return log


def rattler_gaps(Y, V, rat, Xf, rad, c_lo, d):
    """TRUSTED. Every pair (r, j), r a rattler, has <q_r, q_j> < c_lo <= c*, with q_r = Y_r / |Y_r| and q_j within rad of
    Xf_j (j in V) or q_j = Y_j / |Y_j| (j a rattler): decided exactly in rational arithmetic (squares compared)."""
    inV = {int(v): a for a, v in enumerate(V)}
    bad = []
    rs = [int(r) for r in rat]
    for ii, r in enumerate(rs):
        yr = [Fraction(float(x)) for x in Y[r]]
        nr2 = sum(x * x for x in yr)
        for j in range(len(Y)):
            if j == r or (j in rs and rs.index(j) < ii):
                continue
            if j in inV:
                a = inV[j]
                lhs = sum(x * Fraction(float(p)) for x, p in zip(yr, Xf[a])) + \
                    sum(abs(x) * Fraction(float(e)) for x, e in zip(yr, rad[a]))
                rhs2 = c_lo * c_lo * nr2                          # |q_j| = 1 exactly
            else:
                ys = [Fraction(float(x)) for x in Y[j]]
                lhs = sum(x * y for x, y in zip(yr, ys))
                rhs2 = c_lo * c_lo * nr2 * sum(x * x for x in ys)
            # need lhs < c_lo |Y_r| (times |Y_j| for a rattler j)
            if c_lo > 0:
                ok = lhs <= 0 or lhs * lhs < rhs2
            else:
                ok = lhs < 0 and lhs * lhs > rhs2
            if not ok:
                bad.append((r, j))
    return not bad, bad


def rational_orthogonal(Ot, d, bits=60):
    """An exactly orthogonal rational matrix close to the orthogonal matrix Ot: the product of the Householder
    reflections I - 2 v v^T / (v^T v) of a QR reduction of Ot, with every v rounded to the grid 2^-bits (each factor
    is exactly orthogonal whatever the rounding), times the diagonal of signs that remains."""
    A = np.array(Ot, dtype=float)
    Hs = []
    for k in range(d):
        x = A[k:, k].copy()
        v = x.copy()
        v[0] += (1.0 if x[0] >= 0 else -1.0) * np.linalg.norm(x)
        vq = [Fraction(0)] * k + [Fraction(round(t * 2 ** bits), 2 ** bits) for t in v]
        vv = sum(t * t for t in vq)
        if vv == 0:
            continue
        Hs.append((vq, vv))
        vf = np.array([float(t) for t in vq])
        A = A - 2.0 * np.outer(vf, vf @ A) / float(vv)
    one = lambda i, j: flint.fmpq(int(i == j))
    Q = flint.fmpq_mat(d, d, [one(i, j) for i in range(d) for j in range(d)])
    for vq, vv in Hs:                    # Q = H_1 H_2 ... H_k
        H = flint.fmpq_mat(d, d, [flint.fmpq(int(i == j)) - flint.fmpq(2) * flint.fmpq((vq[i] * vq[j] / vv).numerator, (vq[i] * vq[j] / vv).denominator) for i in range(d) for j in range(d)])
        Q = Q * H
    s = [1 if A[i, i] >= 0 else -1 for i in range(d)]
    S_ = flint.fmpq_mat(d, d, [s[i] if i == j else 0 for i in range(d) for j in range(d)])
    return Q * S_


def distance_to_stored(dec, V, Yb, O, d):
    """The least-squares and largest coordinate distance between the exactly normalised stored code (points of V)
    and P* Q, for an exactly orthogonal rational Q close to O^T (rows: points; Xf O^T ~ X_V). The configurations
    P* Q (Q in O(d)) all have the same inner products as P*, so each is an exact strict local minimum as well."""
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    Q = rational_orthogonal(O.T, d)
    ok = Q.transpose() * Q == flint.fmpq_mat(d, d, [int(i == j) for i in range(d) for j in range(d)])
    Qa = [[flint.arb(Q[j, k]) for k in range(d)] for j in range(d)]
    tot = flint.arb(0)
    mx = 0.0
    for a, v in enumerate(V):
        y = Yb[a]
        q = [Fraction(t) for t in dec[v]]
        s2 = sum(x * x for x in q)
        nrm = flint.arb(flint.fmpq(s2.numerator, s2.denominator)).sqrt()
        for k in range(d):
            z = flint.arb(0)
            for j in range(d):
                z = z + y[j] * Qa[j][k]
            diff = z - flint.arb(flint.fmpq(q[k].numerator, q[k].denominator)) / nrm
            tot = tot + diff * diff
            mx = max(mx, math.nextafter(float(abs(diff).upper()), math.inf))
    D2 = tot.sqrt()
    flint.ctx.prec = old
    return {"rotation_exactly_orthogonal": bool(ok), "dist_ls_upper": math.nextafter(float(D2.upper()), math.inf),
            "dist_coord_upper": mx}


def main():
    snapdir, resdir, out = sys.argv[1], sys.argv[2], sys.argv[3]
    args = sys.argv[4:]
    sh, nsh = 0, 1
    only = None
    if "--shard" in args:
        sh, nsh = (int(x) for x in args[args.index("--shard") + 1].split("/"))
    if "--only" in args:
        only = set(args[args.index("--only") + 1].split(","))
    import glob
    import snapshot
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    todo = []
    for f in sorted(glob.glob(os.path.join(resdir, "*.json"))):
        j = json.load(open(f))
        if j.get("status") != "LOCALMIN":
            continue
        if only and f"{j['d']}_{j['n']}" not in only:
            continue
        a = [x for x in j["attempts"] if x.get("status") == "LOCALMIN"][0]
        todo.append((a["unknowns"], j["d"], j["n"], a["tol"]))
    todo.sort()
    os.makedirs(out, exist_ok=True)
    for k, (size, d, n, tol) in enumerate(todo):
        if k % nsh != sh:
            continue
        js = os.path.join(out, f"d{d}_n{n}.json")
        if os.path.exists(js):
            continue
        log = {"d": d, "n": n, "tol": tol, "unknowns": size, "code_sha256": SHA}
        t0, c0 = time.time(), time.process_time()
        try:
            pts = snapshot.read_points(snapdir, idx[(d, n)])
            dec = [[t.strip() for t in p] for p in pts]
            X = np.array([[float(x) for x in p] for p in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            strict_entry(X, dec, d, n, tol, idx[(d, n)]["poly"], log)
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
