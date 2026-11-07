"""Exact verification of the classical criterion (paper, Theorem 2.5) for codes with exact coordinates: no floating-point
arithmetic enters any decision.

  python3 -I certify/census/exact_strict.py SNAPDIR RESULTS_DIR OUTDIR [--only d_n,...] [--with-optimal] [--shard I/K]

Without --only: every entry certified LOCALMIN, and with --with-optimal also every entry marked provably optimal that
the census did not certify by the isostatic certificate (some are above its size limit).

A code qualifies when its stored coordinates are integers, or when, after division by the smallest nonzero absolute
coordinate, each of them is recognised (untrusted: continued fractions or PSLQ on 16 digits) as an element a + b sqrt(D)
of Q(sqrt D), D in {2, 3, 5}, with small denominators; or else when every inner product is recognised as a rational of
denominator at most 10^4 (then check_gram verifies the rational Gram matrix, the uniform stress and bar-rigidity in the
coordinates of a basis of d points, all exactly). The exact configuration is y (one vector per point; the code is
y_i / |y_i|). Then, exactly, in Q(sqrt D):
 (i)   contacts: the pairs with the largest cosine <y_a, y_b> / (|y_a| |y_b|), every other pair strictly smaller
       (cosines compared through their signs and squares); V = the points in a contact pair, the others are rattlers;
 (ii)  span: d points of V with nonzero determinant;
 (iii) a strictly positive equilibrium stress: weights mu_ab > 0 with sum_b mu_ab y_b a multiple of y_a for every a
       (then omega_ab = mu_ab |y_a| |y_b| is a stress of the unit vectors); tried: uniform in the unit vectors (when all
       |y_a| / |y_0| are rational), uniform in y, and otherwise a rational stress found by projecting a numerical one
       onto the exact stress space (exact linear algebra over Q, for rational codes);
 (iv)  bar-rigidity: with u_a = w_a / |y_a|, the infinitesimal motions w are the solutions of <y_a, u_a> = 0 and
       N_b <y_a, u_b> + N_a <y_b, u_a> = 0 (ab in E; N = |y|^2), whose matrix M has rational (or Q(sqrt D)) entries;
       the rank of M over F_p (sqrt D mapped to a square root of D mod p) is a lower bound for its rank, and it is
       required to be |V| d - d(d-1)/2, the largest possible (the infinitesimal rotations are in the kernel).
Then the exact configuration is a strict local minimum of the largest cosine modulo rotations and rattlers.
"""
import hashlib, json, math, os, sys, time, traceback
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/flintlib"))
import numpy as np
import flint
import mpmath

SHA = {"exact_strict.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()}


def prime_for(D):
    """The largest prime p < 2^31 with p = 3 mod 4 in which D is a square (D = 1: any), and a square root of D."""
    p = (1 << 31) - 1
    while True:
        if p % 4 == 3 and flint.fmpz(p).is_prime():
            r = pow(D, (p + 1) // 4, p)
            if r * r % p == D % p:
                return p, r
        p -= 2


# ============================================================================================ the field Q(sqrt D)
class Q:
    """a + b sqrt(D) with rational a, b (D fixed per run; D = 1 means b = 0)."""
    __slots__ = ("a", "b")
    D = 1

    def __init__(self, a, b=0):
        self.a, self.b = Fraction(a), Fraction(b)

    def __add__(s, o):
        return Q(s.a + o.a, s.b + o.b)

    def __sub__(s, o):
        return Q(s.a - o.a, s.b - o.b)

    def __neg__(s):
        return Q(-s.a, -s.b)

    def __mul__(s, o):
        if isinstance(o, (int, Fraction)):
            return Q(s.a * o, s.b * o)
        return Q(s.a * o.a + Q.D * s.b * o.b, s.a * o.b + s.b * o.a)

    def iszero(s):
        return s.a == 0 and s.b == 0

    def sign(s):
        """Exact sign of a + b sqrt(D)."""
        if s.b == 0 or Q.D == 1:
            return (s.a > 0) - (s.a < 0)
        sa, sb = (s.a > 0) - (s.a < 0), (s.b > 0) - (s.b < 0)
        if sa == 0:
            return sb
        if sa == sb:
            return sa
        # opposite signs: compare a^2 with D b^2
        t = s.a * s.a - Q.D * s.b * s.b
        return sa if t > 0 else (-sa if t < 0 else 0)

    def modp(s, p, r):
        f = lambda q: q.numerator % p * pow(q.denominator, -1, p) % p
        return (f(s.a) + f(s.b) * r) % p

    def __float__(s):
        return float(s.a) + float(s.b) * math.sqrt(Q.D)

    def inv(s):
        nrm = s.a * s.a - Q.D * s.b * s.b
        return Q(s.a / nrm, -s.b / nrm)

    def __str__(s):
        return str(s.a) if s.b == 0 else f"{s.a} + {s.b} sqrt({Q.D})"


def rank_mod(A, p):
    """The rank of the int64 matrix A over F_p (p < 2^31; products stay below 2^62), by Gaussian elimination in
    chunks. A lower bound for the rank over the field of the matrix that A represents."""
    A = np.asarray(A, dtype=np.int64) % p
    rank, nr = 0, A.shape[0]
    for c in range(A.shape[1]):
        if rank == nr:
            break
        nz = np.nonzero(A[rank:, c])[0]
        if len(nz) == 0:
            continue
        pr = rank + int(nz[0])
        if pr != rank:
            A[[rank, pr]] = A[[pr, rank]]
        A[rank] = (A[rank] * pow(int(A[rank, c]), p - 2, p)) % p
        col = A[:, c].copy()
        col[rank] = 0
        idx = np.nonzero(col)[0]
        piv = A[rank][None, :]
        for k0 in range(0, len(idx), 512):
            ii = idx[k0:k0 + 512]
            A[ii] = (A[ii] - (col[ii, None] * piv) % p) % p
        rank += 1
    return rank


def rigidity_rows(nV, E, g, d, attempt=0, max_entries=5 * 10 ** 7, seed=2026):
    """The tangency rows and the contact rows to use: all of them, or (when the matrix would be larger than
    max_entries) a random subset of the contact rows, larger at each attempt (a subset of the rows can only have a
    smaller rank, so the test stays a lower bound). Returns None when no larger subset fits."""
    cols = nV * d
    if (nV + len(E)) * cols <= max_entries:
        return (list(E), False) if attempt == 0 else None
    import random
    need = cols - g - nV
    fit = max_entries // cols - nV                      # the most contact rows that fit in memory
    sizes = [min(len(E), need + 2048, fit), min(len(E), 2 * need + 2048, fit), min(len(E), fit)]
    if attempt >= 3 or sizes[attempt] < need or (attempt > 0 and sizes[attempt] <= sizes[attempt - 1]):
        return None                                     # too few rows fit to reach the rank needed, or no larger set
    return random.Random(seed + attempt).sample(list(E), sizes[attempt]), True


def rank_rows_q(V, pos, Es, d, cols, p, tang, cb, ca):
    """The rank over F_p of the rigidity matrix with the tangency rows of V and the contact rows Es: in the row of a
    contact {a, b}, cb(a, b, k) goes to block b and ca(a, b, k) to block a."""
    nV = len(V)
    Mp = np.zeros((nV + len(Es), cols), dtype=np.int64)
    for a in V:
        for k in range(d):
            Mp[pos[a], pos[a] * d + k] = tang(a, k) % p
    for e, (a, b) in enumerate(Es):
        i = nV + e
        for k in range(d):
            Mp[i, pos[b] * d + k] = (Mp[i, pos[b] * d + k] + cb(a, b, k)) % p
            Mp[i, pos[a] * d + k] = (Mp[i, pos[a] * d + k] + ca(a, b, k)) % p
    return rank_mod(Mp, p)


def dot(u, v):
    acc = Q(0)
    for x, y in zip(u, v):
        if not (x.iszero() or y.iszero()):
            acc = acc + x * y
    return acc


def cmp_cos(ip1, n1, ip2, n2):
    """Sign of ip1 / sqrt(n1) - ip2 / sqrt(n2), with n1, n2 > 0 (products of squared norms), exactly."""
    s1, s2 = ip1.sign(), ip2.sign()
    if s1 != s2:
        return (s1 > s2) - (s1 < s2)
    if s1 == 0:
        return 0
    t = (ip1 * ip1 * n2 - ip2 * ip2 * n1).sign()     # compares ip1^2 / n1 with ip2^2 / n2
    return t if s1 > 0 else -t


# ============================================================================================ UNTRUSTED: recognition
def recognise_coords(dec):
    """Exact coordinates y (lists of Q) and D, or None. Integers are taken as they are."""
    fr = [[Fraction(t) for t in row] for row in dec]
    if all(x.denominator == 1 for row in fr for x in row):
        Q.D = 1
        return [[Q(x) for x in row] for row in fr], 1
    vals = sorted({abs(float(x)) for row in fr for x in row if x != 0})
    s0 = vals[0]
    for D in (1, 5, 2, 3):
        rec = {}
        ok = True
        for v in vals:
            r = v / s0
            if D == 1:
                f = Fraction(r).limit_denominator(1000)
                if abs(float(f) - r) > 1e-12 * max(1, r):
                    ok = False
                    break
                rec[v] = (f, Fraction(0))
            else:
                mpmath.mp.dps = 15
                rel = mpmath.pslq([r, 1, mpmath.sqrt(D)], maxcoeff=1000, maxsteps=10000, tol=mpmath.mpf(10) ** -12)
                if not rel or rel[0] == 0:
                    ok = False
                    break
                a, b = Fraction(-rel[1], rel[0]), Fraction(-rel[2], rel[0])
                if abs(float(a) + float(b) * math.sqrt(D) - r) > 1e-12 * max(1, r):
                    ok = False
                    break
                rec[v] = (a, b)
        if ok:
            Q.D = D
            y = []
            for row in fr:
                yr = []
                for x in row:
                    if x == 0:
                        yr.append(Q(0))
                    else:
                        a, b = rec[abs(float(x))]
                        yr.append(Q(a, b) if x > 0 else Q(-a, -b))
                y.append(yr)
            return y, D
    return None, None


# ============================================================================================ TRUSTED: the checks
def check(y, d, log):
    n = len(y)
    N = [dot(v, v) for v in y]
    if any(x.sign() <= 0 for x in N):
        log["status"] = "a zero vector"
        return False
    # (i) contacts: a reference pair with the largest cosine (untrusted choice), then every pair exactly
    yf = np.array([[float(x) for x in v] for v in y])
    yn = yf / np.linalg.norm(yf, axis=1, keepdims=True)
    G = yn @ yn.T
    np.fill_diagonal(G, -2)
    a0, b0 = np.unravel_index(np.argmax(G), G.shape)
    ip0, n0 = dot(y[a0], y[b0]), N[a0] * N[b0]
    E = []
    for a in range(n):
        for b in range(a + 1, n):
            c = cmp_cos(dot(y[a], y[b]), N[a] * N[b], ip0, n0)
            if c > 0:
                log["status"] = "the reference pair is not a largest cosine"
                return False
            if c == 0:
                E.append((a, b))
    V = sorted({a for e in E for a in e})
    pos = {v: k for k, v in enumerate(V)}
    nV, m = len(V), len(E)
    g = d * (d - 1) // 2
    log.update(contacts=m, nV=nV, rattlers=n - nV, isostatic_count=nV * (d - 1) - g + 1,
               cos_sign=ip0.sign(), cos_squared=str(ip0 * ip0 * n0.inv()), c_float=float(ip0) / math.sqrt(float(n0)))
    # (ii) span: d points with nonzero determinant (choice untrusted; determinant exact)
    import scipy.linalg as sl
    _, _, piv = sl.qr(yf[V].T, pivoting=True, mode="economic")
    B = [V[i] for i in sorted(piv[:d])]
    Mx = [[y[b][k] for k in range(d)] for b in B]
    det_nonzero = True
    for c in range(d):
        pr = next((r for r in range(c, d) if not Mx[r][c].iszero()), None)
        if pr is None:
            det_nonzero = False
            break
        Mx[c], Mx[pr] = Mx[pr], Mx[c]
        inv = Mx[c][c].inv()
        for r in range(c + 1, d):
            if not Mx[r][c].iszero():
                f = Mx[r][c] * inv
                Mx[r] = [x - f * z for x, z in zip(Mx[r], Mx[c])]
    log["span"] = det_nonzero
    if not det_nonzero:
        log["status"] = "span failed"
        return False
    # (iii) a strictly positive stress
    nbr = {a: [] for a in V}
    for a, b in E:
        nbr[a].append(b)
        nbr[b].append(a)
    tq = []
    rational_t = True
    for a in range(n):
        r2 = None
        if Q.D == 1 or N[a].b == 0 and N[0].b == 0:
            q = N[a].a / N[V[0]].a
            s = math.isqrt(q.numerator), math.isqrt(q.denominator)
            if s[0] ** 2 == q.numerator and s[1] ** 2 == q.denominator:
                r2 = Fraction(s[0], s[1])
        tq.append(r2)
        if a in pos and r2 is None:
            rational_t = False

    def is_stress(weight):
        for a in V:
            s = [Q(0)] * d
            for b in nbr[a]:
                w = weight(a, b)
                s = [x + y[b][k] * w for k, x in enumerate(s)]
            lam_num = dot(s, y[a])
            # s parallel to y_a: N_a s - <s, y_a> y_a = 0
            for k in range(d):
                if not (s[k] * N[a] - y[a][k] * lam_num).iszero():
                    return False
        return True
    kind = None
    if rational_t and is_stress(lambda a, b: 1 / (tq[a] * tq[b])):
        kind = "uniform in the unit vectors"
    elif is_stress(lambda a, b: Fraction(1)):
        kind = "uniform in y"
    elif Q.D == 1 and m <= 4000:
        kind = rational_stress(y, V, E, nbr, d, log)
    log["stress"] = kind
    if kind is None:
        log["status"] = "no exact positive stress found"
        return False
    # (iv) bar-rigidity: rank of M over F_p (a lower bound for its rank over K)
    p, r = prime_for(Q.D)
    cols = nV * d
    rk, sampled = -1, None
    for attempt in range(3):
        rr = rigidity_rows(nV, E, g, d, attempt)
        if rr is None:
            break
        Es, sampled = rr
        rk = rank_rows_q(V, pos, Es, d, cols, p, lambda a, k: y[a][k].modp(p, r),
                         lambda a, b, k: (N[b] * y[a][k]).modp(p, r), lambda a, b, k: (N[a] * y[b][k]).modp(p, r))
        if rk == nV * d - g:
            break
    log.update(rank_mod_p=int(rk), rank_needed=nV * d - g, prime=p, rigidity_rows_sampled=sampled)
    if rk != nV * d - g:
        log["status"] = ("bar-rigidity not shown (rank mod p %d < %d)" % (rk, nV * d - g)) if rk >= 0 else \
            "rigidity matrix above the memory limit"
        return False
    log["status"] = "STRICT (exact)"
    return True


def rational_stress(y, V, E, nbr, d, log):
    """A strictly positive rational stress (D = 1). Untrusted: the stress of largest smallest entry by a linear program.
    Exact: the reduced row echelon form of [A | b] for the equilibrium equations (sum_b mu_ab y_b - lambda_a y_a = 0,
    sum mu = 1) over Q; the free unknowns are set to the rounded LP values and the pivot unknowns solved exactly.
    Accepted if A z = b holds exactly and every mu is positive."""
    from scipy.optimize import linprog
    pos = {v: k for k, v in enumerate(V)}
    nV, m = len(V), len(E)
    ncol = m + nV
    yf = np.array([[float(x) for x in v] for v in y])
    Aeq = np.zeros((nV * d + 1, ncol + 1))
    for e, (a, b) in enumerate(E):
        Aeq[pos[a] * d:(pos[a] + 1) * d, e] += yf[b]
        Aeq[pos[b] * d:(pos[b] + 1) * d, e] += yf[a]
    for a in V:
        Aeq[pos[a] * d:(pos[a] + 1) * d, m + pos[a]] = -yf[a]
    Aeq[nV * d, :m] = 1.0
    beq = np.r_[np.zeros(nV * d), 1.0]
    Aub = np.hstack([-np.eye(m), np.zeros((m, nV)), np.ones((m, 1))])
    lp = linprog(np.r_[np.zeros(ncol), -1.0], A_ub=Aub, b_ub=np.zeros(m), A_eq=Aeq, b_eq=beq,
                 bounds=[(None, None)] * ncol + [(None, 1.0)], method="highs")
    if lp.status != 0 or lp.x[-1] <= 1e-12:
        log["stress_lp"] = None if lp.status != 0 else float(lp.x[-1])
        return None
    z0 = [Fraction(float(x)).limit_denominator(10 ** 12) for x in lp.x[:ncol]]
    fq = lambda q: flint.fmpq(q.numerator, q.denominator)
    Ab = flint.fmpq_mat(nV * d + 1, ncol + 1)
    for e, (a, b) in enumerate(E):
        for k in range(d):
            Ab[pos[a] * d + k, e] = Ab[pos[a] * d + k, e] + fq(y[b][k].a)
            Ab[pos[b] * d + k, e] = Ab[pos[b] * d + k, e] + fq(y[a][k].a)
        Ab[nV * d, e] = 1
    for a in V:
        for k in range(d):
            Ab[pos[a] * d + k, m + pos[a]] = -fq(y[a][k].a)
    Ab[nV * d, ncol] = 1
    R, rk = Ab.rref()
    piv = []
    for i in range(rk):
        j = next(j for j in range(ncol + 1) if R[i, j] != 0)
        if j == ncol:
            log["stress_projection"] = "inconsistent"
            return None
        piv.append(j)
    pset = set(piv)
    z = [fq(q) for q in z0]
    for i, j in enumerate(piv):
        s = R[i, ncol]
        for jj in range(ncol):
            if jj not in pset and R[i, jj] != 0:
                s = s - R[i, jj] * z[jj]
        z[j] = s
    zm = flint.fmpq_mat(ncol, 1, z)
    A = flint.fmpq_mat(nV * d + 1, ncol, [Ab[i, j] for i in range(nV * d + 1) for j in range(ncol)])
    bv = flint.fmpq_mat(nV * d + 1, 1, [Ab[i, ncol] for i in range(nV * d + 1)])
    if A * zm != bv:
        log["stress_projection"] = "not exact"
        return None
    mu = z[:m]
    if min(mu) <= 0:
        log["stress_projection"] = "not positive"
        return None
    mn = min(mu)
    log["stress_min_exact_rel"] = int(mn.p) / int(mn.q) * m
    return "rational (exact solution near the linear program's)"


def recognise_gram(X):
    """UNTRUSTED. A rational Gram matrix: every inner product of the code (normalised in double precision) rounded to a
    rational with denominator at most 10^4, accepted if within 1e-11. Returns the matrix (Fractions) or None."""
    G = X @ X.T
    n = len(X)
    out = [[Fraction(1) if i == j else None for j in range(n)] for i in range(n)]
    cache = {}
    dev = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            v = float(G[i, j])
            key = round(v, 9)
            if key not in cache:
                q = Fraction(v).limit_denominator(10 ** 4)
                if abs(float(q) - v) > 1e-11:
                    return None, None
                cache[key] = q
            q = cache[key]
            dev = max(dev, abs(float(q) - v))
            out[i][j] = out[j][i] = q
    return out, dev


def check_gram(G, X, d, log, max_entries=6 * 10 ** 7):
    """TRUSTED (exact, rational). G is the Gram matrix of a configuration P* in R^d with the contacts E (pairs at the
    largest entry c), the others below c; P* has a strictly positive stress (uniform, checked in the coordinates of a
    basis of d points) and is bar-rigid (rank over F_q of the rigidity matrix in those coordinates)."""
    n = len(G)
    c = max(G[i][j] for i in range(n) for j in range(i + 1, n))
    E = [(i, j) for i in range(n) for j in range(i + 1, n) if G[i][j] == c]
    V = sorted({a for e in E for a in e})
    pos = {v: k for k, v in enumerate(V)}
    nV, m, g = len(V), len(E), d * (d - 1) // 2
    log.update(contacts=m, nV=nV, rattlers=n - nV, cos=str(c), isostatic_count=nV * (d - 1) - g + 1)
    import scipy.linalg as sl
    _, _, piv = sl.qr(X[V].T, pivoting=True, mode="economic")
    B = [V[i] for i in sorted(piv[:d])]
    fq = lambda q: flint.fmpq(q.numerator, q.denominator)
    GBB = flint.fmpq_mat(d, d, [fq(G[a][b]) for a in B for b in B])
    # G_BB positive definite: the pivots of its LDL^T (exact) are positive
    M = [[G[a][b] for b in B] for a in B]
    for k in range(d):
        if M[k][k] <= 0:
            log["status"] = "G_BB not positive definite"
            return False
        for i in range(k + 1, d):
            f = M[i][k] / M[k][k]
            for j in range(k, d):
                M[i][j] -= f * M[k][j]
    GBV = flint.fmpq_mat(d, nV, [fq(G[b][v]) for b in B for v in V])
    W = GBB.solve(GBV)                                 # exact: p_v = sum_k W[k, v] p_{B_k}
    GVV = flint.fmpq_mat(nV, nV, [fq(G[a][b]) for a in V for b in V])
    if GBV.transpose() * W != GVV:
        log["status"] = "the rational Gram matrix does not have rank d"
        return False
    # every point outside V (a rattler) is in R^d too: its Gram row must be G_{r,B} W
    rat = [r for r in range(n) if r not in pos]
    if rat:
        GrB = flint.fmpq_mat(len(rat), d, [fq(G[r][b]) for r in rat for b in B])
        Wall = GBB.solve(flint.fmpq_mat(d, n, [fq(G[b][v]) for b in B for v in range(n)]))
        if GrB * Wall != flint.fmpq_mat(len(rat), n, [fq(G[r][v]) for r in rat for v in range(n)]):
            log["status"] = "the rational Gram matrix does not have rank d (rattlers)"
            return False
    # the uniform stress: sum of the neighbours of a is a multiple of a, in the coordinates W
    nbr = {a: [] for a in V}
    for a, b in E:
        nbr[a].append(b)
        nbr[b].append(a)
    for a in V:
        s = [sum((W[k, pos[b]] for b in nbr[a]), flint.fmpq(0)) for k in range(d)]
        w = [W[k, pos[a]] for k in range(d)]
        k0 = next(k for k in range(d) if w[k] != 0)
        lam = s[k0] / w[k0]
        if any(s[k] != lam * w[k] for k in range(d)):
            log["status"] = "the uniform stress is not in equilibrium"
            return False
    log["stress"] = "uniform"
    # bar-rigidity: w_a = sum_k u_ak p_{B_k}; rows <p_a, w_a> and <p_a, w_b> + <p_b, w_a>; rank over F_q
    p, _ = prime_for(1)
    cols = nV * d
    gm = {}
    for a in V:
        for k, b in enumerate(B):
            q = G[a][b]
            gm[(a, k)] = q.numerator % p * pow(q.denominator, -1, p) % p
    rk, sampled = -1, None
    for attempt in range(3):
        rr = rigidity_rows(nV, E, g, d, attempt)
        if rr is None:
            break
        Es, sampled = rr
        rk = rank_rows_q(V, pos, Es, d, cols, p, lambda a, k: gm[(a, k)], lambda a, b, k: gm[(a, k)],
                         lambda a, b, k: gm[(b, k)])
        if rk == cols - g:
            break
    log.update(rank_mod_p=int(rk), rank_needed=cols - g, prime=p, rigidity_rows_sampled=sampled)
    if rk != cols - g:
        log["status"] = ("bar-rigidity not shown (rank mod p %d < %d)" % (rk, cols - g)) if rk >= 0 else \
            "rigidity matrix above the memory limit"
        return False
    log["status"] = "STRICT (exact Gram)"
    return True


def main():
    snapdir, resdir, out = sys.argv[1], sys.argv[2], sys.argv[3]
    args = sys.argv[4:]
    only = set(args[args.index("--only") + 1].split(",")) if "--only" in args else None
    import glob
    import snapshot
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    todo = []
    for f in sorted(glob.glob(os.path.join(resdir, "*.json"))):
        j = json.load(open(f))
        key = f"{j['d']}_{j['n']}"
        optimal_open = "--with-optimal" in args and idx[(j["d"], j["n"])]["optimal"] and j.get("status") != "PASS"
        if (only and key in only) or (not only and (j.get("status") == "LOCALMIN" or optimal_open)):
            todo.append((j["d"], j["n"]))
    os.makedirs(out, exist_ok=True)
    sh, nsh = (int(x) for x in args[args.index("--shard") + 1].split("/")) if "--shard" in args else (0, 1)
    for t_, (d, n) in enumerate(todo):
        if t_ % nsh != sh:
            continue
        js = os.path.join(out, f"d{d}_n{n}.json")
        if os.path.exists(js):
            continue
        log = {"d": d, "n": n, "code_sha256": SHA}
        t0 = time.process_time()
        try:
            pts = snapshot.read_points(snapdir, idx[(d, n)])
            dec = [[t.strip() for t in p] for p in pts]
            y, D = recognise_coords(dec)
            if y is None:
                X = np.array([[float(t) for t in row] for row in dec])
                X = X / np.linalg.norm(X, axis=1, keepdims=True)
                G, dev = recognise_gram(X)
                if G is None:
                    log["status"] = "coordinates not recognised as exact"
                else:
                    log["field"] = "Q (Gram matrix)"
                    log["gram_deviation_from_stored"] = dev
                    check_gram(G, X, d, log)
            else:
                log["field"] = "Q" if D == 1 else f"Q(sqrt {D})"
                log["stored_integers"] = all(Fraction(t).denominator == 1 for row in dec for t in row)
                check(y, d, log)
        except Exception:
            log["status"] = "error"
            log["error"] = traceback.format_exc()[-1500:]
        log["cpu_seconds"] = round(time.process_time() - t0, 2)
        json.dump(log, open(js, "w"), default=str)
        if log["status"] != "coordinates not recognised as exact":
            print(f"d{d}_n{n} {log.get('field')} {log['status']} {log['cpu_seconds']}s", flush=True)


if __name__ == "__main__":
    main()
