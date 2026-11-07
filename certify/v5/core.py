"""The exact contact core of a table entry (paper, version 5): an exact Gram matrix for the core points, the stored
rattlers kept where they are, and the certified gaps and distances (used by subspace.py and rigidity.py).

UNTRUSTED: which points form the core (pairs within tol of the largest cosine of the code in double precision, then
the points with fewer than d such pairs removed repeatedly); the recognition of the core's inner products in
K = Q(sqrt DISC); the basis B; the least-squares rotation.

TRUSTED (exact over K, and FLINT ball arithmetic with rigorous radii):
 (G) the recognised core matrix G (indexed by the core) is symmetric with unit diagonal, G_BB is positive definite for a
     set B of r core points (exact LDL^T pivots positive), and G = G_{.B} G_BB^{-1} G_{B.} exactly. So there is a core
     configuration P* (unique up to O(d)) with Gram matrix G spanning an r-dimensional subspace, r <= d.
 (E) c = the largest off-diagonal entry of G; E = {ab : G_ab = c}; every other core pair is strictly below c; every core
     point has at least d contacts in E (so the core is closed under the removal rule).
 (Q) an exactly orthogonal rational d x d matrix Q (product of Householder reflections with rational vectors) and ball
     enclosures of the rows of P*Q (P* = [G_{.B} L^{-T}, 0], L the Cholesky factor of G_BB, in ball arithmetic).
 (R) every pair of a core point of P*Q and a stored rattler (normalised exactly), and every pair of stored rattlers, has
     inner product strictly below c (ball arithmetic; rattler pairs from exact rationals). So in the configuration
     Q* = (P*Q on the core, stored rattlers) the contacts are exactly E, the largest cosine is exactly c, and the
     rattlers touch nothing.
 (D) the least-squares and the largest coordinate distance between Q* and the exactly normalised stored code.
"""
import math, os, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/flintlib"))
import numpy as np
import flint
from kfield import K, ZERO, ONE, ksum, fq, is_pos_def, solve as ksolve
import exact2 as X2
import dist

PREC = 300


def propose_core(Xf, d, tol=1e-10):
    """UNTRUSTED: the core at tolerance tol (indices into the code)."""
    n = len(Xf)
    G = Xf @ Xf.T
    iu = np.triu_indices(n, 1)
    c = G[iu].max()
    sel = G[iu] >= c - tol
    E = list(zip(iu[0][sel].tolist(), iu[1][sel].tolist()))
    V = set(range(n))
    while True:
        deg = dict.fromkeys(V, 0)
        for a, b in E:
            if a in V and b in V:
                deg[a] += 1
                deg[b] += 1
        bad = {i for i in V if deg[i] < d}
        if not bad:
            break
        V -= bad
    return sorted(V)


def recognise_core(Xf, V, disc):
    """UNTRUSTED: the inner products of the core points, recognised in K (exact2.recognise_value)."""
    XV = Xf[V]
    G, dev, nval = X2.recognise_gram(XV, disc)
    return G, dev, nval


def gram_rank_basis(XV, G, d):
    """UNTRUSTED choice of B (pivoted QR on the core points), then the TRUSTED check (G) with r = |B|. The rank r is
    taken as the numerical rank; the exact check then proves rank exactly r (G_BB positive definite, G = G_.B G_BB^-1 G_B.)."""
    import scipy.linalg as sl
    s = np.linalg.svd(XV, compute_uv=False)
    r = int((s > 1e-8 * s[0]).sum())
    _, _, piv = sl.qr(XV.T, pivoting=True, mode="economic")
    B = sorted(int(i) for i in piv[:r])
    nV = len(G)
    for i in range(nV):
        if not (G[i][i] == ONE):
            raise ValueError("diagonal entry not 1")
        for j in range(i + 1, nV):
            if not (G[i][j] == G[j][i]):
                raise ValueError("not symmetric")
    GBB = [[G[a][b] for b in B] for a in B]
    pd, _ = is_pos_def(GBB)
    if not pd:
        raise ValueError("G_BB not positive definite")
    I = [[ONE if i == j else ZERO for j in range(r)] for i in range(r)]
    GBBinv = ksolve(GBB, I)
    W = [[ksum(GBBinv[k][t] * G[B[t]][v] for t in range(r)) for v in range(nV)] for k in range(r)]
    for a in range(nV):
        ga = [G[a][b] for b in B]
        for v in range(a, nV):
            if not (ksum(ga[t] * W[t][v] for t in range(r)) == G[a][v]):
                raise ValueError(f"core Gram matrix has rank above {r}")
    return r, B, W, GBB, GBBinv


def contacts_core(G, d):
    nV = len(G)
    c = None
    for i in range(nV):
        for j in range(i + 1, nV):
            if c is None or G[i][j] > c:
                c = G[i][j]
    E = [(i, j) for i in range(nV) for j in range(i + 1, nV) if G[i][j] == c]
    deg = [0] * nV
    for a, b in E:
        deg[a] += 1
        deg[b] += 1
    return c, E, min(deg)


def karb(x, Ds):
    return flint.arb(x.a) + flint.arb(x.b) * Ds


def embed(G, B, r, d, disc):
    """Ball enclosures of the rows of P* = [G_.B L^-T, 0] (r coordinates, then d - r zeros)."""
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    Ds = flint.arb(disc).sqrt()
    GBB = [[karb(G[a][b], Ds) for b in B] for a in B]
    L = [[flint.arb(0)] * r for _ in range(r)]
    for j in range(r):
        s = GBB[j][j] - sum((L[j][k] * L[j][k] for k in range(j)), flint.arb(0))
        if not s > 0:
            raise ValueError("Cholesky pivot not provably positive")
        L[j][j] = s.sqrt()
        for i in range(j + 1, r):
            L[i][j] = (GBB[i][j] - sum((L[i][k] * L[j][k] for k in range(j)), flint.arb(0))) / L[j][j]
    rows = []
    for a in range(len(G)):
        g = [karb(G[a][b], Ds) for b in B]
        y = []
        for i in range(r):
            y.append((g[i] - sum((L[i][k] * y[k] for k in range(i)), flint.arb(0))) / L[i][i])
        rows.append(y + [flint.arb(0)] * (d - r))
    flint.ctx.prec = old
    return rows


def frame(rows, dec, V, d):
    """An exactly orthogonal rational Q with P* Q close to the stored core (untrusted choice, exact orthogonality)."""
    P = np.array([[float(x.mid()) for x in r] for r in rows])
    Xf = np.array([[float(t) for t in dec[a]] for a in V])
    Xf /= np.linalg.norm(Xf, axis=1)[:, None]
    U, _, Vt = np.linalg.svd(P.T @ Xf)
    Q = dist.rational_orthogonal(U @ Vt, d)
    if not (Q.transpose() * Q == flint.fmpq_mat(d, d, [int(i == j) for i in range(d) for j in range(d)])):
        raise ValueError("Q not orthogonal")
    return Q


def rotate(rows, Q, d):
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    Qa = [[flint.arb(Q[j, k]) for k in range(d)] for j in range(d)]
    out = []
    for r in rows:
        out.append([sum((r[j] * Qa[j][k] for j in range(d) if not r[j].is_zero()), flint.arb(0)) for k in range(d)])
    flint.ctx.prec = old
    return out


def stored_unit(dec_row):
    """The exactly normalised stored point as balls (numerators exact, one square root)."""
    q = [Fraction(t) for t in dec_row]
    s2 = sum(x * x for x in q)
    nrm = flint.arb(flint.fmpq(s2.numerator, s2.denominator)).sqrt()
    return [flint.arb(flint.fmpq(x.numerator, x.denominator)) / nrm for x in q], q, s2


def rattler_gaps(PQ, dec, V, n, c, disc):
    """(R): every core-rattler pair (balls) and rattler-rattler pair (exact) strictly below c. Returns (ok, worst upper
    bound of <.,.> - c, number of pairs not shown below c)."""
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    Ds = flint.arb(disc).sqrt()
    cA = karb(c, Ds)
    Vs = set(V)
    rat = [i for i in range(n) if i not in Vs]
    units = {r: stored_unit(dec[r]) for r in rat}
    bad, worst = 0, None
    for r in rat:
        u = units[r][0]
        for row in PQ:
            ip = sum((x * y for x, y in zip(row, u)), flint.arb(0))
            dlt = ip - cA
            up = float(dlt.upper())
            worst = up if worst is None else max(worst, up)
            if not (dlt < 0):
                bad += 1
    for i, r in enumerate(rat):
        for r2 in rat[i + 1:]:
            q1, s1 = units[r][1], units[r][2]
            q2, s2 = units[r2][1], units[r2][2]
            ip = sum(x * y for x, y in zip(q1, q2))
            # <x_r, x_r2> = ip / sqrt(s1 s2) < c ?  (c in K: decide with balls at high precision, then exactly if needed)
            val = flint.arb(flint.fmpq(ip.numerator, ip.denominator)) / (flint.arb(flint.fmpq(s1.numerator, s1.denominator)) * flint.arb(flint.fmpq(s2.numerator, s2.denominator))).sqrt()
            dlt = val - cA
            up = float(dlt.upper())
            worst = up if worst is None else max(worst, up)
            if not (dlt < 0):
                bad += 1
    flint.ctx.prec = old
    return bad == 0, worst, bad


def distance(PQ, dec, V, n):
    """(D): upper bounds for the least-squares and the largest coordinate distance between Q* and the stored code (the
    rattlers coincide)."""
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    tot = flint.arb(0)
    mx = 0.0
    for row, a in zip(PQ, V):
        u = stored_unit(dec[a])[0]
        for x, y in zip(row, u):
            df = x - y
            tot = tot + df * df
            mx = max(mx, math.nextafter(float(abs(df).upper()), math.inf))
    D2 = tot.sqrt()
    flint.ctx.prec = old
    return math.nextafter(float(D2.upper()), math.inf), mx


def build(snapdir, d, n, disc, tol=1e-10):
    """Everything above; returns a dict with the exact objects and the log entries (raises on a failed check)."""
    K.D = disc
    entry, dec, Xf = X2.read_code(snapdir, d, n)
    V = propose_core(Xf, d, tol)
    G, dev, nval = recognise_core(Xf, V, disc)
    if G is None:
        return {"status": f"core Gram matrix not recognised (value {nval})", "core_size": len(V)}
    r, B, W, GBB, GBBinv = gram_rank_basis(Xf[V], G, d)
    c, E, mindeg = contacts_core(G, d)
    if mindeg < d:
        return {"status": "a core point has fewer than d exact contacts", "core_size": len(V)}
    rows = embed(G, B, r, d, disc)
    Q = frame(rows, dec, V, d)
    PQ = rotate(rows, Q, d)
    ok, worst, bad = rattler_gaps(PQ, dec, V, n, c, disc)
    D2, Dc = distance(PQ, dec, V, n)
    return {"status": "ok" if ok else "a rattler is not shown below c", "entry": entry, "dec": dec, "Xf": Xf, "V": V,
            "G": G, "r": r, "B": B, "W": W, "GBB": GBB, "GBBinv": GBBinv, "c": c, "E": E, "Q": Q, "PQ": PQ,
            "log": {"core_size": len(V), "rattlers": n - len(V), "core_rank": r, "cos": repr(c), "contacts": len(E),
                    "gram_distinct_values": nval, "gram_recognition_deviation": dev, "rattler_pairs_not_below_c": bad,
                    "rattler_max_minus_c_upper": worst, "dist_ls_upper": D2, "dist_coord_upper": Dc}}
