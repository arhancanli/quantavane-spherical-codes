"""Library of the version 5 exact programs (paper, version 5): reading a stored code; UNTRUSTED recognition of inner
products in K = Q or Q(sqrt D) (recognise_value, recognise_gram); and TRUSTED helpers in exact arithmetic over K: the rank
of a sparse matrix over a prime field F_p, possibly after multiplication by a random matrix (rank_mod; either way a lower
bound for the rank over K, since the reduction modulo p is a ring homomorphism on the elements with denominators prime
to p, and a product can only lower the rank), the reduction of vectors modulo p (kvec_rows_modp), and the matrix of the
second-order form H_omega on given flexes (H_matrix). Used by core.py, rigidity.py, subspace.py and jam116.py.
"""
import argparse, hashlib, json, os, random, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/flintlib"))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "census"))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import flint
import mpmath
from fractions import Fraction
from kfield import K, ZERO, ONE, ksum, prime_for, is_pos_def, solve as ksolve, fq

SHA = {}
for _f in ("exact2.py", "kfield.py"):
    SHA[_f] = hashlib.sha256(open(os.path.join(HERE, _f), "rb").read()).hexdigest()


# ===================================================================================== UNTRUSTED: reading, recognition
def read_code(snapdir, d, n):
    import snapshot
    idx = {(r["d"], r["n"]): r for r in snapshot.load_index() if r["listed"]}
    r = idx[(d, n)]
    dec = snapshot.read_points(snapdir, r)
    X = np.array([[float(t) for t in p] for p in dec])
    X /= np.linalg.norm(X, axis=1)[:, None]
    return r, dec, X


def recognise_value(v, disc, maxden=10 ** 4):
    q = Fraction(v).limit_denominator(maxden)
    if abs(float(q) - v) < 1e-11:
        return K(fq(q.numerator, q.denominator))
    if disc == 1:
        return None
    mpmath.mp.dps = 20
    rel = mpmath.pslq([mpmath.mpf(v), 1, mpmath.sqrt(disc)], maxcoeff=10 ** 5, maxsteps=10 ** 5,
                      tol=mpmath.mpf(10) ** -12)
    if rel and rel[0] != 0:
        a, b = fq(-rel[1], rel[0]), fq(-rel[2], rel[0])
        x = K(a, b)
        if abs(float(x) - v) < 1e-11:
            return x
    return None


def recognise_gram(X, disc):
    """UNTRUSTED. Every inner product recognised as an element of K with small height (values within 1e-9 of each
    other share one recognition). Returns the matrix (lists of K) or None, and the largest deviation."""
    n = len(X)
    G = X @ X.T
    keys = np.round(G, 9)
    cache = {}
    out = [[None] * n for _ in range(n)]
    dev = 0.0
    for i in range(n):
        out[i][i] = ONE
        for j in range(i + 1, n):
            k = float(keys[i, j])
            if k not in cache:
                cache[k] = recognise_value(float(G[i, j]), disc)
                if cache[k] is None:
                    return None, None, k
            x = cache[k]
            dev = max(dev, abs(float(x) - float(G[i, j])))
            out[i][j] = out[j][i] = x
    return out, dev, len(cache)


# ===================================================================================== TRUSTED: exact checks
def check_gram(G, B, d):
    """(G) exact. Returns W (d x n, lists of K), GBB (d x d) and GBBinv, or raises."""
    n = len(G)
    for i in range(n):
        if not (G[i][i] == ONE):
            raise ValueError("diagonal entry not 1")
        for j in range(i + 1, n):
            if not (G[i][j] == G[j][i]):
                raise ValueError("not symmetric")
    GBB = [[G[a][b] for b in B] for a in B]
    pd, piv = is_pos_def(GBB)
    if not pd:
        raise ValueError("G_BB not positive definite")
    I = [[ONE if i == j else ZERO for j in range(d)] for i in range(d)]
    GBBinv = ksolve(GBB, I)
    GB = [[G[b][v] for v in range(n)] for b in B]           # d x n
    W = [[ksum(GBBinv[k][t] * GB[t][v] for t in range(d)) for v in range(n)] for k in range(d)]
    # G = G_{.B} W exactly, entry by entry
    for a in range(n):
        ga = [G[a][b] for b in B]
        for v in range(a, n):
            if not (ksum(ga[t] * W[t][v] for t in range(d)) == G[a][v]):
                raise ValueError(f"G has rank above d (entry {a},{v})")
    return W, GBB, GBBinv


def contacts(G):
    n = len(G)
    c = None
    for i in range(n):
        for j in range(i + 1, n):
            if c is None or G[i][j] > c:
                c = G[i][j]
    E = [(i, j) for i in range(n) for j in range(i + 1, n) if G[i][j] == c]
    return c, E


def core(E, n, d):
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
    V = sorted(V)
    Vs = set(V)
    return V, [(a, b) for a, b in E if a in Vs and b in Vs]


def rank_mod(rows, ncols, p, sketch_seed=None):
    """Rank over F_p of the sparse matrix given as a list of dicts {col: value mod p}. With sketch_seed, the rank of
    S R for a random S with ncols rows (a lower bound for the rank of R, equal to it with high probability)."""
    if sketch_seed is None or len(rows) <= ncols:
        M = flint.nmod_mat(len(rows), ncols, [row.get(j, 0) for row in rows for j in range(ncols)], p)
        return M.rank()
    rng = np.random.default_rng(sketch_seed)
    acc = np.zeros((ncols, ncols), dtype=np.int64)
    # acc = S R, S random ncols x len(rows); done row by row of R to keep memory small
    for row in rows:
        s = rng.integers(0, p, size=ncols, dtype=np.int64)
        for j, v in row.items():
            acc[:, j] = (acc[:, j] + s * v) % p
    M = flint.nmod_mat(ncols, ncols, [int(x) for x in acc.reshape(-1)], p)
    return M.rank()


def kvec_rows_modp(vecs, p, r):
    return [{j: x.modp(p, r) for j, x in enumerate(v) if not x.iszero()} for v in vecs]


def apply_R(G, B, V, E, u, d):
    """R u exactly (u: list over the |V| d coordinates, K elements): tangency then contact components."""
    pos = {v: k for k, v in enumerate(V)}
    out = []
    for a in V:
        ia = pos[a]
        out.append(ksum(G[a][B[k]] * u[ia * d + k] for k in range(d) if not u[ia * d + k].iszero()))
    for a, b in E:
        ia, ib = pos[a], pos[b]
        out.append(ksum(G[a][B[k]] * u[ib * d + k] for k in range(d) if not u[ib * d + k].iszero()) +
                   ksum(G[b][B[k]] * u[ia * d + k] for k in range(d) if not u[ia * d + k].iszero()))
    return out


def rotations(V, W, GBBinv, d, comp=None):
    """The infinitesimal rotations u_a = G_BB^{-1} S W_{.a} (S = e_s e_t^T - e_t e_s^T), restricted to the points of V
    in comp (others zero) when comp is given."""
    out = []
    for s in range(d):
        for t in range(s + 1, d):
            u = []
            for a in V:
                if comp is not None and a not in comp:
                    u.extend([ZERO] * d)
                    continue
                Sw = [ZERO] * d
                Sw[s] = W[t][a]
                Sw[t] = -W[s][a]
                u.extend(ksum(GBBinv[k][l] * Sw[l] for l in range(d) if not Sw[l].iszero()) for k in range(d))
            out.append(u)
    return out


def check_stress(G, W, V, E, d, omega):
    """(S) exact: sum_b omega_ab W_{.b} is a multiple of W_{.a} for every a in V."""
    nb = {a: [] for a in V}
    for e, (a, b) in enumerate(E):
        nb[a].append((b, omega[e]))
        nb[b].append((a, omega[e]))
    for a in V:
        s = [ksum(w * W[k][b] for b, w in nb[a]) for k in range(d)]
        wa = [W[k][a] for k in range(d)]
        k0 = next(k for k in range(d) if not wa[k].iszero())
        lam = s[k0] / wa[k0]
        if any(not (s[k] == lam * wa[k]) for k in range(d)):
            return False
    return True


def H_matrix(GBB, V, E, c, omega, flexes, d):
    pos = {v: k for k, v in enumerate(V)}
    f = len(flexes)

    def ip(u, v):
        return ksum(u[i] * GBB[i][j] * v[j] for i in range(d) for j in range(d)
                    if not u[i].iszero() and not v[j].iszero())

    def blk(y, a):
        return y[pos[a] * d:(pos[a] + 1) * d]

    H = [[ZERO] * f for _ in range(f)]
    for i in range(f):
        for j in range(i, f):
            acc = ZERO
            for e, (a, b) in enumerate(E):
                ya, yb, za, zb = blk(flexes[i], a), blk(flexes[i], b), blk(flexes[j], a), blk(flexes[j], b)
                t = ip(ya, zb) + ip(yb, za) - c * ip(ya, za) - c * ip(yb, zb)
                if not t.iszero():
                    acc = acc + omega[e] * t
            H[i][j] = H[j][i] = acc
    return H


# ===================================================================================== UNTRUSTED: candidate stress
def lp_stress(X, V, E, d):
    """A floating-point stress maximising its smallest entry (sum 1)."""
    from scipy.optimize import linprog
    pos = {v: k for k, v in enumerate(V)}
    nV, m = len(V), len(E)
    XV = X[V]
    Aeq = np.zeros((nV * d + 1, m + nV + 1))
    for e, (a, b) in enumerate(E):
        ia, ib = pos[a], pos[b]
        Aeq[ia * d:(ia + 1) * d, e] += XV[ib]
        Aeq[ib * d:(ib + 1) * d, e] += XV[ia]
    for k in range(nV):
        Aeq[k * d:(k + 1) * d, m + k] = XV[k]
    Aeq[nV * d, :m] = 1.0
    beq = np.zeros(nV * d + 1)
    beq[-1] = 1.0
    Aub = np.zeros((m, m + nV + 1))
    Aub[np.arange(m), np.arange(m)] = -1
    Aub[:, -1] = 1
    cobj = np.zeros(m + nV + 1)
    cobj[-1] = -1
    res = linprog(cobj, A_ub=Aub, b_ub=np.zeros(m), A_eq=Aeq, b_eq=beq, bounds=[(None, None)] * (m + nV + 1),
                  method="highs")
    return (res.x[:m], -res.fun) if res.status == 0 else (None, None)


def components(V, E):
    par = {v: v for v in V}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b in E:
        ra, rb = find(a), find(b)
        if ra != rb:
            par[ra] = rb
    comp = {}
    for v in V:
        comp.setdefault(find(v), set()).add(v)
    return sorted(comp.values(), key=lambda s: min(s))
