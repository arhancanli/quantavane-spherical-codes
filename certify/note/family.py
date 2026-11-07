# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""Does the construction of X* generalise to PGammaL(2,q)? A search in floating point (UNTRUSTED: nothing here is a
certificate; it reports what a numerical search finds).

  python3 -I certify/note/family.py OUT.json q [q ...] [--starts K] [--seed S]

For each q = p^e the program builds, from the group alone, the objects of the recipe of the note's model:
  A   the halves: the two orbits on PG(1,q) of each cyclic subgroup of order (q+1)/2 of PGL(2,q) (q odd; for q = 9
      these are the 72 pentads);
  Bk  (q = r^2) the pairs {pi, pi'} of pencils of circles at a point (circles = images of PG(1,r), the inversive plane
      of order r) whose directions differ by an element of order k of F_q^*/F_r^* (for q = 9 and k = 4 these are the
      40 pairs of non-perpendicular pencils);
  P,Q two points each, swapped by the elements outside an index-2 subgroup M (for q = p^2 the subgroup M(q) with
      chi(det) (-1)^frob = 1, which is M10 for q = 9; otherwise the kernel of chi(det), that is PSL(2,q) x <frob>).
For a union Omega of these orbits it decomposes the permutation module R^Omega into irreducible subrepresentations
(eigenspaces of a random symmetric element of the centraliser algebra, identified by their characters on a sample of
group elements), and for a profile (a list of irreducible types with multiplicities, e.g. 1 + 2 eps + psi + psi') it
minimises the largest cosine over the Gamma-invariant configurations x_{g a} = rho(g) x_a in that representation
(x_a in the fixed space of the stabiliser of a): SLSQP from random starts, then a Gauss-Newton polish of the active
equalities. The result is analysed: largest cosine, contact orbits, rattlers, rank, the flexes modulo rotations
(dense SVD when N d is small; otherwise in each linear character), invariant strictly positive stresses (LP) and the
second-order energy H on the flexes. With q = 9 and the profile of X* it must find X* (largest cosine 0.11736).
"""
import argparse, itertools, json, math, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [p for p in os.environ.get("QV_LIBS", "").split(os.pathsep) if p]   # optional library folders
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.seterr(all="ignore")
from scipy.optimize import minimize, linprog


def svd(M, compute_uv=True, full_matrices=True):
    """numpy's SVD (gesdd), falling back to LAPACK gesvd when gesdd does not converge"""
    try:
        return np.linalg.svd(M, compute_uv=compute_uv, full_matrices=full_matrices)
    except np.linalg.LinAlgError:
        import scipy.linalg
        return scipy.linalg.svd(M, compute_uv=compute_uv, full_matrices=full_matrices, lapack_driver="gesvd")


# ===================================================================================================== F_q, q = p^e
class Field:
    def __init__(self, p, e):
        self.p, self.e, self.q = p, e, p ** e
        q = self.q
        def to_poly(x):
            return [(x // p ** i) % p for i in range(e)]
        def to_int(c):
            return sum(int(ci) % p * p ** i for i, ci in enumerate(c))
        def pmulmod(a, b, f):           # f monic of degree e, given by its e low coefficients
            prod = [0] * (2 * e - 1)
            for i, ai in enumerate(a):
                for j, bj in enumerate(b):
                    prod[i + j] = (prod[i + j] + ai * bj) % p
            for k in range(2 * e - 2, e - 1, -1):     # x^k = x^(k-e) * x^e, x^e = -sum f_i x^i
                c = prod[k]
                if c:
                    prod[k] = 0
                    for i in range(e):
                        prod[k - e + i] = (prod[k - e + i] - c * f[i]) % p
            return prod[:e]
        for f in itertools.product(range(p), repeat=e):
            if f[0] == 0:
                continue
            x = [0] * e
            if e == 1:
                x = [(-f[0]) % p]
            else:
                x[1] = 1
            cur, k = list(x), 1
            one = [1] + [0] * (e - 1)
            while cur != one and k <= q:
                cur, k = pmulmod(cur, x, f), k + 1
            if k == q - 1:
                break
        self.f = f
        self.mul = np.zeros((q, q), dtype=np.int64)
        self.add = np.zeros((q, q), dtype=np.int64)
        for a in range(q):
            for b in range(q):
                pa, pb = to_poly(a), to_poly(b)
                self.add[a, b] = to_int([u + v for u, v in zip(pa, pb)])
                self.mul[a, b] = to_int(pmulmod(pa, pb, f))
        self.w = to_int(x)               # a primitive element
        self.neg = np.array([next(b for b in range(q) if self.add[a, b] == 0) for a in range(q)])
        self.inv = np.array([0] + [next(b for b in range(q) if self.mul[a, b] == 1) for a in range(1, q)])
        self.frob = np.array([self.pw(a, p) for a in range(q)])
        sq = set(int(self.mul[a, a]) for a in range(1, q))
        self.chi = np.array([0] + [1 if a in sq else -1 for a in range(1, q)])

    def pw(self, a, k):
        r = 1
        for _ in range(k):
            r = int(self.mul[r, a])
        return r


# ===================================================================================================== the group
def moebius_perm(F, a, b, c, d, k):
    """x -> (a phi^k(x) + b)/(c phi^k(x) + d) on PG(1,q) = {0..q-1, q = inf}."""
    q, INF = F.q, F.q
    out = []
    for x in range(q + 1):
        y = x
        if y != INF:
            for _ in range(k):
                y = int(F.frob[y])
        if y == INF:
            out.append(INF if c == 0 else int(F.mul[a, F.inv[c]]))
            continue
        num = int(F.add[F.mul[a, y], b])
        den = int(F.add[F.mul[c, y], d])
        out.append(INF if den == 0 else int(F.mul[num, F.inv[den]]))
    return np.array(out, dtype=np.int64)


def group_generators(F, which="PGammaL"):
    """Generators of PGammaL(2,q), or of its subgroups PGL, PSL, PSigmaL = PSL<frob>, M = PSL<w x^p> (M10 for q = 9),
    with their linear-character data (chi(det), frobenius exponent)."""
    one, zero, w = 1, 0, F.w
    w2 = int(F.mul[w, w])
    mone = int(F.neg[1])
    if which in ("PGammaL", "PGL"):
        gens = [(moebius_perm(F, one, one, zero, one, 0), (1, 0)),               # x + 1
                (moebius_perm(F, w, zero, zero, one, 0), (int(F.chi[w]), 0)),    # w x
                (moebius_perm(F, zero, one, one, zero, 0), (int(F.chi[mone]), 0))]   # 1/x, det -1
    else:
        gens = [(moebius_perm(F, one, one, zero, one, 0), (1, 0)),               # x + 1
                (moebius_perm(F, w2, zero, zero, one, 0), (1, 0)),               # w^2 x
                (moebius_perm(F, zero, mone, one, zero, 0), (1, 0))]             # -1/x
    if which in ("PGammaL", "PSigmaL") and F.e > 1:
        gens.append((moebius_perm(F, one, zero, zero, one, 1), (1, 1)))     # x^p
    if which == "M" and F.e > 1:
        gens.append((moebius_perm(F, w, zero, zero, one, 1), (int(F.chi[w]), 1)))   # w x^p
    return gens


def compose(p1, p2):
    """(p1 o p2)(x) = p1[p2[x]]"""
    return p1[p2]


def closure_perms(gens, limit=200000):
    e = tuple(range(len(gens[0])))
    seen, todo = {e}, [e]
    while todo:
        x = todo.pop()
        for g in gens:
            y = tuple(int(g[i]) for i in x)
            if y not in seen:
                seen.add(y)
                todo.append(y)
                if len(seen) > limit:
                    raise RuntimeError("group too large")
    return [np.array(x) for x in seen]


def eps_value(chdata, F, kind):
    """linear characters: 'det' = chi(det); 'frob' = (-1)^k (e even); 'M' = chi(det)(-1)^k."""
    cdet, k = chdata
    if kind == "det":
        return cdet
    if kind == "frob":
        return (-1) ** (k % 2) if F.e % 2 == 0 else 1
    if kind == "M":
        return cdet * ((-1) ** (k % 2) if F.e % 2 == 0 else 1)
    raise ValueError(kind)


# ===================================================================================================== the objects
def orbit_sets(start, gens, act):
    seen, todo = {start}, [start]
    while todo:
        x = todo.pop()
        for g in gens:
            y = act(g, x)
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return seen


def build_objects(F, gens):
    q, INF = F.q, F.q
    perms = [g for g, _ in gens]
    objs = {}
    act_set = lambda g, S: frozenset(int(g[i]) for i in S)
    # ---- halves (q odd): an element of PGL of order (q+1)/2 acting with two orbits of that size
    if F.p != 2:
        h = (q + 1) // 2
        pgl = closure_perms([g for g, cd in gens if cd[1] == 0])
        halves = None
        for g in pgl:
            k, x, ident = 1, g.copy(), np.arange(q + 1)
            while not (x == ident).all():
                x, k = g[x], k + 1
            if k != h:
                continue
            orbs, seen = [], set()
            for i in range(q + 1):
                if i in seen:
                    continue
                o, j = [], i
                while j not in o:
                    o.append(j)
                    j = int(g[j])
                seen.update(o)
                orbs.append(frozenset(o))
            if len(orbs) == 2 and all(len(o) == h for o in orbs):
                halves = orbs
                break
        allh = set()
        for o in halves:
            allh |= orbit_sets(o, perms, act_set)
        objs["A"] = sorted(allh, key=sorted)
    # ---- circles and pencils (q = r^2)
    if F.e % 2 == 0:
        r = F.p ** (F.e // 2)
        sub = [x for x in range(q) if F.pw(x, r) == x]
        assert len(sub) == r
        base = frozenset(sub + [INF])
        circles = sorted(orbit_sets(base, perms, act_set), key=sorted)
        cidx = {c: i for i, c in enumerate(circles)}
        pencils_at = {}
        for z in range(q + 1):
            through = [c for c in circles if z in c]
            groups = []
            for c in through:
                for grp in groups:
                    if all(len(c & c2) == 1 for c2 in grp):
                        grp.append(c)
                        break
                else:
                    groups.append([c])
            pencils_at[z] = [frozenset(cidx[c] for c in gp) for gp in groups]
        # directions of the pencils at inf, and the order of their ratio in F_q^* / F_r^*
        def direction(pen):
            c = circles[next(iter(pen))]
            xs = sorted(x for x in c if x != INF)
            return int(F.add[xs[1], F.neg[xs[0]]])
        def ratio_order(d1, d2):
            u = int(F.mul[d2, F.inv[d1]])
            k, v = 1, u
            while F.pw(v, r) != v:
                v, k = int(F.mul[v, u]), k + 1
            return k
        pen_inf = pencils_at[INF]
        dirs = [direction(pn) for pn in pen_inf]
        def act_pencil_pair(g, pp):
            return frozenset(frozenset(cidx[act_set(g, circles[c])] for c in pn) for pn in pp)
        by_k = {}
        for i, j in itertools.combinations(range(len(pen_inf)), 2):
            by_k.setdefault(ratio_order(dirs[i], dirs[j]), []).append(frozenset([pen_inf[i], pen_inf[j]]))
        for k, pairs in sorted(by_k.items()):
            allp = set()
            for pp in pairs:
                allp |= orbit_sets(pp, perms, act_pencil_pair)
            objs["B%d" % k] = sorted(allp, key=lambda pp: sorted(sorted(x) for x in pp))
        objs["_circles"] = circles
        objs["_cidx"] = cidx
    return objs


# ===================================================================================================== Omega and its orbits
class Action:
    """Gamma acting on Omega (a list of objects of several kinds), by generator permutations."""

    def __init__(self, F, gens, objs, kinds, eps_kind):
        self.F = F
        self.kinds = kinds
        pts = []
        for k in kinds:
            if k in ("P", "Q"):
                pts += [(k, 1), (k, -1)]
            else:
                pts += [(k, o) for o in objs[k]]
        self.pts = pts
        self.N = N = len(pts)
        index = {x: i for i, x in enumerate(pts)}
        circles, cidx = objs.get("_circles"), objs.get("_cidx")
        self.gchar = [cd for _, cd in gens]
        self.geps = [eps_value(cd, F, eps_kind) for cd in self.gchar]
        self.gperm = []
        for (g, cd), ge in zip(gens, self.geps):
            img = []
            for k, o in pts:
                if k in ("P", "Q"):
                    img.append(index[(k, o * ge)])
                elif k == "A":
                    img.append(index[(k, frozenset(int(g[i]) for i in o))])
                else:
                    img.append(index[(k, frozenset(frozenset(cidx[frozenset(int(g[i]) for i in circles[c])]
                                                                for c in pn) for pn in o))])
            self.gperm.append(np.array(img, dtype=np.int64))
        # orbits, coset representatives (as permutations), their characters (chi(det), frob exponent mod e)
        self.orbit_of = -np.ones(N, dtype=np.int64)
        self.orbits, self.rep_perm, self.rep_char = [], {}, {}
        ident = np.arange(N)
        for s in range(N):
            if self.orbit_of[s] >= 0:
                continue
            oi = len(self.orbits)
            orb = [s]
            self.orbit_of[s] = oi
            self.rep_perm[s], self.rep_char[s] = ident, (1, 0)
            todo = [s]
            while todo:
                b = todo.pop()
                for gi, gp in enumerate(self.gperm):
                    c = int(gp[b])
                    if self.orbit_of[c] < 0:
                        self.orbit_of[c] = oi
                        self.rep_perm[c] = compose(gp, self.rep_perm[b])
                        cd0, cd1 = self.rep_char[b], self.gchar[gi]
                        self.rep_char[c] = (cd0[0] * cd1[0], (cd0[1] + cd1[1]) % F.e)
                        orb.append(c)
                        todo.append(c)
            self.orbits.append(sorted(orb))
        self.reps = [o[0] for o in self.orbits]
        # Schreier generators of the stabilisers of the representatives
        self.stab = []
        for o, a in zip(self.orbits, self.reps):
            S, seen = [], set()
            for b in o:
                for gi, gp in enumerate(self.gperm):
                    c = int(gp[b])
                    inv = np.argsort(self.rep_perm[c])
                    sp = compose(inv, compose(gp, self.rep_perm[b]))
                    key = sp.tobytes()
                    if key in seen or (sp == ident).all():
                        continue
                    seen.add(key)
                    cb, cg, cc = self.rep_char[b], self.gchar[gi], self.rep_char[c]
                    S.append((sp, (cb[0] * cg[0] * cc[0], (cb[1] + cg[1] - cc[1]) % F.e)))
            assert all(int(sp[a]) == a for sp, _ in S)
            self.stab.append(S)
        # orbitals: L[b, c] = label of the orbit of (b, c)
        self.L = -np.ones((N, N), dtype=np.int64)
        self.orbital_rep = []                       # (orbit index, rep point a, c)
        off = 0
        self.suborbits = []
        for oi, (o, a) in enumerate(zip(self.orbits, self.reps)):
            lab = np.arange(N)
            perms = [sp for sp, _ in self.stab[oi]]
            invs = [np.argsort(sp) for sp in perms]
            while True:
                old = lab.copy()
                for sp, si in zip(perms, invs):
                    lab = np.minimum(lab, lab[sp])
                    lab = np.minimum(lab, lab[si])
                while True:
                    l2 = lab[lab]
                    if (l2 == lab).all():
                        break
                    lab = l2
                if (lab == old).all():
                    break
            roots = sorted(set(int(x) for x in lab))
            rid = {r_: i for i, r_ in enumerate(roots)}
            labs = np.array([rid[int(x)] for x in lab]) + off
            self.suborbits.append([np.where(labs == off + i)[0] for i in range(len(roots))])
            for i, r_ in enumerate(roots):
                self.orbital_rep.append((oi, a, r_))
            for b in o:
                self.L[b, self.rep_perm[b]] = labs
            off += len(roots)
        self.K = off
        assert (self.L >= 0).all()
        # unordered classes: {b, c} -> min(L[b,c], L[c,b])
        self.U = np.minimum(self.L, self.L.T)

    def eps_of_rep(self, b, F, kind):
        return eps_value(self.rep_char[b], F, kind)

    def sample_elements(self, rng, k=40):
        out = []
        for _ in range(k):
            p = np.arange(self.N)
            cd = (1, 0)
            for _ in range(int(rng.integers(1, 9))):
                gi = int(rng.integers(len(self.gperm)))
                p = compose(self.gperm[gi], p)
                cd = (cd[0] * self.gchar[gi][0], (cd[1] + self.gchar[gi][1]) % self.F.e)
            out.append((p, cd))
        return out


# ===================================================================================================== representations
def irreducibles(act, rng):
    """Eigenspaces of a random symmetric element of the centraliser algebra, grouped into isomorphism types by their
    characters on a sample of elements."""
    N = act.N
    wts = rng.standard_normal(act.K)
    Z = wts[act.L]
    Z = Z + Z.T
    vals, vecs = np.linalg.eigh(Z)
    groups, i = [], 0
    scale = max(1.0, float(np.abs(vals).max()))
    while i < N:
        j = i + 1
        while j < N and vals[j] - vals[j - 1] < 1e-7 * scale:
            j += 1
        groups.append(vecs[:, i:j])
        i = j
    sample = [(np.arange(N), (1, 0))] + [(g, (1, 0)) for g in act.gperm] + act.sample_elements(rng, 60)
    types, out = [], []
    for Q in groups:
        ch = np.array([float(np.sum(Q[p] * Q)) for p, _ in sample])
        for t, (dim, chv) in enumerate(types):
            if dim == Q.shape[1] and np.abs(chv - ch).max() < 1e-6:
                out.append((t, Q))
                break
        else:
            types.append((Q.shape[1], ch))
            out.append((len(types) - 1, Q))
    return types, out, sample


def describe_linear(types, sample, act, F):
    """name the one-dimensional types: 1, det (chi(det)), frob ((-1)^k), M (chi(det)(-1)^k)."""
    names = {}
    # the characters on the sample need the (chi(det), k) of the sample elements: recompute from words is not
    # stored for gperm (index 1..), so compare with the generator values only (they determine a linear character)
    for t, (dim, ch) in enumerate(types):
        if dim != 1:
            continue
        gv = np.round(ch[1:1 + len(act.gperm)]).astype(int)
        for kind in ("triv", "det", "frob", "M"):
            ref = [1] * len(act.gperm) if kind == "triv" else [eps_value(cd, F, kind) for cd in act.gchar]
            if list(gv) == list(ref):
                names[t] = kind
                break
        else:
            names[t] = "lin?"
    return names


# ===================================================================================================== the optimisation
class Model:
    def __init__(self, act, comps, profile):
        """comps: list of (type, Q) eigenspaces; profile: list of (type, r). One eigenspace per type is used."""
        self.act = act
        N = act.N
        chosen = {}
        for t, Q in comps:
            if t not in chosen:
                chosen[t] = Q
        self.blocks = []          # (type, copy j, Q (N x m), Phi per orbit (N x f))
        for t, r in profile:
            Q = chosen[t]
            P = Q @ Q.T
            Phis = []
            for oi, a in enumerate(act.reps):
                ind = np.zeros((N, len(act.suborbits[oi])))
                for s, members in enumerate(act.suborbits[oi]):
                    ind[members, s] = 1.0
                V = P @ ind
                if V.size == 0:
                    Phis.append(np.zeros((N, 0)))
                    continue
                U, sv, _ = svd(V, full_matrices=False)
                f = int((sv > 1e-8 * max(1.0, sv.max())).sum())
                Phis.append(U[:, :f])
            for j in range(r):
                self.blocks.append((t, j, Q, Phis))
        # parameters: for each orbit and block, a vector of length f
        self.slices, pos = {}, 0
        for oi in range(len(act.reps)):
            for bi, (t, j, Q, Phis) in enumerate(self.blocks):
                f = Phis[oi].shape[1]
                self.slices[(oi, bi)] = slice(pos, pos + f)
                pos += f
        self.P = pos
        self.dim = sum(Q.shape[1] for (t, j, Q, Phis) in self.blocks)
        # quadratic forms of the orbitals: <x_a, x_c> = w^T Qk w for the representative pair (a, c)
        self.Qk = np.zeros((act.K, pos, pos))
        for k, (oi, a, c) in enumerate(act.orbital_rep):
            oj = int(act.orbit_of[c])
            gc = act.rep_perm[c]
            ginv = np.argsort(gc)
            for bi, (t, j, Q, Phis) in enumerate(self.blocks):
                A_ = Phis[oi]
                B_ = Phis[oj][ginv, :]                 # (P_g Phi)[x] = Phi[g^-1 x]
                M = A_.T @ B_
                si, sj = self.slices[(oi, bi)], self.slices[(oj, bi)]
                self.Qk[k][si, sj] += M / 2
                self.Qk[k][sj, si] += M.T / 2
        self.diag_k = [k for k, (oi, a, c) in enumerate(act.orbital_rep) if a == c]
        self.off_k = [k for k, (oi, a, c) in enumerate(act.orbital_rep) if a != c]
        self.orbit_params = [np.concatenate([np.arange(self.slices[(oi, bi)].start, self.slices[(oi, bi)].stop)
                                             for bi in range(len(self.blocks))]).astype(int)
                             for oi in range(len(act.reps))]

    def values(self, w):
        return np.einsum("i,kij,j->k", w, self.Qk, w)

    def normalise(self, w):
        w = w.copy()
        for idx in self.orbit_params:
            nrm = np.linalg.norm(w[idx])
            if nrm > 0:
                w[idx] /= nrm
        return w

    def optimise(self, rng, starts=8, iters=400):
        best = None
        Qoff = self.Qk[self.off_k]
        Qd = self.Qk[self.diag_k]
        for s in range(starts):
            w0 = self.normalise(rng.standard_normal(self.P))
            if any(np.linalg.norm(w0[idx]) == 0 for idx in self.orbit_params):
                return None                       # some orbit has no admissible vector
            t0 = float(self.values(w0)[self.off_k].max())
            z0 = np.concatenate([w0, [t0]])
            cons = [{"type": "ineq", "fun": lambda z: z[-1] - np.einsum("i,kij,j->k", z[:-1], Qoff, z[:-1]),
                     "jac": lambda z: np.hstack([-2 * np.einsum("kij,j->ki", Qoff, z[:-1]), np.ones((len(Qoff), 1))])},
                    {"type": "eq", "fun": lambda z: np.einsum("i,kij,j->k", z[:-1], Qd, z[:-1]) - 1,
                     "jac": lambda z: np.hstack([2 * np.einsum("kij,j->ki", Qd, z[:-1]), np.zeros((len(Qd), 1))])}]
            res = minimize(lambda z: z[-1], z0, jac=lambda z: np.concatenate([np.zeros(len(z) - 1), [1.0]]),
                           constraints=cons, method="SLSQP", options={"maxiter": iters, "ftol": 1e-15})
            w = self.normalise(res.x[:-1])
            mu = float(self.values(w)[self.off_k].max())
            if best is None or mu < best[0] - 1e-12:
                best = (mu, w)
        return best

    def polish(self, w, tol=1e-6, steps=60):
        """Gauss-Newton on: the orbitals within tol of the maximum share one value s; norms 1."""
        v = self.values(w)
        mu = v[self.off_k].max()
        act_k = [k for k in self.off_k if v[k] > mu - tol]
        z = np.concatenate([w, [mu]])
        for _ in range(steps):
            ww, s = z[:-1], z[-1]
            r = np.concatenate([np.einsum("i,kij,j->k", ww, self.Qk[act_k], ww) - s,
                                np.einsum("i,kij,j->k", ww, self.Qk[self.diag_k], ww) - 1])
            J = np.vstack([np.hstack([2 * np.einsum("kij,j->ki", self.Qk[act_k], ww), -np.ones((len(act_k), 1))]),
                           np.hstack([2 * np.einsum("kij,j->ki", self.Qk[self.diag_k], ww),
                                      np.zeros((len(self.diag_k), 1))])])
            dz = np.linalg.lstsq(J, -r, rcond=None)[0]
            z = z + dz
            if np.abs(r).max() < 1e-15:
                break
        w = self.normalise(z[:-1])
        v = self.values(w)
        return w, float(np.abs(v[act_k] - v[act_k].max()).max()), act_k

    def coordinates(self, w):
        """X (N x dim): the configuration in coordinates of V = sum of r copies of the chosen eigenspaces."""
        act, N = self.act, self.act.N
        cols = []
        for bi, (t, j, Q, Phis) in enumerate(self.blocks):
            Xb = np.zeros((N, Q.shape[1]))
            for oi, o in enumerate(act.orbits):
                v = Phis[oi] @ w[self.slices[(oi, bi)]]           # in R^Omega, fixed by the stabiliser
                for b in o:
                    vb = v[np.argsort(act.rep_perm[b])]
                    Xb[b] = Q.T @ vb
            cols.append(Xb)
        return np.hstack(cols)


# ===================================================================================================== analysis
def rotation_basis(d):
    idx = [(i, j) for i in range(d) for j in range(i + 1, d)]
    return idx


def analyse(model, w, F, eps_kind, full_limit=7000):
    act = model.act
    N = act.N
    X = model.coordinates(w)
    d = X.shape[1]
    G = X @ X.T
    off = ~np.eye(N, dtype=bool)
    mu = float(G[off].max())
    # unordered orbitals and their values
    cls = act.U[np.triu_indices(N, 1)]
    gv = G[np.triu_indices(N, 1)]
    uvals = {}
    for c_, v_ in zip(cls, gv):
        uvals.setdefault(int(c_), []).append(float(v_))
    spread = max(max(v) - min(v) for v in uvals.values())
    srt = sorted(((float(np.mean(v)), c_, len(v)) for c_, v in uvals.items()), reverse=True)
    gaps = [srt[i][0] - srt[i + 1][0] for i in range(len(srt) - 1)]
    contact_cls = [c_ for v_, c_, n_ in srt if v_ > mu - 1e-9]
    gap = (mu - max([v_ for v_, c_, n_ in srt if v_ <= mu - 1e-9], default=-1.0))
    E = [(a, b) for a, b in zip(*np.triu_indices(N, 1)) if act.U[a, b] in contact_cls]
    E = np.array(E, dtype=np.int64).reshape(-1, 2)
    deg = np.bincount(E.ravel(), minlength=N)
    rank = int(np.linalg.matrix_rank(X, tol=1e-9))
    rep = {"N": N, "dim": d, "largest_cosine": mu, "orbital_spread": spread, "contact_orbits": len(contact_cls),
           "contacts": int(len(E)), "gap_below_contacts": float(gap), "rank": rank,
           "degrees_by_orbit": [int(deg[o[0]]) for o in act.orbits], "orbit_sizes": [len(o) for o in act.orbits],
           "rattler_orbits": [oi for oi, o in enumerate(act.orbits) if deg[o[0]] == 0]}
    # ---------------- invariant stresses: equilibrium at the orbit representatives (enough, by invariance)
    ccls = sorted(contact_cls)
    m = len(ccls)
    rows = []
    for oi, a in enumerate(act.reps):
        R = np.zeros((d, m))
        for j, c_ in enumerate(ccls):
            nb = np.where(act.U[a] == c_)[0]
            nb = nb[nb != a]
            if len(nb):
                R[:, j] = (X[nb] - mu * X[a]).sum(axis=0)
        rows.append(R)
    EQ = np.vstack(rows)
    sizes = np.array([sum(1 for (a, b) in E if act.U[a, b] == c_) for c_ in ccls], dtype=float)
    # the invariant stresses: the null space of EQ (numerically), omega = Ns z; LP: maximise t subject to
    # Ns z >= t, sizes . Ns z = 1
    _, sv, Vt = svd(EQ)
    rkq = int((sv > 1e-8 * sv.max()).sum())
    Ns = Vt[rkq:].T
    rep["invariant_stress_dimension"] = int(Ns.shape[1])
    omega = None
    rep["strictly_positive_stress"] = False
    if Ns.shape[1]:
        nz = Ns.shape[1]
        cvec = np.zeros(nz + 1); cvec[-1] = -1
        A_ub = np.hstack([-Ns, np.ones((m, 1))]); b_ub = np.zeros(m)
        A_eq = np.concatenate([sizes @ Ns, [0]])[None, :]
        lp = linprog(cvec, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=[(None, None)] * (nz + 1),
                     method="highs")
        if lp.status == 0:
            omega = Ns @ lp.x[:nz]
            rep["maxmin_invariant_stress_times_E"] = float(lp.x[-1] * len(E))
            rep["strictly_positive_stress"] = bool(lp.x[-1] > 1e-9)
            rep["equilibrium_residual"] = float(np.abs(EQ @ omega).max())
    om_e = None
    if omega is not None:
        om_e = np.array([omega[ccls.index(act.U[a, b])] for a, b in E])
    # ---------------- flexes, full (dense) when small
    if N * d <= full_limit:
        Rm = np.zeros((len(E) + N, N * d))
        for r_, (a, b) in enumerate(E):
            Rm[r_, a * d:(a + 1) * d] = X[b]
            Rm[r_, b * d:(b + 1) * d] = X[a]
        for a in range(N):
            Rm[len(E) + a, a * d:(a + 1) * d] = X[a]
        sv = svd(Rm, compute_uv=False)
        rkR = int((sv > 1e-8 * sv.max()).sum())
        nullity = N * d - rkR
        dimK = d * (d - 1) // 2 if rank == d else None
        rep["flex_nullity"] = nullity
        rep["flexes_mod_rotations"] = (nullity - dimK) if dimK is not None else None
        svs = np.sort(sv)
        rep["singular_values_around_the_gap"] = [float(svs[nullity - 1]) if 0 < nullity <= len(svs) else None,
                                                 float(svs[nullity]) if nullity < len(svs) else None]
        if rep["flexes_mod_rotations"] and rep["flexes_mod_rotations"] > 0 and omega is not None:
            _, s_, Vt = svd(Rm)
            Nul = Vt[rkR:].T                                  # (N d) x nullity
            # remove the rotations: K = {(A x_a)_a}
            Kb = []
            for (i, j) in rotation_basis(d):
                A = np.zeros((d, d)); A[i, j], A[j, i] = 1, -1
                Kb.append((X @ A.T).ravel())
            Kb = np.array(Kb).T
            Qk_, _ = np.linalg.qr(Kb)
            Fl = Nul - Qk_ @ (Qk_.T @ Nul)
            U_, s2, _ = svd(Fl, full_matrices=False)
            Fl = U_[:, :int((s2 > 1e-6).sum())]
            # H_omega on the flexes for the max-min invariant stress
            Hm = np.zeros((Fl.shape[1], Fl.shape[1]))
            for i in range(Fl.shape[1]):
                for j in range(i, Fl.shape[1]):
                    yi, yj = Fl[:, i].reshape(N, d), Fl[:, j].reshape(N, d)
                    h = np.sum(om_e * (np.sum(yi[E[:, 0]] * yj[E[:, 1]], 1) + np.sum(yj[E[:, 0]] * yi[E[:, 1]], 1)
                                       - mu * np.sum(yi[E[:, 0]] * yj[E[:, 0]], 1) - mu * np.sum(yi[E[:, 1]] * yj[E[:, 1]], 1)))
                    Hm[i, j] = Hm[j, i] = h
            rep["H_on_flexes_maxmin_stress_eigs"] = [float(x) for x in np.linalg.eigvalsh(Hm)]
    # ---------------- linear characters: the forced-flex count and the flexes of each linear type
    rep["linear_sectors"] = linear_sectors(model, X, E, mu, F, om_e)
    return rep, X


def linear_sectors(model, X, E, mu, F, om_e):
    act = model.act
    N, d = X.shape
    # R_g on V for the generators and for the Schreier generators: blocks Q^T P_g Q
    def Rmat(perm):
        # (P_g v)[g(a)] = v[a]  =>  coordinates: Q^T P_g Q, with (P_g Q)[g(a)] = Q[a]: P_g Q = Q[g^-1]
        out = np.zeros((d, d))
        pos = 0
        for (t, j, Q, Phis) in model.blocks:
            m = Q.shape[1]
            out[pos:pos + m, pos:pos + m] = Q.T @ Q[np.argsort(perm)]
            pos += m
        return out
    res = {}
    edges = {(int(a), int(b)) for a, b in E}
    for kind in ("triv", "det", "frob", "M"):
        if kind == "frob" and F.e % 2:
            continue
        if kind in ("det", "M") and F.p == 2:
            continue
        psi = (lambda cd: 1) if kind == "triv" else (lambda cd, k_=kind: eps_value(cd, F, k_))
        # tangent vectors at the representatives in the psi-eigenspace of the stabiliser
        Ys = []
        for oi, a in enumerate(act.reps):
            rows = [X[a][None, :]]
            st = act.stab[oi]
            pick = np.random.default_rng(oi).permutation(len(st))[:300]
            for sp, cd in (st[int(i)] for i in pick):
                rows.append(Rmat(sp) - psi(cd) * np.eye(d))
            M = np.vstack(rows)
            _, s, Vt = svd(M)
            rk = int((s > 1e-8).sum())
            Ys.append(Vt[rk:].T)
        tdim = sum(Y.shape[1] for Y in Ys)
        offs = np.cumsum([0] + [Y.shape[1] for Y in Ys])
        # the rigidity map on contact-orbit representatives
        ccls = sorted({int(act.U[a, b]) for a, b in E})
        rows = []
        for c_ in ccls:
            a, b = next((int(a), int(b)) for a, b in E if act.U[a, b] == c_)
            row = np.zeros(tdim)
            for (u, v) in ((a, b), (b, a)):
                ou = int(act.orbit_of[u])
                gu = act.rep_perm[u]
                Ru = Rmat(gu) * psi(act.rep_char[u])
                row[offs[ou]:offs[ou + 1]] += X[v] @ Ru @ Ys[ou]
            rows.append(row)
        Mr = np.array(rows) if rows else np.zeros((0, tdim))
        sv = svd(Mr, compute_uv=False) if Mr.size else np.zeros(0)
        rk = int((sv > 1e-8 * max(1.0, sv.max() if sv.size else 1.0)).sum())
        # psi-rotations: skew A with R_g A R_g^T = psi(g) A for the generators
        idx = rotation_basis(d)
        cols = []
        for (i, j) in idx:
            A = np.zeros((d, d)); A[i, j], A[j, i] = 1, -1
            col = []
            for gp, cd in zip(act.gperm, act.gchar):
                Rg = Rmat(gp)
                col.append((Rg @ A @ Rg.T - psi(cd) * A).ravel())
            cols.append(np.concatenate(col))
        Ms = np.array(cols).T
        s2 = svd(Ms, compute_uv=False)
        kdim = len(idx) - int((s2 > 1e-8).sum())
        # edges: psi-equivariant functions on the contacts (sign-consistent orbits)
        ecount = 0
        for c_ in ccls:
            a, b = next((int(a), int(b)) for a, b in E if act.U[a, b] == c_)
            sign = {(a, b): 1}
            todo, ok = [(a, b)], True
            while todo and ok:
                e_ = todo.pop()
                for gp, cd in zip(act.gperm, act.gchar):
                    f_ = (int(gp[e_[0]]), int(gp[e_[1]]))
                    f_ = f_ if f_ in edges else (f_[1], f_[0])
                    sg = sign[e_] * psi(cd)
                    if f_ in sign:
                        if sign[f_] != sg:
                            ok = False
                            break
                    else:
                        sign[f_] = sg
                        todo.append(f_)
            ecount += 1 if ok else 0
        res[kind] = {"m_T": tdim, "m_K": kdim, "m_E": ecount, "count_T_minus_K_minus_E": tdim - kdim - ecount,
                     "flex_dim_mod_rotations": tdim - rk - kdim}
        if tdim - rk - kdim > 0 and om_e is not None:
            # the flexes of this sector, on all points, modulo the rotations of this sector; H for the stress om_e
            _, s_, Vt = svd(Mr) if Mr.size else (None, np.zeros(0), np.eye(tdim))
            Zk = Vt[rk:].T
            full = []
            for z in Zk.T:
                Yall = np.zeros((N, d))
                for oi, o in enumerate(act.orbits):
                    y0 = Ys[oi] @ z[offs[oi]:offs[oi + 1]]
                    for b in o:
                        Yall[b] = psi(act.rep_char[b]) * (Rmat(act.rep_perm[b]) @ y0)
                full.append(Yall.ravel())
            full = np.array(full).T
            _, s3, V3 = svd(Ms)
            kb = V3[int((s3 > 1e-8).sum()):].T
            rots = []
            for col in kb.T:
                A = np.zeros((d, d))
                for c_, (i, j) in zip(col, idx):
                    A[i, j], A[j, i] = c_, -c_
                rots.append((X @ A.T).ravel())
            if rots:
                Qr, _ = np.linalg.qr(np.array(rots).T)
                full = full - Qr @ (Qr.T @ full)
            U_, s4, _ = svd(full, full_matrices=False)
            Fl = U_[:, :int((s4 > 1e-6 * max(1.0, s4.max())).sum())]
            Hm = np.zeros((Fl.shape[1], Fl.shape[1]))
            for i in range(Fl.shape[1]):
                for j in range(i, Fl.shape[1]):
                    yi, yj = Fl[:, i].reshape(N, d), Fl[:, j].reshape(N, d)
                    Hm[i, j] = Hm[j, i] = float(np.sum(om_e * (np.sum(yi[E[:, 0]] * yj[E[:, 1]], 1)
                                                               + np.sum(yj[E[:, 0]] * yi[E[:, 1]], 1)
                                                               - mu * np.sum(yi[E[:, 0]] * yj[E[:, 0]], 1)
                                                               - mu * np.sum(yi[E[:, 1]] * yj[E[:, 1]], 1))))
            res[kind]["H_eigs_maxmin_stress"] = [float(x) for x in np.linalg.eigvalsh(Hm)]
    return res


# ===================================================================================================== driver
def setup(q, which, kinds, eps_kind, rng):
    p = next(pp for pp in range(2, q + 1) if q % pp == 0)
    e = round(math.log(q, p))
    F = Field(p, e)
    gens = group_generators(F, which)
    objs = build_objects(F, group_generators(F, "PGammaL"))      # the objects of the whole group; act by a subgroup
    if any(k not in ("P", "Q") and k not in objs for k in kinds):
        return None
    act = Action(F, gens, objs, kinds, eps_kind)
    types, comps, sample = irreducibles(act, rng)
    lin = describe_linear(types, sample, act, F)
    mult = {}
    for t, Q in comps:
        mult[t] = mult.get(t, 0) + 1
    return F, act, types, comps, lin, mult


def profiles(types, lin, eps_kind, mode):
    """1 + 2 eps + (two, or up to three) non-linear types; mode 'pairs', ('dims', D) or ('upto3', dmax)."""
    t1 = [t for t, n in lin.items() if n == "triv"]
    te = [t for t, n in lin.items() if n == eps_kind]
    if not t1 or not te:
        return []
    nonlin = [t for t, (dim, _) in enumerate(types) if dim > 1]
    combos = []
    if mode == "pairs":
        combos = list(itertools.combinations(nonlin, 2))
    elif mode[0] == "dims":
        combos = [c for c in itertools.combinations(nonlin, 2) if sorted(types[t][0] for t in c) == sorted(mode[1])]
    elif mode[0] == "pairs_max":
        combos = [c for c in itertools.combinations(nonlin, 2) if all(types[t][0] <= mode[1] for t in c)]
    elif mode[0] == "upto3":
        for k in (2, 3):
            combos += [c for c in itertools.combinations(nonlin, k) if 3 + sum(types[t][0] for t in c) <= mode[1]]
    return [[(t1[0], 1), (te[0], 2)] + [(t, 1) for t in c] for c in combos]


CASES = [
    # q, group, kinds, eps (the character by which P and Q are swapped), profiles, how many to analyse
    (9, "PGammaL", ["A", "B4", "P", "Q"], "M", "pairs", 3),
    (9, "PGammaL", ["A", "B2", "P", "Q"], "M", "pairs", 3),
    (9, "PGammaL", ["A", "P", "Q"], "M", "pairs", 3),
    (9, "PGL", ["A", "B4", "P", "Q"], "det", ("upto3", 22), 3),
    (9, "PSigmaL", ["A", "B4", "P", "Q"], "frob", ("upto3", 22), 3),
    (25, "PGammaL", ["A", "B6", "P", "Q"], "det", ("dims", (25, 26)), 3),
    (25, "PGammaL", ["A", "B3", "P", "Q"], "det", ("dims", (25, 26)), 3),
    (25, "PGammaL", ["A", "P", "Q"], "det", ("dims", (25, 26)), 3),
    (25, "PGammaL", ["A", "B6", "P", "Q"], "M", ("dims", (25, 26)), 3),
    (25, "PGammaL", ["A", "B3", "P", "Q"], "M", ("dims", (25, 26)), 3),
    (7, "PGammaL", ["A", "P", "Q"], "det", "pairs", 3),
    (11, "PGammaL", ["A", "P", "Q"], "det", "pairs", 3),
    (13, "PGammaL", ["A", "P", "Q"], "det", "pairs", 3),
    (17, "PGammaL", ["A", "P", "Q"], "det", ("dims", (17, 18)), 3),
    (19, "PGammaL", ["A", "P", "Q"], "det", ("dims", (19, 20)), 3),
    (27, "PGammaL", ["A", "P", "Q"], "det", ("pairs_max", 60), 3),
    (16, "PGammaL", ["B5", "P", "Q"], "frob", "pairs", 3),
]


def xstar_restarts(out, starts, seed):
    """q = 9, Omega = A u B4 u P u Q, the profile 1 + 2 eps + 9 + 10 of X*: run SLSQP from many random starts; record
    every local minimum found and compare the best with the thirteen values of X* (data/values.gp, numerically)."""
    rng = np.random.default_rng(seed)
    F, act, types, comps, lin, mult = setup(9, "PGammaL", ["A", "B4", "P", "Q"], "M", rng)
    best = None
    for prof in profiles(types, lin, "M", "pairs"):
        model = Model(act, comps, prof)
        r = model.optimise(rng, starts=4)
        if r and (best is None or r[0] < best[0]):
            best = (r[0], prof, model)
    mu0, prof, model = best
    tv_c = (math.sqrt(145) - 5) / 60
    finals, smallest, keep = [], [], {}
    for s in range(starts):
        r = model.optimise(rng, starts=1)
        w, resid, _ = model.polish(r[1])
        vv = model.values(w)[model.off_k]
        if abs(vv.max() - tv_c) < 1e-9:
            keep.setdefault("X*" if abs(vv.min() + 0.709041) < 1e-5 else "X'" if abs(vv.min() + 0.676040) < 1e-5 else "?", w)
        finals.append(float(vv.max()))
        smallest.append(float(vv.min()))          # t0 = -0.70904 for X*, t0' = -0.67604 for X'
    kinds = ["X*" if abs(f - tv_c) < 1e-9 and abs(m + 0.709041) < 1e-5 else
             "X'" if abs(f - tv_c) < 1e-9 and abs(m + 0.676040) < 1e-5 else "other" for f, m in zip(finals, smallest)]
    # the values t_k of X*, numerically, against the polished optimum
    import re
    from fractions import Fraction
    vt = open(os.path.join(HERE, "data", "values.gp")).read()
    vec = lambda name: [Fraction(x.strip()) for x in re.search(r"^%s = \[(.*?)\];" % name, vt, re.M).group(1).split(",")]
    D1, D2 = vec("D1"), vec("D2")
    TC = [[Fraction(x.strip()) for x in row.split(",")] for row in re.search(r"^TC = \[(.*?)\];", vt, re.M).group(1).split(";")]
    r145 = math.sqrt(145)
    s1 = -math.sqrt(float(D1[0]) + float(D1[1]) * r145)
    s2 = -math.sqrt(float(D2[0]) + float(D2[1]) * r145)
    bas = [1, r145, s1, r145 * s1, s2, r145 * s2, s1 * s2, r145 * s1 * s2]
    tv = sorted(sum(float(c) * b for c, b in zip(row, bas)) for row in TC)
    r = model.optimise(rng, starts=4)
    w, resid, _ = model.polish(r[1])
    v = model.values(w)
    found = sorted(set(round(float(x), 9) for k, x in enumerate(v) if k in model.off_k))
    dev = max(min(abs(x - t) for t in tv) for x in v[model.off_k])
    rep = {"profile_dims": [types[t][0] for t, rr in prof for _ in range(rr)], "params": model.P, "dim": model.dim,
           "starts": starts, "final_largest_cosines": sorted(finals),
           "starts_reaching_c": int(sum(abs(x - tv[-1]) < 1e-9 for x in finals)),
           "starts_reaching_Xstar": kinds.count("X*"), "starts_reaching_Xprime": kinds.count("X'"),
           "other_local_minima": sorted(set(round(f, 6) for f, k in zip(finals, kinds) if k == "other")),
           "c": tv[-1], "values_found": found, "max_deviation_from_the_t_k": float(dev),
           "distinct_values_found": len(found)}
    # the analysis of one configuration of each kind (full SVD of the rigidity matrix, LP for the invariant stresses)
    rep["analysis"] = {}
    for kind, w in sorted(keep.items()):
        an, X = analyse(model, w, F, "M")
        vv = model.values(w)
        an["values"] = sorted(set(round(float(x), 9) for k, x in enumerate(vv) if k in model.off_k))
        rep["analysis"][kind] = an
    json.dump(rep, open(out, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in rep.items() if k != "analysis"}))


def main():
    if len(sys.argv) > 2 and sys.argv[2] == "--xstar":
        xstar_restarts(sys.argv[1], int(sys.argv[3]) if len(sys.argv) > 3 else 100, 20261107)
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--cases", default="all", help="comma-separated indices into CASES")
    ap.add_argument("--starts", type=int, default=6)
    ap.add_argument("--seed", type=int, default=20261107)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    sel = range(len(CASES)) if args.cases == "all" else [int(x) for x in args.cases.split(",")]
    results = []
    t0 = time.time()
    def dump():
        json.dump({"program": "family.py", "seconds": round(time.time() - t0, 1), "results": results},
                  open(args.out, "w"), indent=1, default=str)
    for ci in sel:
        q, which, kinds, eps_kind, mode, n_an = CASES[ci]
        st = setup(q, which, kinds, eps_kind, rng)
        case = {"case": ci, "q": q, "group": which, "kinds": kinds, "eps": eps_kind,
                "mode": mode if isinstance(mode, str) else list(mode)}
        if st is None:
            case["skipped"] = "an object kind is not defined for this q"
            results.append(case); dump(); continue
        F, act, types, comps, lin, mult = st
        case.update(N=act.N, orbit_sizes=[len(o) for o in act.orbits], orbitals=act.K,
                    types=[[dim, mult[t], lin.get(t, "")] for t, (dim, _) in enumerate(types)])
        profs = profiles(types, lin, eps_kind, mode)
        case["profiles"] = len(profs)
        opt = []
        for prof in profs:
            model = Model(act, comps, prof)
            best = model.optimise(rng, starts=args.starts)
            if best is None:
                continue
            mu, w = best
            w, resid, act_k = model.polish(w)
            mu = float(model.values(w)[model.off_k].max())
            opt.append((mu, prof, model, w, resid))
        opt.sort(key=lambda x: x[0])
        case["optimised"] = [{"profile_dims": [types[t][0] for t, r in prof for _ in range(r)], "dim": model.dim,
                              "largest_cosine": mu, "polish_residual": resid} for mu, prof, model, w, resid in opt]
        case["analysed"] = []
        print(json.dumps({k: case[k] for k in ("case", "q", "group", "kinds", "eps", "N", "orbit_sizes", "profiles")}),
              "best", [round(o[0], 6) for o in opt[:6]], round(time.time() - t0), flush=True)
        for mu, prof, model, w, resid in opt[:n_an]:
            try:
                rep, X = analyse(model, w, F, eps_kind)
            except Exception as ex:
                import traceback
                rep = {"analysis_error": repr(ex), "traceback": traceback.format_exc()[-1500:]}
            rep["profile_dims"] = [types[t][0] for t, r in prof for _ in range(r)]
            rep["polish_residual"] = resid
            case["analysed"].append(rep)
            short = {k: rep.get(k) for k in ("profile_dims", "dim", "largest_cosine", "rank", "contact_orbits",
                                             "rattler_orbits", "degrees_by_orbit", "strictly_positive_stress",
                                             "maxmin_invariant_stress_times_E", "flexes_mod_rotations",
                                             "gap_below_contacts", "analysis_error")}
            h = rep.get("H_on_flexes_maxmin_stress_eigs")
            short["H_eigs_min_max"] = [min(h), max(h)] if h else None
            print("   ", json.dumps(short), flush=True)
            print("    linear", json.dumps(rep.get("linear_sectors")), flush=True)
        results.append(case)
        dump()
    dump()


if __name__ == "__main__":
    main()
