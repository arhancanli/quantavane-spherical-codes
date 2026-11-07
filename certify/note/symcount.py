# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""Symmetry-forced flexes across the archived table (UNTRUSTED floating point: a screening, not a certificate).

  python3 -I certify/note/symcount.py SNAPDIR OUT.jsonl [--shard I/K] [--only d_n,d_n,...]

For a code X spanning R^d with contacts E, tangent space T, infinitesimal rotations K and a finite group H of orthogonal
maps permuting X, the rigidity map T -> R^E is H-equivariant with kernel the flexes Fl (containing K) and cokernel
dual to the stresses Str. So, as characters of H,
    chi_{Fl/K} - chi_{Str} = chi_T - chi_K - chi_E =: D,
    chi_T(g) = #{a : g a = a} (tr g - 1),  chi_K(g) = ((tr g)^2 - tr g^2)/2,  chi_E(g) = #{contacts fixed by g},
and for every irreducible character psi of H, m_psi(Fl/K) >= <D, psi>: if <D, psi> > 0, the symmetry forces a flex of
type psi (Proposition 4.4 of the note). For the full symmetry group the count forces a flex exactly when Fl/K is not
isomorphic to a subrepresentation of Str, and then it does so for no larger group and every subgroup forces at most as
much; so one group per entry suffices.

For every entry of the table this program reads the stored decimals and decides, in double precision:
  contacts: the pairs within TOL = 1e-9 of the largest inner product (and the gap to the next inner product);
  rattlers (points in no contact) and the rank of X (entries with rattlers or not spanning are only counted);
  m0 = N(d-1) - d(d-1)/2 - |E|, the count for the trivial group (Maxwell);
  the symmetry group H, only where the trivial group forces nothing and a flex exists or the rank was not computed:
    with few distinct inner products, individualisation and refinement on the class matrix (aut116.automorphisms,
    exhaustive, the order from a stabiliser chain), then the elements by closure; otherwise every orthogonal map
    permuting the points to 1e-7, found by backtracking on the images of a base of d independent points (pruned by
    point invariants and inner products, TLIM seconds); groups of more than CAP elements are not enumerated;
  the character table of H by Burnside's method (eigenvectors of the class multiplication matrices, from the conjugacy
    classes and the class multiplication coefficients counted over the whole group), and <D, psi> for every psi,
    rounded to integers (the largest rounding error is recorded); forced = sum over psi of max(0, <D, psi>) psi(1);
  when the size allows (min(N + |E|, N d) <= RLIM) and the result can matter (m0 <= 0, or the stress test below applies):
    f = dim Fl/K from the rank of the rigidity matrix (eigenvalues of R R^T or R^T R), and s = f - m0; and when f > 0, s > 0, N d <= HLIM and the table does not mark the code as needing
    more optimisation: a stress whose smallest entry is largest
    (linear programming), the matrix of H_omega on Fl/K (Section 2 of the note) with its eigenvalues, and for a strictly
    positive stress the modes of the stressed contact graph softer than the coordinate functions (Proposition 2.7) and
    the energy of each flex direction mode by mode; for f = 1 also the largest H_omega(F) over the stresses with every
    entry at least half that smallest entry.
Writes one JSON line per entry.
"""
import argparse, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [p for p in os.environ.get("QV_LIBS", "").split(os.pathsep) if p]   # optional library folders
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "census"))
sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.seterr(all="ignore")
from scipy.optimize import linprog
from scipy.spatial import cKDTree
import scipy.sparse as sp
from snapshot import load_index, read_points

TOL = 1e-9          # contacts: inner products within TOL of the largest
TOLG = 1e-7         # symmetry: images within TOLG of a point, inner products within 1e2*TOLG
CAP = 20000         # largest group enumerated
TLIM = 120.0        # seconds for the backtracking search of one entry
RLIM = 7000         # largest min(N + |E|, N d) for the rank of the rigidity matrix
HLIM = 2600         # largest N d for the flex basis and the stress energy
LPTIME = float(os.environ.get("QV_LP_TIME", "120"))   # seconds for each linear program
NMAX = 20000        # entries with more points are not read (the Gram matrix would need more than 3 GB)


def contacts(X):
    G = X @ X.T
    n = len(X)
    iu = np.triu_indices(n, 1)
    v = G[iu]
    M = float(v.max())
    sel = v >= M - TOL
    below = v[~sel]
    gap = float(M - below.max()) if len(below) else None
    E = np.stack([iu[0][sel], iu[1][sel]], 1)
    return G, M, E, gap


def point_labels(G, E, n):
    """Invariants of the points, equal on orbits: degree and power sums of the inner products (clustered)."""
    deg = np.bincount(E.ravel(), minlength=n).astype(float)
    off = G - np.eye(n)
    inv = np.stack([deg] + [(off ** k).sum(1) for k in (1, 2, 3, 4)], 1)
    order = np.lexsort(inv.T[::-1])
    lab = np.zeros(n, dtype=np.int64)
    cur = 0
    for a, b in zip(order[:-1], order[1:]):
        if np.max(np.abs(inv[a] - inv[b]) / (1 + np.abs(inv[a]))) > 1e-6:
            cur += 1
        lab[b] = cur
    return lab


def value_classes(G, gap=1e-8):
    """The inner products clustered into classes (consecutive sorted values closer than gap)."""
    n = len(G)
    iu = np.triu_indices(n, 1)
    v = G[iu]
    order = np.argsort(v)
    lab_sorted = np.concatenate([[0], np.cumsum(np.diff(v[order]) > gap)])
    lab = np.empty_like(lab_sorted)
    lab[order] = lab_sorted
    C = -np.ones((n, n), dtype=np.int64)
    C[iu] = lab
    C.T[iu] = lab
    return C, int(lab_sorted[-1]) + 1


def group(X, G, E, cap=CAP, tlim=TLIM):
    """All orthogonal maps permuting X (to TOLG), as (perms, matrices); None if more than cap or out of time.
    Few distinct inner products: individualisation and refinement on the class matrix (aut116.automorphisms, exhaustive,
    with the order from a stabiliser chain), then the elements by closure and their matrices by least squares, each
    checked. Otherwise: backtracking on the images of a base of d independent points."""
    n, d = X.shape
    C, k = value_classes(G)
    if k <= 64 and n <= 400:
        from aut116 import automorphisms
        order, _, _, gens, _ = automorphisms(C)
        if order > cap:
            return None, "order %d > %d" % (order, cap)
        ident = np.arange(n, dtype=np.int32)
        seen = {ident.tobytes(): ident}
        todo = [ident]
        gl = [np.asarray(g, dtype=np.int32) for g in gens]
        while todo:
            p = todo.pop()
            for g in gl:
                q = g[p]
                kq = q.tobytes()
                if kq not in seen:
                    seen[kq] = q
                    todo.append(q)
        if len(seen) != order:
            return None, "closure %d != order %d" % (len(seen), order)
        perms = np.array(list(seen.values()))
        P = np.linalg.solve(X.T @ X, X.T)                      # left inverse of X
        mats = np.array([(P @ X[p]).T for p in perms])        # R with R x_a = x_{p(a)}
        dev = max(float(np.abs(X @ R.T - X[p]).max()) for R, p in zip(mats, perms))
        if dev > 1e-6:
            return None, "a class automorphism is not orthogonal (deviation %.1e)" % dev
        return (perms, mats), None
    import scipy.linalg as sl
    _, _, piv = sl.qr(X.T, pivoting=True, mode="economic")
    base = [int(i) for i in piv[:d]]
    A = X[base]
    Ainv_T = np.linalg.inv(A).T                 # R = B^T A^{-T}
    lab = point_labels(G, E, n)
    tree = cKDTree(X)
    perms, mats = [], []
    t0 = time.time()
    tolip = 1e2 * TOLG
    stack = [[]]
    timed_out = False
    while stack:
        img = stack.pop()
        k = len(img)
        if k == d:
            B = X[img]
            R = B.T @ Ainv_T
            if np.max(np.abs(R.T @ R - np.eye(d))) > 1e-6:
                continue
            Y = X @ R.T
            dist, idx = tree.query(Y)
            if dist.max() > TOLG or len(set(idx.tolist())) != n:
                continue
            perms.append(idx.astype(np.int32))
            mats.append(R)
            if len(perms) > cap:
                return None, "more than %d elements" % cap
            continue
        if time.time() - t0 > tlim:
            timed_out = True
            break
        b = base[k]
        mask = lab == lab[b]
        if k:
            mask &= np.all(np.abs(G[img] - G[base[:k], b][:, None]) < tolip, axis=0)
            mask[img] = False
        for j in np.flatnonzero(mask)[::-1]:
            stack.append(img + [int(j)])
    if timed_out:
        return None, "search out of time (%d elements found)" % len(perms)
    return (np.array(perms), np.array(mats)), None


def characters(perms, mats, E, d):
    """Conjugacy classes, Burnside's character table, and the virtual character D on the classes."""
    order, n = perms.shape
    index = {p.tobytes(): i for i, p in enumerate(perms)}
    ident = np.arange(n, dtype=np.int32)
    if ident.tobytes() not in index:
        raise RuntimeError("identity missing")
    inv = np.empty_like(perms)
    rows = np.arange(order)[:, None]
    inv[rows, perms] = np.arange(n, dtype=np.int32)[None, :]
    rng = np.random.default_rng(1)

    def generated(gens):
        seen = {index[ident.tobytes()]}
        todo = list(seen)
        while todo:
            x = todo.pop()
            for s in gens:
                y = index[perms[s][perms[x]].tobytes()]
                if y not in seen:
                    seen.add(y)
                    todo.append(y)
        return seen
    gens = [int(rng.integers(order))]
    sub = generated(gens)
    while len(sub) < order:
        gens.append(int(rng.choice([i for i in range(order) if i not in sub])))
        sub = generated(gens)
    cls = -np.ones(order, dtype=np.int64)
    reps, sizes = [], []
    for x in range(order):
        if cls[x] >= 0:
            continue
        c = len(reps)
        cls[x] = c
        members, todo = [x], [x]
        while todo:
            y = todo.pop()
            for s in gens:
                z = index[perms[s][perms[y]][inv[s]].tobytes()]
                if cls[z] < 0:
                    cls[z] = c
                    members.append(z)
                    todo.append(z)
        reps.append(x)
        sizes.append(len(members))
    r = len(reps)
    sizes = np.array(sizes, dtype=float)
    # class multiplication coefficients: a[i, k, j] = #{(x, y) in C_i x C_k : x y = z_j}
    a = np.zeros((r, r, r))
    for j, z in enumerate(reps):
        Y = inv[:, perms[z]]                     # x^{-1} z for every x (rows)
        ky = np.array([cls[index[row.tobytes()]] for row in Y])
        np.add.at(a, (cls, ky, np.full(order, j)), 1.0)
    Bm = np.einsum("i,ikj->kj", rng.standard_normal(r), a)
    w, V = np.linalg.eig(Bm)
    e = int(np.flatnonzero(cls[reps] == cls[index[ident.tobytes()]])[0])
    omega = V / V[e][None, :]                    # columns: central characters, 1 at the identity
    deg2 = order / np.sum(np.abs(omega) ** 2 / sizes[:, None], 0)
    degs = np.sqrt(deg2)
    chi = (omega * degs[None, :]) / sizes[:, None]   # chi[j, psi] = psi(z_j)
    # D on the classes, averaged over each class
    fix = (perms == ident[None, :]).sum(1)
    tr = np.trace(mats, axis1=1, axis2=2)
    tr2 = np.einsum("gij,gji->g", mats, mats)
    pa, pb = perms[:, E[:, 0]], perms[:, E[:, 1]]
    fixE = (((pa == E[:, 0]) & (pb == E[:, 1])) | ((pa == E[:, 1]) & (pb == E[:, 0]))).sum(1)
    Dg = fix * (tr - 1) - (tr ** 2 - tr2) / 2 - fixE
    Dc = np.bincount(cls, weights=Dg, minlength=r) / sizes
    mult = (sizes[:, None] * Dc[:, None] * np.conj(chi)).sum(0) / order
    return {"order": int(order), "classes": r, "generators": len(gens), "degrees": degs, "mult": mult,
            "deg_check": float(np.sum(deg2)), "D1": float(Dc[e])}


def rank_rigidity(X, E):
    n, d = X.shape
    m = len(E)
    rows = np.concatenate([np.repeat(np.arange(m), d), np.repeat(np.arange(m), d), m + np.repeat(np.arange(n), d)])
    cols = np.concatenate([(E[:, 0][:, None] * d + np.arange(d)).ravel(), (E[:, 1][:, None] * d + np.arange(d)).ravel(),
                           (np.arange(n)[:, None] * d + np.arange(d)).ravel()])
    vals = np.concatenate([X[E[:, 1]].ravel(), X[E[:, 0]].ravel(), X.ravel()])
    R = sp.csr_matrix((vals, (rows, cols)), shape=(m + n, n * d))
    Gm = (R @ R.T) if m + n <= n * d else (R.T @ R)
    ev = np.linalg.eigvalsh(Gm.toarray())
    thr = 1e-12 * ev.max()
    rk = int((ev > thr).sum())
    pos = ev[ev > thr]
    zero = ev[ev <= thr]
    return R, rk, (float(zero.max()) if len(zero) else None, float(pos.min()) if len(pos) else None)


def energy(X, E, R, f, orb=None):
    """Flex basis of Fl/K, a max-min stress (a sparse linear program in the stresses; over the invariant stresses, one
    value per orbit of contacts, when the symmetry group is known: averaging over the group keeps a stress and does not
    lower its smallest entry, so the optimum is the same), and H_omega on Fl/K."""
    import scipy.linalg as sl
    n, d = X.shape
    a_, b_ = E[:, 0], E[:, 1]
    M = float(np.einsum("ij,ij->i", X[a_], X[b_]).max())
    k0 = f + d * (d - 1) // 2
    _, nul = sl.eigh((R.T @ R).toarray(), subset_by_index=[0, k0 - 1])
    rots = []
    for i in range(d):
        for j in range(i + 1, d):
            A = np.zeros((d, d)); A[i, j], A[j, i] = 1, -1
            rots.append((X @ A.T).ravel())
    Q, _ = np.linalg.qr(np.array(rots).T)
    Fl = nul - Q @ (Q.T @ nul)
    U, s3, _ = np.linalg.svd(Fl, full_matrices=False)
    Fl = U[:, :f]
    m = len(E)
    # equilibrium sum_b omega_ab (x_b - M x_a) = 0 at every point, as a sparse matrix (n d) x m
    rows = np.concatenate([(a_[:, None] * d + np.arange(d)).ravel(), (b_[:, None] * d + np.arange(d)).ravel()])
    cols = np.concatenate([np.repeat(np.arange(m), d), np.repeat(np.arange(m), d)])
    vals = np.concatenate([(X[b_] - M * X[a_]).ravel(), (X[a_] - M * X[b_]).ravel()])
    Eq = sp.csr_matrix((vals, (rows, cols)), shape=(n * d, m))
    A_eq = sp.vstack([sp.hstack([Eq, sp.csr_matrix((n * d, 1))]),
                      sp.csr_matrix(np.concatenate([np.ones(m), [0.0]])[None, :])]).tocsr()
    b_eq = np.concatenate([np.zeros(n * d), [1.0]])
    out = {}
    if orb is not None:
        r_ = int(orb.max()) + 1
        P = sp.csr_matrix((np.ones(m), (np.arange(m), orb)), shape=(m, r_))
        EqP = (Eq @ P).toarray()
        sizes_ = np.bincount(orb, minlength=r_).astype(float)
        A_eqr = np.vstack([np.hstack([EqP, np.zeros((n * d, 1))]), np.concatenate([sizes_, [0.0]])[None, :]])
        A_ubr = np.hstack([-np.eye(r_), np.ones((r_, 1))])
        cr = np.zeros(r_ + 1); cr[-1] = -1
        lp = linprog(cr, A_ub=A_ubr, b_ub=np.zeros(r_), A_eq=A_eqr, b_eq=b_eq, bounds=[(None, None)] * (r_ + 1),
                     method="highs", options={"time_limit": LPTIME})
        out["contact_orbits"] = r_
        if lp.status != 0:
            out["lp_status"] = int(lp.status)
            return out
        om = P @ lp.x[:r_]
    else:
        A_ub = sp.hstack([-sp.identity(m), sp.csr_matrix(np.ones((m, 1)))]).tocsr()
        c = np.zeros(m + 1); c[-1] = -1
        lp = linprog(c, A_ub=A_ub, b_ub=np.zeros(m), A_eq=A_eq, b_eq=b_eq, bounds=[(None, None)] * (m + 1),
                     method="highs", options={"time_limit": LPTIME})
        if lp.status != 0:
            out["lp_status"] = int(lp.status)
            return out
        om = lp.x[:m]
    out["maxmin_times_E"] = float(om.min() * m)
    out["equilibrium_residual"] = float(np.abs(Eq @ om).max())
    Y = [Fl[:, i].reshape(n, d) for i in range(f)]

    def hvec(yi, yj):          # the coefficient of each contact in H(yi, yj)
        return (np.sum(yi[a_] * yj[b_], 1) + np.sum(yj[a_] * yi[b_], 1)
                - M * np.sum(yi[a_] * yj[a_], 1) - M * np.sum(yi[b_] * yj[b_], 1))
    H = np.zeros((f, f))
    for i in range(f):
        for j in range(i, f):
            H[i, j] = H[j, i] = float(om @ hvec(Y[i], Y[j]))
    he = np.linalg.eigvalsh(H)
    out["H_maxmin_eigs"] = [float(he.min()), float(he.max())]
    if om.min() > 0:
        # Proposition 2.7: the modes of the stressed contact graph softer than the coordinate functions
        w = np.bincount(np.concatenate([a_, b_]), weights=np.concatenate([om, om]), minlength=n)
        L = sp.csr_matrix((np.concatenate([-om, -om]), (np.concatenate([a_, b_]), np.concatenate([b_, a_]))),
                          shape=(n, n)).toarray() + np.diag(w)
        mu, Phi = sl.eigh(L, np.diag(w))
        soft = mu < (1 - M) * (1 - 1e-7)
        out["soft_modes"] = int(soft.sum())
        out["mu2_over_1mc"] = float(mu[1] / (1 - M))
        parts = []
        for i in range(f):
            contrib = (1 - M - mu) * ((Phi.T @ (w[:, None] * Y[i])) ** 2).sum(1)
            parts.append([float(contrib[0]), float(contrib[1:][soft[1:]].sum()), float(contrib[~soft].sum())])
        out["H_parts_const_soft_rest"] = parts
        out["H_parts_check"] = float(max(abs(sum(pp) - H[i, i]) for i, pp in enumerate(parts)))   # must be ~0
        if f == 1:
            hv = hvec(Y[0], Y[0])
            lo = 0.5 * om.min()
            lp2 = linprog(-hv, A_ub=-sp.identity(m).tocsr(), b_ub=-lo * np.ones(m), A_eq=A_eq[:, :m], b_eq=b_eq,
                          bounds=[(None, None)] * m, method="highs", options={"time_limit": LPTIME})
            if lp2.status == 0:
                out["H_best_positive_stress"] = float(hv @ lp2.x)
    return out


def analyse(snapdir, e):
    t0 = time.time()
    d, n = e["d"], e["n"]
    if n > NMAX:                                    # the Gram matrix alone would not fit in memory
        return {"d": d, "n": n, "flagged": int(e["flagged"]), "optimal": int(e["optimal"]), "category": "too large"}
    X = np.array([[float(t) for t in row] for row in read_points(snapdir, e)])
    X /= np.linalg.norm(X, axis=1)[:, None]
    G, M, E, gap = contacts(X)
    rep = {"d": d, "n": n, "flagged": int(e["flagged"]), "optimal": int(e["optimal"]), "contacts": int(len(E)),
           "gap": gap}
    deg = np.bincount(E.ravel(), minlength=n)
    rep["rattlers"] = int((deg == 0).sum())
    if rep["rattlers"]:
        rep["category"] = "rattlers"
        return rep
    sv = np.linalg.svd(X, compute_uv=False)
    if sv[-1] < 1e-8 * sv[0]:
        rep["category"] = "not spanning"
        return rep
    rep["m0"] = m0 = int(n * (d - 1) - d * (d - 1) // 2 - len(E))
    R = None
    # the rank can change the outcome only when the trivial group forces nothing, or for the stress test below
    need_rank = m0 <= 0 or (not e["flagged"] and n * d <= HLIM)
    if need_rank and min(len(E) + n, n * d) <= RLIM:
        R, rk, gapR = rank_rigidity(X, E)
        rep["f"] = int(n * d - rk - d * (d - 1) // 2)
        rep["s"] = rep["f"] - m0
        rep["rank_gap"] = gapR
    f = rep.get("f")
    # the group is needed only where a flex exists (or may exist) and the trivial group forces none
    group_perms = None
    if m0 <= 0 and (f is None or f > 0):
        gr, why = group(X, G, E)
        if gr is None:
            rep["group"] = why
        else:
            perms, mats = gr
            if len(perms) == 1:
                rep.update(group_order=1, forced=0, forced_types=[])
            else:
                try:
                    ch = characters(perms, mats, E, d)
                except KeyError:                       # the maps found are not closed under products
                    ch = None
                    rep["group"] = "%d maps found, not closed under products" % len(perms)
            if len(perms) > 1 and ch is not None:
                group_perms = perms
                mult = ch["mult"]
                mr = np.round(mult.real).astype(int)
                degs = np.round(ch["degrees"]).astype(int)
                rep.update(group_order=ch["order"], classes=ch["classes"], generators=ch["generators"],
                           rounding=float(max(np.abs(mult - mr).max(), np.abs(ch["degrees"] - degs).max())),
                           deg_check=round(ch["deg_check"], 6), D1=round(ch["D1"], 6))
                rep["forced"] = int(sum(int(m) * int(g) for m, g in zip(mr, degs) if m > 0))
                rep["forced_types"] = sorted([[int(g), int(m)] for m, g in zip(mr, degs) if m > 0])
                rep["counts_by_degree"] = sorted([[int(g), int(m)] for m, g in zip(mr, degs)])
    # the stress test, only for codes the table does not mark as needing more optimisation (their contacts are not
    # resolved); without a stress there is nothing to test
    if f is not None and f > 0 and rep["s"] > 0 and n * d <= HLIM and not e["flagged"]:
        orb = None
        if group_perms is not None:                    # the orbits of the group on the contacts
            eidx = {(int(a), int(b)): k for k, (a, b) in enumerate(E)}
            orb = -np.ones(len(E), dtype=np.int64)
            nxt = 0
            for k in range(len(E)):
                if orb[k] >= 0:
                    continue
                a, b = int(E[k, 0]), int(E[k, 1])
                for p_ in group_perms:
                    x_, y_ = int(p_[a]), int(p_[b])
                    orb[eidx[(min(x_, y_), max(x_, y_))]] = nxt
                nxt += 1
        rep.update(energy(X, E, R, f, orb))
    forced = rep.get("forced")
    if f == 0:
        rep["category"] = "bar-rigid"
    elif m0 > 0:
        rep["category"] = "forced, trivial group"
    elif forced is not None and forced > 0:
        rep["category"] = "forced by symmetry"
    elif forced == 0:
        rep["category"] = "flexes not forced" if f is not None else "not forced (rank not computed)"
    else:
        rep["category"] = "undecided"
    rep["seconds"] = round(time.time() - t0, 1)
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("snapdir")
    ap.add_argument("out")
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    i, k = (int(x) for x in a.shard.split("/"))
    only = {tuple(int(t) for t in s.split("_")) for s in a.only.split(",") if s}
    done = set()
    if os.path.exists(a.out):
        for line in open(a.out):
            r = json.loads(line)
            done.add((r["d"], r["n"]))
    entries = [e for e in load_index() if (not only or (e["d"], e["n"]) in only)]
    entries = [e for j, e in enumerate(entries) if j % k == i and (e["d"], e["n"]) not in done]
    with open(a.out, "a") as fo:
        for e in entries:
            try:
                rep = analyse(a.snapdir, e)
            except Exception as ex:                                   # recorded, not hidden
                rep = {"d": e["d"], "n": e["n"], "category": "error", "error": repr(ex)[:300]}
            fo.write(json.dumps(rep) + "\n")
            fo.flush()


if __name__ == "__main__":
    main()
