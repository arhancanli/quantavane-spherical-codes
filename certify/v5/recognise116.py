"""UNTRUSTED recognition for the entry (22,116), for the note "A jammed spherical code that is not infinitesimally
jammed" (certify/note/; written for version 5 of the census paper, which no longer contains the result): nothing here
is trusted; jam116.py checks everything it uses.

  python3 -I certify/v5/recognise116.py SNAPDIR OUT.json

1. The pairs of points of the stored code are grouped into classes by their inner product (double precision, rounded
   to 8 decimals); there are 13 classes.
2. The code is refined to 300 digits on the class system (unit norms, and <x_a, x_b> = t_k for every pair of class k,
   the t_k unknown) by chord steps: residuals in 320-digit arithmetic (mpmath), corrections from one pseudo-inverse
   in double precision.
3. Each t_k is recognised by PSLQ in the basis 1, sqrt145, s1, sqrt145 s1, s2, sqrt145 s2, s1 s2, sqrt145 s1 s2,
   where s1 = t_7 - t_9 and s2 = t_3 - t_10 (differences of the two real roots of two quartics, found to have their
   squares in Q(sqrt145) by PSLQ).
4. A numerical infinitesimal flex of the contact framework (the null space of the rigidity matrix modulo rotations,
   which is one-dimensional numerically) and its first-order change S of the Gram matrix; the sign of S on every pair.
5. Edge types (vertex classes by colour refinement of the coloured Gram matrix, and the number of common contact
   neighbours) and a linear-programming stress constant on edge types, of largest smallest entry.
"""
import json, os, sys, collections, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
import numpy as np
np.seterr(all="ignore")
import mpmath
from mpmath import mp, mpf
import exact2 as X2


def main():
    snapdir, out = sys.argv[1], sys.argv[2]
    d, n = 22, 116
    t0 = time.time()
    entry, dec, Xf = X2.read_code(snapdir, d, n)
    G = Xf @ Xf.T
    iu = np.triu_indices(n, 1)
    keys = np.round(G, 8)
    vals = sorted(set(keys[iu].tolist()))
    cls = {v: k for k, v in enumerate(vals)}
    pairs = list(zip(iu[0].tolist(), iu[1].tolist()))
    pc = [cls[keys[a, b]] for a, b in pairs]
    nc = len(vals)
    # ---- 2. refinement
    J = np.zeros((n + len(pairs), n * d + nc))
    for a in range(n):
        J[a, a * d:(a + 1) * d] = 2 * Xf[a]
    for r, (a, b) in enumerate(pairs):
        J[n + r, a * d:(a + 1) * d] = Xf[b]
        J[n + r, b * d:(b + 1) * d] = Xf[a]
        J[n + r, n * d + pc[r]] = -1.0
    Jp = np.linalg.pinv(J, rcond=1e-10)
    mp.dps = 320
    X = [[mpf(t) for t in row] for row in dec]
    X = [[x / mpmath.sqrt(mpmath.fsum(y * y for y in row)) for x in row] for row in X]
    t = [mpf(v) for v in vals]
    res = None
    for it in range(40):
        F = [mpmath.fsum(x * x for x in X[a]) - 1 for a in range(n)]
        F += [mpmath.fdot(X[a], X[b]) - t[pc[r]] for r, (a, b) in enumerate(pairs)]
        res = max(abs(f) for f in F)
        if res < mpf(10) ** -300:
            break
        scale = float(res)
        step = Jp @ (np.array([float(f / res) for f in F]))
        for a in range(n):
            for k in range(d):
                X[a][k] -= mpf(step[a * d + k]) * res
        for k in range(nc):
            t[k] -= mpf(step[n * d + k]) * res
    # ---- 3. recognition in the tower
    mp.dps = 300
    r145 = mpmath.sqrt(145)
    s1, s2 = t[7] - t[9], t[3] - t[10]
    rel1 = mpmath.pslq([s1 ** 2, 1, r145], maxcoeff=10 ** 30, maxsteps=10 ** 6)
    rel2 = mpmath.pslq([s2 ** 2, 1, r145], maxcoeff=10 ** 30, maxsteps=10 ** 6)
    basis = [mpf(1), r145, s1, r145 * s1, s2, r145 * s2, s1 * s2, r145 * s1 * s2]
    rels = [[int(v) for v in mpmath.pslq([x] + basis, maxcoeff=10 ** 30, maxsteps=10 ** 6)] for x in t]
    # ---- 4. numerical flex
    c = G[iu].max()
    E = [(int(a), int(b)) for a, b in zip(*iu) if G[a, b] >= c - 1e-10]
    R = np.zeros((n + len(E), n * d))
    for a in range(n):
        R[a, a * d:(a + 1) * d] = Xf[a]
    for e, (a, b) in enumerate(E):
        R[n + e, a * d:(a + 1) * d] = Xf[b]
        R[n + e, b * d:(b + 1) * d] = Xf[a]
    _, s, Vt = np.linalg.svd(R, full_matrices=False)
    rank = int((s > 1e-9 * s[0]).sum())
    null = Vt[rank:].T
    rot = []
    for p in range(d):
        for q in range(p + 1, d):
            A = np.zeros((d, d))
            A[p, q], A[q, p] = 1, -1
            rot.append((Xf @ A.T).reshape(-1))
    Qr = np.linalg.qr(np.array(rot).T)[0]
    P = null - Qr @ (Qr.T @ null)
    u, ss, _ = np.linalg.svd(P, full_matrices=False)
    f = u[:, 0].reshape(n, d)
    S = Xf @ f.T + f @ Xf.T
    sgn = {}
    flex_classes = set()
    for (a, b), k in zip(pairs, pc):
        v = S[a, b]
        if abs(v) > 1e-9:
            sgn[f"{a},{b}"] = 1 if v > 0 else -1
            flex_classes.add(k)
    # ---- 5. edge types and an LP stress constant on them
    col = {i: 0 for i in range(n)}
    for it in range(10):
        sig = {i: (col[i], tuple(sorted(collections.Counter((keys[i, j], col[j]) for j in range(n) if j != i).items())))
               for i in range(n)}
        ss_ = sorted(set(sig.values()))
        new = {i: ss_.index(sig[i]) for i in range(n)}
        if len(set(new.values())) == len(set(col.values())):
            col = new
            break
        col = new
    nb = collections.defaultdict(set)
    for a, b in E:
        nb[a].add(b)
        nb[b].add(a)
    types = sorted({(min(col[a], col[b]), max(col[a], col[b]), len(nb[a] & nb[b])) for a, b in E})
    etype = [types.index((min(col[a], col[b]), max(col[a], col[b]), len(nb[a] & nb[b]))) for a, b in E]
    from scipy.optimize import linprog
    A = np.zeros((n * d, len(types)))
    for (a, b), j in zip(E, etype):
        A[a * d:(a + 1) * d, j] += Xf[b] - c * Xf[a]
        A[b * d:(b + 1) * d, j] += Xf[a] - c * Xf[b]
    _, sv, vt = np.linalg.svd(A)
    nul = vt[int((sv > 1e-9 * sv[0]).sum()):].T
    cnt = np.bincount(etype, minlength=len(types)).astype(float)
    k = nul.shape[1]
    resl = linprog(np.r_[np.zeros(k), -1.0], A_ub=np.hstack([-nul, np.ones((len(types), 1))]),
                   b_ub=np.zeros(len(types)), A_eq=np.hstack([(cnt @ nul)[None, :], [[0.0]]]), b_eq=[1.0],
                   bounds=[(None, None)] * (k + 1), method="highs")
    w_lp = (nul @ resl.x[:k]).tolist()
    import scipy.linalg as sl
    _, _, piv = sl.qr(Xf.T, pivoting=True, mode="economic")
    B = sorted(int(i) for i in piv[:d])
    json.dump({"d": d, "n": n, "classes": vals, "pair_class": pc, "refine_residual": mpmath.nstr(res, 5),
               "t_300": [mpmath.nstr(x, 60) for x in t], "Delta1_rel": [int(v) for v in rel1],
               "Delta2_rel": [int(v) for v in rel2], "sigma1": 1 if s1 > 0 else -1, "sigma2": 1 if s2 > 0 else -1,
               "class_rel": rels, "flex_rank": rank, "flex_sign": sgn, "flex_classes": sorted(flex_classes),
               "contacts": E, "edge_types": [list(x) for x in types], "edge_type": etype, "lp_stress_by_type": w_lp,
               "lp_min": float(-resl.fun), "basis": B, "seconds": round(time.time() - t0, 1)}, open(out, "w"))
    print("ok", round(time.time() - t0, 1), "s; residual", mpmath.nstr(res, 3), "; rank", rank, "; LP min", -resl.fun)


if __name__ == "__main__":
    main()
