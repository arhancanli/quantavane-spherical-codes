"""UNTRUSTED numerical analysis of the entries that the table notes as continuous families (paper, version 5,
Section 5.5). Nothing here is proved; the paper reports it as numerical evidence only.

  python3 -I certify/v5/families.py SNAPDIR OUT.json D,N [D,N ...]

For each entry: the contacts within 1e-10 of the largest cosine (double precision) and the core; the numerical rank of
the bar rigidity matrix (singular values, threshold 1e-9 relative) and so the number of nontrivial infinitesimal flexes;
a stress of largest smallest entry (linear program, total 1); and, when there is exactly one flex f (unit, orthogonal to
the rotations), the largest and the smallest value of H_omega(f) over the stresses omega >= 0 of total 1 (two linear
programs).
"""
import json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
import numpy as np
np.seterr(all="ignore")
from scipy.optimize import linprog
import exact2 as X2
import core


def analyse(snapdir, d, n):
    t0 = time.time()
    entry, dec, Xf = X2.read_code(snapdir, d, n)
    V = core.propose_core(Xf, d)
    XV = Xf[V]
    nV = len(V)
    G = XV @ XV.T
    iu = np.triu_indices(nV, 1)
    c = G[iu].max()
    E = [(int(a), int(b)) for a, b in zip(*iu) if G[a, b] >= c - 1e-10]
    m = len(E)
    g = d * (d - 1) // 2
    R = np.zeros((nV + m, nV * d))
    for a in range(nV):
        R[a, a * d:(a + 1) * d] = XV[a]
    for e, (a, b) in enumerate(E):
        R[nV + e, a * d:(a + 1) * d] = XV[b]
        R[nV + e, b * d:(b + 1) * d] = XV[a]
    _, s, Vt = np.linalg.svd(R, full_matrices=False)
    rank = int((s > 1e-9 * s[0]).sum())
    out = {"d": d, "n": n, "core": nV, "contacts": m, "rank": rank, "flexes": nV * d - rank - g,
           "sv_last": float(s[rank - 1]), "sv_next": float(s[rank]) if rank < len(s) else None}
    # stress LP: variables omega (m), mu (nV), t
    Aeq = np.zeros((nV * d + 1, m + nV + 1))
    for e, (a, b) in enumerate(E):
        Aeq[a * d:(a + 1) * d, e] += XV[b]
        Aeq[b * d:(b + 1) * d, e] += XV[a]
    for a in range(nV):
        Aeq[a * d:(a + 1) * d, m + a] = XV[a]
    Aeq[nV * d, :m] = 1.0
    beq = np.zeros(nV * d + 1)
    beq[-1] = 1.0
    Aub = np.zeros((m, m + nV + 1))
    Aub[np.arange(m), np.arange(m)] = -1
    Aub[:, -1] = 1
    res = linprog(np.r_[np.zeros(m + nV), -1.0], A_ub=Aub, b_ub=np.zeros(m), A_eq=Aeq, b_eq=beq,
                  bounds=[(None, None)] * (m + nV + 1), method="highs")
    out["lp_maxmin_stress"] = float(-res.fun) if res.status == 0 else None
    if out["flexes"] == 1:
        rot = []
        for p in range(d):
            for q in range(p + 1, d):
                A = np.zeros((d, d))
                A[p, q], A[q, p] = 1, -1
                rot.append((XV @ A.T).reshape(-1))
        Qr = np.linalg.qr(np.array(rot).T)[0]
        null = Vt[rank:].T
        P = null - Qr @ (Qr.T @ null)
        f = np.linalg.svd(P, full_matrices=False)[0][:, 0].reshape(nV, d)
        h = np.array([2 * f[a] @ f[b] - c * (f[a] @ f[a] + f[b] @ f[b]) for a, b in E])
        for sense, key in ((-1.0, "H_max_over_stresses"), (1.0, "H_min_over_stresses")):
            r2 = linprog(np.r_[sense * h, np.zeros(nV)], A_eq=Aeq[:, :m + nV], b_eq=beq,
                         bounds=[(0, None)] * m + [(None, None)] * nV, method="highs")
            out[key] = float(sense * r2.fun) if r2.status == 0 else None
    out["seconds"] = round(time.time() - t0, 1)
    return out


def main():
    snapdir, outp = sys.argv[1], sys.argv[2]
    res = [analyse(snapdir, *map(int, s.split(","))) for s in sys.argv[3:]]
    for r in res:
        print(json.dumps(r))
    json.dump({"program": "families.py", "entries": res}, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
