"""The worked example of the paper (Section 5): every quantity of the exact certificate of one entry, with more detail
than exactcheck.py records. Same trusted code (exactcheck.py); this script only prints more of it.

  nice -n 10 python3 -I example.py D N OUT.json
"""
import json, os, sys
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "verify"))
import exactcheck as XC
import numpy as np
import snapshot, cq


def main():
    d, n, outp = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    e = idx[(d, n)]
    r4 = json.load(open(os.path.join(XC.CERT, "census", "results_v4", f"d{d}_n{n}.json")))
    st = r4["status"]
    tol = [a for a in r4["attempts"] if a.get("status") == st][0]["tol"]
    pts = snapshot.read_points(XC.SNAP, e)
    dec = [[t.strip() for t in p] for p in pts]
    X = np.array([[float(x) for x in p] for p in dec])
    X = X / np.linalg.norm(X, axis=1, keepdims=True)
    log = {}
    XC.check_entry(X, dec, d, n, tol, st, log)
    keep = {}
    cq.certify(X, d, n, tol, dec=dec, keep=keep)
    V, E, Esq, extra = keep["V"], keep["E"], keep["Esq"], keep["extra"]
    sq = XC.Square(len(V), d, Esq, keep["xs"], keep["Z"], keep["C"], keep["K"])
    inv = XC.Inverse(np.asarray(keep["R"], dtype=np.float64))
    bspec = (extra, keep["beta_f"]) if extra else None
    res = XC.square_test(sq, inv, "krawczyk", bspec=bspec, want_constants=True)
    eps = res["eps"]
    G = X @ X.T
    np.fill_diagonal(G, -2)
    deg = [0] * len(V)
    for a, b in E:
        deg[a] += 1
        deg[b] += 1
    out = {"d": d, "n": n, "tol": tol, "table_cos": e["cos"], "table_poly": e["poly"], "flagged": e["flagged"],
           "optimal": e["optimal"], "max_cos_double": float(G.max()), "nV": len(V), "contacts": len(E),
           "extra": len(extra), "iso_count": len(V) * (d - 1) - d * (d - 1) // 2 + 1, "degrees": sorted(deg),
           "unknowns": sq.U, "eps_max": XC.up_float(max(eps)), "eps_gamma": XC.up_float(eps[-1]),
           "stress_lower": [XC.dn_float(x) for x in res["stress_lo"]],
           "stress_upper": [XC.up_float(x) for x in res["stress_hi"]],
           "normC": XC.up_float(res["normC"]), "Jinv_entry_bound": XC.up_float(res["Jinv"]),
           "next_gap_double": float(G.max() - np.sort(G[np.triu_indices(n, 1)])[::-1][len(E)]) if len(E) < n * (n - 1) // 2 else None,
           "exactcheck": log}
    json.dump(out, open(outp, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k not in ("exactcheck", "stress_lower", "stress_upper")}, default=str))


if __name__ == "__main__":
    main()
