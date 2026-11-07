"""Contact structure of every code of the archived table (UNTRUSTED: it only proposes contacts; nothing is proved here).

  python3 -I certify/census/classify.py SNAPDIR OUT.jsonl [--shard I/K]

For each code: the largest cosine M of the code normalised in double precision, and for each tolerance tau in
TOLS the contacts (pairs within tau of M) and the non-rattlers (points with at least d contacts among the remaining
points, removed repeatedly), with the isostatic count |V|(d-1) - d(d-1)/2 + 1. Also the cosines of the pairs just
below the contacts (the gap that separates the contacts from the other pairs), and the number of antipodal pairs.
"""
import json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, HERE)
import numpy as np
import snapshot

TOLS = [1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14]


def structure(Xn, d, G, M, tol):
    A = G >= M - tol
    alive = np.ones(len(Xn), bool)
    while True:
        deg = (A & alive[None, :] & alive[:, None]).sum(1)
        drop = alive & (deg < d)
        if not drop.any():
            break
        alive &= ~drop
    nV = int(alive.sum())
    m = int(np.triu(A & alive[None, :] & alive[:, None], 1).sum())
    iso = nV * (d - 1) - d * (d - 1) // 2 + 1 if nV else 0
    degs = (A & alive[None, :] & alive[:, None]).sum(1)[alive]
    return {"tol": tol, "nV": nV, "m": m, "iso": iso,
            "deg_min": int(degs.min()) if nV else 0, "deg_max": int(degs.max()) if nV else 0}


def main():
    snapdir, out = sys.argv[1], sys.argv[2]
    shard = (0, 1)
    if "--shard" in sys.argv:
        shard = tuple(int(x) for x in sys.argv[sys.argv.index("--shard") + 1].split("/"))
    done = set()
    if os.path.exists(out):
        for line in open(out):
            r = json.loads(line)
            done.add((r["d"], r["n"]))
    idx = [e for e in snapshot.load_index() if e["d"] <= 32 and e["n"] <= 1024]
    with open(out, "a") as fo:
        for k, e in enumerate(idx):
            if k % shard[1] != shard[0] or (e["d"], e["n"]) in done:
                continue
            t0 = time.process_time()
            d, n = e["d"], e["n"]
            X = np.array([[float(x) for x in p] for p in snapshot.read_points(snapdir, e)])
            Xn = X / np.linalg.norm(X, axis=1, keepdims=True)
            G = Xn @ Xn.T
            np.fill_diagonal(G, -2.0)
            M = float(G.max())
            iu = np.triu_indices(n, 1)
            top = np.sort(G[iu])[::-1]
            gaps = M - top
            rec = {"d": d, "n": n, "cos": e["cos"], "flagged": e["flagged"], "optimal": e["optimal"],
                   "poly": e["poly"], "M": M, "norm_dev": float(np.abs(np.linalg.norm(X, axis=1) - 1).max()),
                   "pairs_within": {f"{t:.0e}": int((gaps <= t).sum()) for t in TOLS},
                   "next_gap_after": {}, "antipodal": int((G <= -1 + 1e-12).sum() // 2),
                   "structs": [structure(Xn, d, G, M, t) for t in TOLS]}
            # the first gap above 1e-6 after the pairs within each tolerance (how far the next pair is)
            for t in (1e-10, 1e-13):
                c = int((gaps <= t).sum())
                rec["next_gap_after"][f"{t:.0e}"] = float(gaps[c]) if c < len(gaps) else None
            rec["cpu"] = round(time.process_time() - t0, 3)
            fo.write(json.dumps(rec) + "\n")
            fo.flush()


if __name__ == "__main__":
    main()
