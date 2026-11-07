"""The two halves of the split contact graphs (8,72) and (9,330) are antipodal, exactly (MCOM version).

certify/v5/rigidity.py proves that the exact core of each entry has two contact components, each bar-rigid with a
strictly positive stress. Proposition 2.9 of the paper also needs each half to span R^d. Here the exact core Gram matrix
(recognised and verified by certify/v5/core.py: rank d, contacts, rattlers below c) is used to check, exactly, that
every point of the first component has an antipode (inner product -1) in the second; so the second half is minus the
first, and both span the span of the core, which is R^d (rank d, proved by core.py).

  nice -n 10 python3 -I check_split_halves.py SNAPDIR OUT.json
"""
import hashlib, json, os, sys

CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
sys.path.insert(0, os.path.join(CERT, "v5"))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
import core


def main():
    snap, out = sys.argv[1], sys.argv[2]
    res = {"program": "check_split_halves.py",
           "code_sha256": {"check_split_halves.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
                           "core.py": hashlib.sha256(open(os.path.join(CERT, "v5", "core.py"), "rb").read()).hexdigest()},
           "entries": []}
    for d, n in ((8, 72), (9, 330)):
        C = core.build(snap, d, n, 1)
        G, E, r = C["G"], C["E"], C["r"]
        m = len(G)
        adj = [[] for _ in range(m)]
        for a, b in E:
            adj[a].append(b)
            adj[b].append(a)
        comp = [-1] * m
        k = 0
        for s in range(m):
            if comp[s] >= 0:
                continue
            st = [s]
            comp[s] = k
            while st:
                u = st.pop()
                for w in adj[u]:
                    if comp[w] < 0:
                        comp[w] = k
                        st.append(w)
            k += 1
        h1 = [a for a in range(m) if comp[a] == 0]
        h2 = [a for a in range(m) if comp[a] == 1]
        minus_one = core.K(-1)
        anti = {}
        for a in h1:
            bs = [b for b in h2 if G[a][b] == minus_one]
            if len(bs) == 1:
                anti[a] = bs[0]
        ok = (k == 2 and len(h1) == len(h2) and len(anti) == len(h1) and len(set(anti.values())) == len(h2)
              and r == d and C["status"] == "ok")
        rec = {"d": d, "n": n, "core_status": C["status"], "core_rank": r, "components": k,
               "sizes": [len(h1), len(h2)], "antipodal_pairs": len(anti),
               "status": "HALVES ANTIPODAL, CORE RANK d" if ok else "FAIL"}
        print(rec, flush=True)
        res["entries"].append(rec)
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
