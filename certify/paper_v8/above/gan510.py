"""Details of the exact core of (10,510) for Section 6 (census version 8): the distinct exact inner products and the
largest one below the contacts' cosine, from certify/v5/core.py (the same exact core as run_above.py).

  nice -n 10 python3 -I gan510.py OUT.json
"""
import hashlib, json, os, sys
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, os.path.join(CERT, "v5"))
import core

C = core.build(SNAP, 10, 510, 1)
G, V, E, c = C["G"], C["V"], C["E"], C["c"]
vals = {}
n = len(V)
deg = [0] * n
for a, b in E:
    deg[a] += 1
    deg[b] += 1
for a in range(n):
    for b in range(a + 1, n):
        key = repr(G[a][b])
        vals[key] = vals.get(key, 0) + 1
below = [G[a][b] for a in range(n) for b in range(a + 1, n) if not (G[a][b] == c)]
nxt = max(below, key=lambda x: float(x))
assert all(float(x) <= float(nxt) for x in below)
out = {"program_sha256": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(), "core_size": n,
       "contacts": len(E), "cos": repr(c), "next_value": repr(nxt), "gap": repr(c - nxt),
       "degree_min": min(deg), "degree_max": max(deg), "distinct_values": vals}
json.dump(out, open(sys.argv[1], "w"), indent=1)
print({k: v for k, v in out.items() if k != "distinct_values"})
