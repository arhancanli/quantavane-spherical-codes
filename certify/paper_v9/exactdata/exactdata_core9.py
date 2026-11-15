"""Census version 9, item A5: the exact descriptions of the strict codes certified through their exact cores
(above the size limit, run_above.py; the exact entries of version 8), in the format of strict9.py --data.

  nice -n 10 python3 -I exactdata_core9.py OUTDIR [--only d_n,...]

For every entry of class strict_above, strict_exact or strict_new in paper_v8/v8_verdicts.tsv, certify/v5/core.py
(released, unchanged; the field recorded by run_above.py) rebuilds the exact core and checks it exactly (unit diagonal,
rank d, the contacts at the exact cosine, every other core pair and every stored rattler below it). Written:
OUTDIR/d{d}_n{n}.json.gz with the field (D: K = Q(sqrt D), D = 1 for Q), the distinct inner products as pairs (a, b)
meaning a + b sqrt(D), the index matrix of the core's Gram matrix, the basis B, the contacts, the cosine, the core and
the rattlers. The certificate of strictness itself is the version 8 record (paper_v8/above/out*, paper_v8 verdicts).
"""
import gzip, hashlib, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, os.path.join(CERT, "v5"))
import warnings
warnings.filterwarnings("ignore")
import core
from kfield import K

SHA = {"exactdata_core9.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
       "core.py": hashlib.sha256(open(os.path.join(CERT, "v5", "core.py"), "rb").read()).hexdigest()}


def entries():
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(CERT, "paper_v8", "v8_verdicts.tsv")) if not l.startswith("#")]
    head = rows[0]
    out = []
    for r in rows[1:]:
        r = dict(zip(head, r))
        if r["certificate"] in ("strict_above", "strict_exact", "strict_new"):
            out.append((int(r["d"]), int(r["n"]), r["certificate"]))
    return out


def disc_of(d, n):
    for sub in ("out", "out_more"):
        p = os.path.join(CERT, "paper_v8", "above", sub, f"d{d}_n{n}.json")
        if os.path.exists(p):
            return int(json.load(open(p)).get("disc", 1))
    return None


def main():
    out = sys.argv[1]
    only = set(sys.argv[sys.argv.index("--only") + 1].split(",")) if "--only" in sys.argv else None
    os.makedirs(out, exist_ok=True)
    log = []
    for d, n, cert in entries():
        if only and f"{d}_{n}" not in only:
            continue
        js = os.path.join(out, f"d{d}_n{n}.json.gz")
        if os.path.exists(js):
            continue
        t0 = time.time()
        discs = [disc_of(d, n)] if disc_of(d, n) else [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15]
        C = None
        for D in discs:
            K.D = D
            C = core.build(SNAP, d, n, D)
            if C["status"] == "ok":
                break
        rec = {"d": d, "n": n, "certificate": cert, "status": C["status"]}
        if C["status"] == "ok":
            G, V, E, c, r, B = (C[k] for k in ("G", "V", "E", "c", "r", "B"))
            vals, idxm = {}, []
            for i in range(len(V)):
                row = []
                for j in range(len(V)):
                    x = G[i][j]
                    key = (str(x.a), str(x.b))
                    if key not in vals:
                        vals[key] = len(vals)
                    row.append(vals[key])
                idxm.append(row)
            desc = {"d": d, "n": n, "field_D": K.D, "cosine": [str(c.a), str(c.b)], "core_rank": r,
                    "nonrattlers": [int(v) for v in V], "B": [int(b) for b in B],
                    "elements": [list(k) for k in sorted(vals, key=vals.get)], "index_matrix": idxm,
                    "contacts": [[int(a), int(b)] for a, b in E], "program_sha256": SHA}
            with gzip.open(js, "wt") as fo:
                json.dump(desc, fo)
            rec.update(field_D=K.D, distinct=len(vals), contacts=len(E), core=len(V))
        rec["seconds"] = round(time.time() - t0, 1)
        log.append(rec)
        print(rec, flush=True)
    json.dump({"program_sha256": SHA, "entries": log}, open(os.path.join(out, "exactdata_core9_log.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
