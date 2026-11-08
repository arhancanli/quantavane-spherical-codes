"""The census rule on codes fetched from the live table, and the exact decision of every positive certificate (census
version 7, the currency check of Section 6).

  nice -n 10 python3 -I run_live.py FETCHDIR SELECT.tsv OUTDIR [--shard I/K]

For each entry of SELECT.tsv with a fetched file FETCHDIR/dD_nN.txt: the contact structure at every tolerance
(census/classify.py), the attempts of the census rule (census/run_census.py: isostatic, then hyperstatic, then the
negative certificate, at most 3,100 unknowns), each by the census certifier (cq.py, unchanged), and for a PASS, LOCALMIN
or NONOPT the exact decision of exactcheck.py. Also, untrusted: the largest cosine of the fetched code and of the archived
one in double precision, and the lower bound (1/2) max_k |a_k - b_k| for the distance of a point (modulo rotations and
relabelling) between the two codes, where a and b are their sorted lists of inner products (an inner product moves by at
most 2 D when every point moves by at most D).
"""
import hashlib, json, os, sys, time, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "verify"))
import exactcheck as XC
import numpy as np
import snapshot, cq, classify, run_census

SHA = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
CTOLS = [1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14]


def read_fetched(path):
    dec = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        dec.append([t.strip() for t in line.replace(" ", ",").split(",") if t.strip()])
    return dec


def structure_record(Xn, d, n):
    G = Xn @ Xn.T
    np.fill_diagonal(G, -2.0)
    M = float(G.max())
    iu = np.triu_indices(n, 1)
    gaps = M - np.sort(G[iu])[::-1]
    return {"d": d, "n": n, "M": M,
            "pairs_within": {f"{t:.0e}": int((gaps <= t).sum()) for t in CTOLS},
            "structs": [classify.structure(Xn, d, G, M, t) for t in CTOLS]}, G


def main():
    fdir, sel, out = sys.argv[1], sys.argv[2], sys.argv[3]
    sh, nsh = 0, 1
    if "--shard" in sys.argv:
        sh, nsh = (int(x) for x in sys.argv[sys.argv.index("--shard") + 1].split("/"))
    os.makedirs(out, exist_ok=True)
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    todo = []
    for line in open(sel):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        todo.append((int(p[0]), int(p[1]), p[2], p[3] if len(p) > 3 else ""))
    for k, (d, n, cls5, why) in enumerate(todo):
        if k % nsh != sh:
            continue
        js = os.path.join(out, f"d{d}_n{n}.json")
        src = os.path.join(fdir, f"d{d}_n{n}.txt")
        if os.path.exists(js) or not os.path.exists(src):
            continue
        t0 = time.process_time()
        rec = {"d": d, "n": n, "census_class": cls5, "selected_because": why, "program_sha256": SHA,
               "fetched_sha256": hashlib.sha256(open(src, "rb").read()).hexdigest(), "attempts": []}
        try:
            dec = read_fetched(src)
            if len(dec) != n or any(len(r) != d for r in dec):
                raise ValueError(f"file has {len(dec)} rows")
            X = np.array([[float(x) for x in r] for r in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            srec, G = structure_record(X, d, n)
            rec["max_cos_double"] = srec["M"]
            # the archived code, for the comparison (untrusted)
            pts = snapshot.read_points(XC.SNAP, idx[(d, n)])
            Y = np.array([[float(x) for x in p] for p in pts])
            Y = Y / np.linalg.norm(Y, axis=1, keepdims=True)
            H = Y @ Y.T
            np.fill_diagonal(H, -2.0)
            rec["archived_max_cos_double"] = float(H.max())
            iu = np.triu_indices(n, 1)
            a, b = np.sort(G[iu]), np.sort(H[iu])
            rec["dist_lower_point"] = float(np.max(np.abs(a - b)) / 2)
            del G, H
            o, f, sk = run_census.plan(srec, run_census.CAP)
            rec["skipped_above_cap"] = [{"kind": x, "tol": t, "size": z} for x, t, z in sk]
            done = False
            for kind, t, z in o:
                lg = cq.certify(X, d, n, tol=t, dec=dec)
                lg["planned_kind"] = kind
                rec["attempts"].append({x: lg.get(x) for x in ("status", "kind", "tol", "unknowns", "contacts",
                                                               "dist_ls", "planned_kind")})
                if lg.get("status") in ("PASS", "LOCALMIN", "NONOPT"):
                    ex = {}
                    XC.check_entry(X, dec, d, n, t, lg.get("status"), ex)
                    XC._ABS.clear()
                    rec["attempts"][-1]["exact"] = {x: ex.get(x) for x in (
                        "status", "dist_ls", "c_interval_30", "gap_lower", "stress_min_rel", "tube", "rho",
                        "code_inside_rho", "unknowns", "contacts", "E_sha256")}
                if lg.get("status") in ("PASS", "LOCALMIN"):
                    done = True
                    break
            if not done:
                for kind, t, z in f:
                    lg = cq.certify_flex(X, d, n, tol=t, dec=dec, max_rows=run_census.CAP)
                    rec["attempts"].append({x: lg.get(x) for x in ("status", "kind", "tol", "contacts",
                                                                   "points_in_contacts", "flex_radius_l1_stored")})
                    if lg.get("status") == "FLEX":
                        break
        except Exception:
            rec["error"] = traceback.format_exc()[-1500:]
        ex_st = [a.get("exact", {}).get("status") for a in rec["attempts"] if a.get("exact")]
        st = [a.get("status") for a in rec["attempts"]]
        first = next((s for s in ("PASS", "LOCALMIN", "FLEX", "NONOPT") if s in st), "none")
        rec["status_first_pass"] = first
        rec["status_exact"] = ex_st[-1] if ex_st else None
        rec["cpu_seconds"] = round(time.process_time() - t0, 2)
        with open(js + ".part", "w") as fo:
            json.dump(rec, fo, default=str)
        os.replace(js + ".part", js)
        print(d, n, cls5, first, rec["status_exact"], rec["cpu_seconds"], flush=True)


if __name__ == "__main__":
    main()
