"""Two scans of the positive entries (census version 8).

  nice -n 10 python3 -I tolscan.py other OUT.jsonl [--shard i/k]   every positive entry, every other tolerance with
                                                                  another proposal
  nice -n 10 python3 -I tolscan.py close OUT.jsonl [--shard i/k]   the certificates whose gap is below 1e-9: the next
                                                                  pairs added

'other': the threshold proposals are nested in tau, so a tolerance proposes another contact graph exactly when the
number of contacts or of non-rattlers differs from that of the certified tolerance (census/classify.jsonl). For each
such tolerance at which the proposal is isostatic or hyperstatic, the census certificate (cq.certify) is run, and a
PASS, LOCALMIN or NONOPT is decided again exactly (verify/exactcheck2.py).
'close': for each PASS or LOCALMIN certificate whose gap (first pass, rattlers included) is below 1e-9, the tolerance is
raised to include the next pair or pairs (the pairs within 10 times the gap of the next one), and the certificate is run
on that proposal and decided exactly.
Every certificate decided here is reported; none changes a verdict of the census.
"""
import json, os, sys, time, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "verify"))
import exactcheck2 as XC
import numpy as np
import snapshot, cq

TOLS = [1e-10, 1e-11, 1e-12, 1e-13, 1e-9, 1e-8, 1e-7]


def load(d, n, idx):
    pts = snapshot.read_points(XC.SNAP, idx[(d, n)])
    dec = [[t.strip() for t in p] for p in pts]
    X = np.array([[float(x) for x in p] for p in dec])
    return X / np.linalg.norm(X, axis=1, keepdims=True), dec


def attempt(X, dec, d, n, tol):
    lg = cq.certify(X, d, n, tol, dec=dec)
    rec = {"tol": tol, "first_pass": lg.get("status"), "kind": lg.get("kind"), "unknowns": lg.get("unknowns"),
           "contacts": lg.get("contacts")}
    if lg.get("status") in ("PASS", "LOCALMIN", "NONOPT"):
        ex = {}
        XC.check_entry(X, dec, d, n, tol, lg.get("status"), ex)
        XC._ABS.clear()
        rec["exact"] = {k: ex.get(k) for k in ("status", "c_interval_30", "dist_ls", "gap_lower", "E_sha256", "contacts",
                                               "unknowns", "tube", "stress_min_rel")}
    return rec


def main():
    mode, outp = sys.argv[1], sys.argv[2]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    cls = {}
    for line in open(os.path.join(XC.CERT, "census", "classify.jsonl")):
        r = json.loads(line)
        cls[(r["d"], r["n"])] = r
    done = set()
    if os.path.exists(outp):
        done = {(r["d"], r["n"]) for r in map(json.loads, open(outp))}
    import glob
    todo = []
    for f in sorted(glob.glob(os.path.join(XC.CERT, "census", "results_v4", "*.json"))):
        j = json.load(open(f))
        if j["status"] not in ("PASS", "LOCALMIN"):
            continue
        a = [x for x in j["attempts"] if x.get("status") == j["status"]][0]
        k = (j["d"], j["n"])
        if mode == "other":
            s = {x["tol"]: x for x in cls[k]["structs"]}
            t0 = a["tol"]
            others = [t for t in TOLS if t != t0 and (s[t]["m"], s[t]["nV"]) != (s[t0]["m"], s[t0]["nV"])
                      and s[t]["nV"] > 0 and s[t]["m"] >= s[t]["iso"]
                      and s[t]["nV"] * j["d"] + 1 <= 3100]
            # one tolerance per distinct proposal
            seen, uniq = set(), []
            for t in others:
                key = (s[t]["m"], s[t]["nV"])
                if key not in seen:
                    seen.add(key)
                    uniq.append(t)
            if uniq:
                todo.append((k, j["status"], t0, uniq))
        else:
            g = (a.get("gaps") or {}).get("min_gap_float")
            if g is not None and g < 1e-9:
                todo.append((k, j["status"], a["tol"], None))
    if "--shard" in sys.argv:
        i_, k_ = map(int, sys.argv[sys.argv.index("--shard") + 1].split("/"))
        todo = todo[i_::k_]
    print(len(todo), "entries", flush=True)
    for k, st, t0, extra in todo:
        if k in done:
            continue
        d, n = k
        rec = {"d": d, "n": n, "census_status": st, "census_tol": t0}
        t1 = time.process_time()
        try:
            X, dec = load(d, n, idx)
            if mode == "other":
                rec["attempts"] = [attempt(X, dec, d, n, t) for t in extra]
            else:
                G = X @ X.T
                np.fill_diagonal(G, -2.0)
                M = float(G.max())
                v = np.sort(G[np.triu_indices(n, 1)])[::-1]
                m0 = int((v >= M - t0).sum())          # the certified contacts
                gap = M - v[m0]
                # the next pairs: those within 10 times the gap of the next one
                m1 = int((v >= v[m0] - 10 * gap).sum())
                tol = (M - v[m1 - 1] + M - v[m1]) / 2 if m1 < len(v) else M - v[-1] + 1e-12
                rec.update(next_gap=float(gap), added_pairs=m1 - m0, tol_close=float(tol))
                rec["attempts"] = [attempt(X, dec, d, n, float(tol))]
        except Exception:
            rec["error"] = traceback.format_exc()[-1200:]
        rec["cpu_seconds"] = round(time.process_time() - t1, 1)
        with open(outp, "a") as fo:
            fo.write(json.dumps(rec, default=str) + "\n")
        print(d, n, st, [(a_["tol"], a_["first_pass"], (a_.get("exact") or {}).get("status")) for a_ in rec.get("attempts", [])],
              rec["cpu_seconds"], flush=True)


if __name__ == "__main__":
    main()
