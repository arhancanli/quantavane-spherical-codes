"""Gap-based proposals for the entries without a positive verdict (census version 7; the referee's question 1).

  nice -n 10 python3 -I gapcut.py scan OUT.jsonl [--shard I/K]
  nice -n 10 python3 -I gapcut.py certify SAMPLE.tsv OUT.jsonl

scan (untrusted, a proposal only): for every entry whose census record has status FLEX or failed, the cosines of the
double-precision code sorted downwards, c_1 >= c_2 >= ...; for a cut after k pairs, the spread s_k = c_1 - c_k and the gap
g_k = c_k - c_(k+1). Among the cuts with s_k <= 1e-4, taken in order of decreasing g_k / max(s_k, 1e-15), the first one
whose core (points with at least d proposed contacts, removed repeatedly) is isostatic or hyperstatic is recorded: the
most cleanly separated proposal that the certificates could take. Nothing is proved here.

certify: for the entries of SAMPLE.tsv, the census certificate (cq.certify, isostatic or hyperstatic as the class says)
at a tolerance inside the gap, tol = s_k + g_k / 2; every PASS or LOCALMIN is then decided again by exactcheck.py.
"""
import json, os, sys, time, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "verify"))
import exactcheck as XC
import numpy as np
import snapshot

S_MAX = 1e-4
NCAND = 40


def core_class(n, d, I_, J_):
    alive = np.ones(n, bool)
    I_, J_ = np.asarray(I_), np.asarray(J_)
    while True:
        keep = alive[I_] & alive[J_]
        deg = np.bincount(I_[keep], minlength=n) + np.bincount(J_[keep], minlength=n)
        drop = alive & (deg < d)
        if not drop.any():
            break
        alive &= ~drop
    nV = int(alive.sum())
    m = int((alive[I_] & alive[J_]).sum())
    iso = nV * (d - 1) - d * (d - 1) // 2 + 1 if nV else 0
    cls = "empty" if nV == 0 else ("iso" if m == iso else ("hyper" if m > iso else "hypo"))
    return cls, nV, m, iso


def scan_entry(X, d, n):
    G = X @ X.T
    iu = np.triu_indices(n, 1)
    g = G[iu]
    K = min(len(g), 6 * n * d)
    top = np.argpartition(-g, K - 1)[:K] if K < len(g) else np.arange(len(g))
    top = top[np.argsort(-g[top], kind="stable")]
    c = g[top]
    s = c[0] - c
    gaps = np.empty(len(c))
    gaps[:-1] = c[:-1] - c[1:]
    gaps[-1] = 0.0
    ok = np.nonzero((s <= S_MAX) & (gaps > 0))[0]
    if len(ok) == 0:
        return {"found": False}
    ratio = gaps[ok] / np.maximum(s[ok], 1e-15)
    order = ok[np.argsort(-ratio)][:NCAND]
    I_all, J_all = iu[0][top], iu[1][top]
    for k0 in order:
        k = int(k0) + 1                                  # cut after k pairs
        cls, nV, m, iso = core_class(n, d, I_all[:k], J_all[:k])
        if cls in ("iso", "hyper"):
            return {"found": True, "k": k, "spread": float(s[k - 1]), "gap": float(gaps[k - 1]),
                    "ratio": float(gaps[k - 1] / max(s[k - 1], 1e-15)), "class": cls, "nV": nV, "m": m, "iso": iso,
                    "unknowns": nV * d + 1, "tol": float(s[k - 1] + gaps[k - 1] / 2)}
    return {"found": False}


def load(d, n, idx):
    pts = snapshot.read_points(XC.SNAP, idx[(d, n)])
    dec = [[t.strip() for t in p] for p in pts]
    X = np.array([[float(x) for x in p] for p in dec])
    return X / np.linalg.norm(X, axis=1, keepdims=True), dec


def main():
    mode = sys.argv[1]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    if mode == "scan":
        outp = sys.argv[2]
        sh, nsh = 0, 1
        if "--shard" in sys.argv:
            sh, nsh = (int(x) for x in sys.argv[sys.argv.index("--shard") + 1].split("/"))
        done = set()
        if os.path.exists(outp):
            done = {(r["d"], r["n"]) for r in map(json.loads, open(outp))}
        import glob
        files = sorted(glob.glob(os.path.join(XC.CERT, "census", "results_v4", "*.json")))
        for k, f in enumerate(files):
            if k % nsh != sh:
                continue
            j = json.load(open(f))
            if j["status"] not in ("FLEX", "failed") or (j["d"], j["n"]) in done:
                continue
            d, n = j["d"], j["n"]
            rec = {"d": d, "n": n, "census_status": j["status"]}
            try:
                X, _ = load(d, n, idx)
                rec.update(scan_entry(X, d, n))
            except Exception:
                rec["error"] = traceback.format_exc()[-800:]
            with open(outp, "a") as fo:
                fo.write(json.dumps(rec) + "\n")
    elif mode == "certify":
        import cq
        sample, outp = sys.argv[2], sys.argv[3]
        done = set()
        if os.path.exists(outp):
            done = {(r["d"], r["n"]) for r in map(json.loads, open(outp))}
        for line in open(sample):
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            d, n, tol = int(p[0]), int(p[1]), float(p[2])
            if (d, n) in done:
                continue
            rec = {"d": d, "n": n, "tol": tol, "stratum": p[3] if len(p) > 3 else ""}
            t0 = time.process_time()
            try:
                X, dec = load(d, n, idx)
                lg = cq.certify(X, d, n, tol, dec=dec)
                rec["first_pass"] = lg.get("status")
                rec["kind"] = lg.get("kind")
                if lg.get("status") in ("PASS", "LOCALMIN", "NONOPT"):
                    ex = {}
                    XC.check_entry(X, dec, d, n, tol, lg.get("status"), ex)
                    _ = XC._ABS.clear()
                    rec["exact"] = {k: ex.get(k) for k in ("status", "dist_ls", "c_interval_30", "gap_lower",
                                                           "stress_min_rel", "tube", "unknowns", "contacts")}
            except Exception:
                rec["error"] = traceback.format_exc()[-800:]
            rec["cpu_seconds"] = round(time.process_time() - t0, 2)
            with open(outp, "a") as fo:
                fo.write(json.dumps(rec, default=str) + "\n")
            print(d, n, rec.get("first_pass"), (rec.get("exact") or {}).get("status"), rec["cpu_seconds"], flush=True)


if __name__ == "__main__":
    main()
