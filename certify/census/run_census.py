"""Run the census certifier (cq.py, version 4) over the archived table.

  python3 -I certify/census/run_census.py SNAPDIR CLASSIFY.jsonl OUTDIR --shard I/K [--cap 3100]

The contacts and rattlers proposed at a tolerance tau (classify.jsonl; untrusted) define the class of an entry at
tau: isostatic, hyperstatic or hypostatic (|E| equal to, above or below |V|(d-1) - d(d-1)/2 + 1 on the non-rattlers
V), or empty (V is empty). One rule for every entry (paper, Section 4.2), with TOLS = 1e-10, 1e-11, 1e-12, 1e-13,
1e-9, 1e-8, 1e-7:
 1. the isostatic certificate at every tau of TOLS where the class is isostatic, in this order; then the hyperstatic
    certificate at every tau where it is hyperstatic, in this order; until one gives PASS or LOCALMIN (a NONOPT is
    recorded and the search goes on);
 2. if none gives PASS or LOCALMIN: the flexibility certificate (all pairs within tau, no rattlers removed) at the
    largest tau, from 1e-7 down to 1e-13, at which the class is hypostatic or empty, until one gives FLEX (it then
    holds at every smaller tau as well: the flexibility proposition of the paper).
An attempt is made only if the system it solves has at most CAP unknowns (the square system of E or E': |V| d + 1;
the flexibility test: the number of its rows, the points in a pair of E plus |E|). Whatever passes is proved,
whichever tolerance proposed it. One JSON file per entry in OUTDIR, stamped with the SHA-256 of the certifier and of
this driver; a stopped run restarts where it stopped.

With --mem-limit-mb M, every hyperstatic attempt runs in a child process whose resident memory is checked every
second; above M megabytes the child is stopped and the attempt recorded as "stopped: memory above M MB" (the
linear program that chooses the combined stress of a hyperstatic graph with many extra contacts can need more than
2 GB). The certifier is the same; only where it runs changes.
"""
import hashlib, json, os, sys, time, traceback, warnings

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
warnings.filterwarnings("ignore", message=".*encountered in matmul")  # Accelerate's spurious flags (untrusted calls)
import cq
import snapshot
import numpy as np

TOLS = [1e-10, 1e-11, 1e-12, 1e-13, 1e-9, 1e-8, 1e-7]
FLEX_TOLS = [1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13]
CAP = 3100
DRIVER_SHA256 = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()


def klass(s):
    """The class at one tolerance, from a classify.jsonl structure record."""
    if s["nV"] == 0:
        return "empty"
    if s["m"] == s["iso"]:
        return "iso"
    return "hyper" if s["m"] > s["iso"] else "hypo"


def plan(rec, cap=CAP):
    """The attempts for one entry: a list of (certificate, tau, size) in the order of the rule; and the attempts
    skipped because their system is above the cap."""
    d, n = rec["d"], rec["n"]
    s = {x["tol"]: x for x in rec["structs"]}
    opt, flex, skipped = [], [], []
    for kind in ("iso", "hyper"):
        for t in TOLS:
            if klass(s[t]) == kind:
                size = s[t]["nV"] * d + 1
                (opt if size <= cap else skipped).append((kind, t, size))
    for t in FLEX_TOLS:
        if klass(s[t]) in ("hypo", "empty"):
            m = rec["pairs_within"][f"{t:.0e}"]
            size = min(n, 2 * m) + m          # an upper bound for the rows; certify_flex applies the cap exactly
            (flex if size <= cap else skipped).append(("flex", t, size))
    return opt, flex, skipped


def _child(conn, X, d, n, t, dec):
    """A hyperstatic attempt in a child process (see --mem-limit-mb)."""
    log = cq.certify(X, d, n, tol=t, dec=dec)
    conn.send(json.dumps(log, default=str))
    conn.close()


def certify_watched(X, d, n, t, dec, limit_mb):
    """cq.certify in a child process, stopped when its resident memory exceeds limit_mb."""
    import multiprocessing, subprocess
    ctx = multiprocessing.get_context("spawn")
    rd, wr = ctx.Pipe(duplex=False)
    p = ctx.Process(target=_child, args=(wr, X, d, n, t, dec))
    p.start()
    wr.close()
    peak, log = 0, None
    while True:
        if rd.poll(1.0):
            log = json.loads(rd.recv())
            break
        if not p.is_alive():
            break
        out = subprocess.run(["ps", "-o", "rss=", "-p", str(p.pid)], capture_output=True, text=True).stdout.strip()
        rss = int(out) // 1024 if out else 0
        peak = max(peak, rss)
        if rss > limit_mb:
            p.kill()
            break
    p.join()
    if log is None:
        log = {"d": d, "n": n, "tol": t, "kind": "hyperstatic", "code_sha256": cq.CODE_SHA256,
               "status": f"stopped: memory above {limit_mb} MB" if peak > limit_mb else "stopped: child process ended"}
    log["peak_rss_mb_sampled"] = peak
    return log


def main():
    snapdir, cls, out = sys.argv[1], sys.argv[2], sys.argv[3]
    args = sys.argv[4:]
    opt = {"--shard": "0/1", "--cap": str(CAP), "--mem-limit-mb": "0"}
    for i, a in enumerate(args):
        if a in opt:
            opt[a] = args[i + 1]
    sh, nsh = (int(x) for x in opt["--shard"].split("/"))
    cap = int(opt["--cap"])
    mem = int(opt["--mem-limit-mb"])
    os.makedirs(out, exist_ok=True)
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    todo = []
    for line in open(cls):
        r = json.loads(line)
        if not idx[(r["d"], r["n"])]["listed"]:
            continue
        o, f, sk = plan(r, cap)
        size = min([x[2] for x in o + f] or [0])
        todo.append((size, r["d"], r["n"], o, f, sk))
    todo.sort()
    for k, (size, d, n, o, f, sk) in enumerate(todo):
        if k % nsh != sh:
            continue
        js = os.path.join(out, f"d{d}_n{n}.json")
        if os.path.exists(js):
            continue
        t0 = time.time()
        res = {"d": d, "n": n, "attempts": [], "skipped_above_cap": [{"kind": a, "tol": t, "size": z} for a, t, z in sk],
               "cap": cap, "mem_limit_mb": mem, "code_sha256": dict(cq.CODE_SHA256, **{"run_census.py": DRIVER_SHA256}),
               "certifier_version": cq.VERSION}
        try:
            pts = snapshot.read_points(snapdir, idx[(d, n)])
            dec = [[t.strip() for t in p] for p in pts]
            X = np.array([[float(x) for x in p] for p in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            done = False
            for kind, t, z in o:
                if kind == "hyper" and mem > 0:
                    log = certify_watched(X, d, n, t, dec, mem)
                else:
                    log = cq.certify(X, d, n, tol=t, dec=dec)
                log["planned_kind"] = kind
                res["attempts"].append(log)
                if log.get("status") in ("PASS", "LOCALMIN"):
                    done = True
                    break
            if not done:
                for kind, t, z in f:
                    log = cq.certify_flex(X, d, n, tol=t, dec=dec, max_rows=cap)
                    res["attempts"].append(log)
                    if log.get("status") == "FLEX":
                        break
        except MemoryError:
            res["error"] = "MemoryError"
        except Exception:
            res["error"] = traceback.format_exc()[-2000:]
        st = [a.get("status") for a in res["attempts"]]
        res["status"] = next((s for s in ("PASS", "LOCALMIN", "FLEX", "NONOPT") if s in st),
                             "failed" if res["attempts"] else ("not attempted" if not res.get("error") else "error"))
        res["wall"] = round(time.time() - t0, 2)
        with open(js + ".part", "w") as fo:
            json.dump(res, fo, default=str)
        os.replace(js + ".part", js)
        print(f"d{d}_n{n} {size} {res['status']} {res['wall']}s", flush=True)


if __name__ == "__main__":
    main()
