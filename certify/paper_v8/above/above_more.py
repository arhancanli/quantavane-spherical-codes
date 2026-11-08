"""A second pass over the structured codes above the size limit for which run_above.py found no strictly positive stress
(census version 8): run_above.py, imported unchanged, with a larger cap on the number of classes of contacts.

  nice -n 10 python3 -I above_more.py ABOVE_OUTDIR OUTDIR [--max-classes 2000] [--mem-mb 1500] [--seconds 900]

For every record of ABOVE_OUTDIR (run_above.py) whose status is 'undecided: no ... stress ...': run_above.run_entry in a
child process, with the cap of run_above.class_stress raised from 300 to --max-classes, the same field K (the record's
D) and the same column limit. The cap only decides which UNTRUSTED class structures are tried; everything that enters a
verdict is run_above.py's code and checks (the exact core, the rank over F_p, the exact stress check). Above the 4,000
columns up to which run_above.py tries a linear-programming stress, a bar-rigid core without a stress also gets one here
(lp_more). Each child is stopped above --mem-mb of resident memory or after --seconds.
"""
import hashlib, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SHA = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()


def _child(conn, d, n, disc, max_cols, max_classes):
    import traceback
    import run_above as RA
    RA.class_stress.__defaults__ = (None, max_classes)        # (target, max_classes)
    log = {}
    try:
        RA.run_entry(d, n, disc, max_cols, log)
        st = str(log.get("status", ""))
        if st.startswith("undecided") and "stress" in st and log.get("flex_space_exact") and len(log.get("components", [])) == 1 \
                and int(log.get("cols", 0)) > 4000:
            lp_more(d, n, int(log.get("disc", disc)), log)
    except MemoryError:
        log["status"] = "stopped: MemoryError"
    except Exception:
        log["status"] = "error"
        log["error"] = traceback.format_exc()[-1500:]
    log.pop("_proj_w", None)
    log.pop("_class_stress", None)
    log["run_above_sha256"] = RA.SHA
    conn.send(json.dumps(log, default=str))
    conn.close()


def lp_more(d, n, disc, log):
    """For a bar-rigid core (proved by run_entry in this process) without a stress, above the 4,000 columns up to which
    run_above.py tries it: the max-min stress of the linear program (UNTRUSTED, v5/rigidity.py lp_stress) and, if it is
    positive, a stress constant on the classes of its values, solved and checked exactly by run_above.class_stress."""
    import numpy as np
    import core
    import rigidity as RG
    import run_above as RA
    C = core.build(RA.SNAP, d, n, disc)
    if C.get("status") != "ok" or len(C["E"]) != int(log.get("contacts", -1)):
        log["lp_more"] = "core not reproduced"
        return
    V, E, W = C["V"], C["E"], C["W"]
    wf, tmax = RG.lp_stress(C["Xf"][V], E, d)
    log["lp_more_maxmin_stress"] = tmax
    if wf is None or not tmax > 1e-9:
        log["lp_more"] = "the linear program finds no strictly positive stress (untrusted)"
        return
    key, target = RA.cluster(np.asarray(wf), 1e-7)
    om = RA.class_stress(W, E, d, len(V), log, key, "lp_more", target)
    if om is not None:
        log.pop("_class_stress", None)
        log["stress"] = "constant on the classes of a linear-programming stress, solved exactly (second pass)"
        log["status"] = "STRICT (bar-rigid core with a strictly positive stress; Theorem 2.4)"


def main():
    import multiprocessing, subprocess
    src, outd = sys.argv[1], sys.argv[2]
    args = sys.argv[3:]
    max_classes = int(args[args.index("--max-classes") + 1]) if "--max-classes" in args else 2000
    mem = int(args[args.index("--mem-mb") + 1]) if "--mem-mb" in args else 1500
    secs = int(args[args.index("--seconds") + 1]) if "--seconds" in args else 900
    os.makedirs(outd, exist_ok=True)
    todo = []
    for f in sorted(os.listdir(src)):
        if f.endswith(".json"):
            j = json.load(open(os.path.join(src, f)))
            st = str(j.get("status", ""))
            if st.startswith("undecided") and "stress" in st:
                todo.append((int(j.get("cols") or 0), j["d"], j["n"], int(j.get("disc", 1))))
    todo.sort()
    print(len(todo), "entries", flush=True)
    for cols, d, n, disc in todo:
        js = os.path.join(outd, f"d{d}_n{n}.json")
        if os.path.exists(js):
            continue
        log = {"d": d, "n": n, "disc": disc, "program_sha256": SHA, "max_classes": max_classes}
        t0 = time.time()
        ctx = multiprocessing.get_context("spawn")
        rd, wr = ctx.Pipe(duplex=False)
        pr = ctx.Process(target=_child, args=(wr, d, n, disc, 20000, max_classes))
        pr.start()
        wr.close()
        peak, got, why = 0, None, None
        while True:
            if rd.poll(1.0):
                got = json.loads(rd.recv())
                break
            if not pr.is_alive():
                break
            o = subprocess.run(["ps", "-o", "rss=", "-p", str(pr.pid)], capture_output=True, text=True).stdout.strip()
            rss = int(o) // 1024 if o else 0
            peak = max(peak, rss)
            if rss > mem:
                why = f"stopped: memory above {mem} MB"
                pr.kill()
                break
            if time.time() - t0 > secs:
                why = f"stopped: time above {secs} s"
                pr.kill()
                break
        pr.join()
        if got is None:
            got = {"status": why or "stopped: child process ended"}
        log.update(got)
        log["peak_rss_mb_sampled"] = peak
        log["wall"] = round(time.time() - t0, 1)
        with open(js + ".part", "w") as f:
            json.dump(log, f)
        os.replace(js + ".part", js)
        print(d, n, log.get("status"), log["wall"], flush=True)


if __name__ == "__main__":
    main()
