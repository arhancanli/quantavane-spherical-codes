"""Run sphcert.certify on a list of codes, one after the other, and write one JSON file per code.

  python3 -I certify/run_batch.py LIST OUTDIR [--shard I/K] [--max-unknowns M] [--tol T]
LIST: lines "path d n" (paths relative to the repository root). A code whose JSON file exists in OUTDIR is skipped,
so a stopped run can be restarted. With --shard I/K only every K-th line from the I-th (0-based) is run.
--retry DIR: only the codes whose result in DIR is not PASS, with the contact tolerance taken, in this order, from
1e-11, 1e-12, 1e-13, 1e-9, 1e-8, 1e-7 (each that makes the system square, until one passes; the others are skipped);
a code that failed a check with 1e-10 is first run again with 1e-10. The certificate does not depend on how the contacts were chosen: whatever is proved is proved.
--duty F (default 0.15): the certifier sleeps between its steps so that its CPU time is at most the fraction F of the wall
time, so that it stays a light background job on a loaded machine.
"""
import json, os, sys, time, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sphcert

ROOT = os.path.dirname(HERE)


def main():
    lst, out = sys.argv[1], sys.argv[2]
    shard = (0, 1)
    maxu, tol, duty, retry = 10 ** 9, 1e-10, 0.15, None
    args = sys.argv[3:]
    for i, a in enumerate(args):
        if a == "--shard":
            shard = tuple(int(x) for x in args[i + 1].split("/"))
        if a == "--max-unknowns":
            maxu = int(args[i + 1])
        if a == "--tol":
            tol = float(args[i + 1])
        if a == "--duty":
            duty = float(args[i + 1])
        if a == "--retry":
            retry = args[i + 1]
    os.makedirs(out, exist_ok=True)
    sphcert.PACE = sphcert.Pace(duty)  # light steps inside each code too (see sphcert.Pace)
    lines = [l.split() for l in open(lst) if l.strip() and not l.startswith("#")]
    for k, (path, d, n) in enumerate(lines):
        if k % shard[1] != shard[0]:
            continue
        d, n = int(d), int(n)
        name = os.path.basename(path).replace(".txt", "")
        js = os.path.join(out, name + ".json")
        if os.path.exists(js):
            continue
        if n * d > maxu:
            print(f"{name}: skipped (n*d = {n * d} > {maxu})", flush=True)
            continue
        tols = [tol]
        if retry is not None:
            prev = os.path.join(retry, name + ".json")
            if not os.path.exists(prev) or json.load(open(prev))["status"] == "PASS":
                continue
            X = sphcert.read_code(os.path.join(ROOT, path), d, n)
            tols = []
            first = json.load(open(prev))["status"]
            for t_ in ((1e-10,) if first.startswith("FAIL") else ()) + (1e-11, 1e-12, 1e-13, 1e-9, 1e-8, 1e-7):
                sphcert.PACE.tick()
                V, E, _, _ = sphcert.contact_structure(X, d, t_)
                if len(V) and len(E) == len(V) * (d - 1) - d * (d - 1) // 2 + 1:
                    tols.append(t_)
            if not tols:
                log = {"path": path, "status": "no tolerance makes the system square", "wall": 0.0}
                with open(js, "w") as f:
                    json.dump(log, f, indent=1)
                print(f"{name}: {log['status']}", flush=True)
                continue
        t, cpu = time.time(), time.process_time()
        for t_ in tols:
            log = {"path": path}
            try:
                sphcert.certify(os.path.join(ROOT, path), d, n, tol=t_, log=log)
            except Exception as e:  # recorded, never hidden
                log["status"] = "ERROR " + repr(e)
                log["traceback"] = traceback.format_exc()
            if log["status"] == "PASS":
                break
        if retry is not None:
            log["retry_tolerances"] = tols
        log["wall"] = round(time.time() - t, 2)
        with open(js + ".tmp", "w") as f:
            json.dump(log, f, indent=1, default=str)
        os.replace(js + ".tmp", js)
        print(f"{name}: {log['status']} ({log['wall']} s)", flush=True)
        time.sleep(1.0)


if __name__ == "__main__":
    main()
