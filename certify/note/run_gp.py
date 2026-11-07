# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""Run a GP script of this directory with the PARI library (cypari) and record what it printed.

  python3 -I certify/note/run_gp.py SCRIPT.gp OUT_BASENAME [FILES_TO_REMOVE_FIRST ...]

The script is executed by PARI's own interpreter (read()), with DATA set to the data/ directory (or to the directory
in the environment variable QV_NOTE_DATA, e.g. data_xprime/ for the second code of Section 7); nothing in this runner
computes anything. Standard output of the run goes to out/OUT_BASENAME.log; the lines "KEY = value" are collected in
out/OUT_BASENAME.json together with the PARI version, the wall time, the SHA-256 of the script (taken before the
run) and of every file in data/ before and after the run. The same scripts run unchanged in the gp program (gp -q, after DATA = "path/to/data/").
"""
import hashlib, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [p for p in os.environ.get("QV_LIBS", "").split(os.pathsep) if p]   # optional library folders


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    script, base = sys.argv[1], sys.argv[2]
    data = os.environ.get("QV_NOTE_DATA") or os.path.join(HERE, "data")     # data_xprime/ for the code X'
    data = os.path.abspath(data)
    for f in sys.argv[3:]:
        p = os.path.join(data, f)
        if os.path.exists(p):
            os.remove(p)
    os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
    logp = os.path.join(HERE, "out", base + ".log")
    script_sha = sha(os.path.join(HERE, script))     # hashed before the run
    data_sha_before = {f: sha(os.path.join(data, f)) for f in sorted(os.listdir(data))}
    from cypari import pari
    pari.allocatemem(5 * 10 ** 8, 4 * 10 ** 9)          # stack 0.5 GB, may grow to 4 GB
    version = str(pari("Str(version())"))
    t0 = time.time()
    sys.stdout.flush()
    fd = os.open(logp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    saved = os.dup(1)
    os.dup2(fd, 1)
    err = None
    try:
        pari('DATA = "%s/"' % data)
        pari('read("%s")' % os.path.join(HERE, script))
    except Exception as e:      # a GP error() ends the run; it is recorded as a failure
        err = repr(e) + ": " + (e.errtext() if hasattr(e, "errtext") else "")
    finally:
        sys.stdout.flush()
        os.dup2(saved, 1)
        os.close(fd)
        os.close(saved)
    secs = round(time.time() - t0, 1)
    vals = {}
    for line in open(logp):
        if " = " in line:
            k, v = line.rstrip("\n").split(" = ", 1)
            vals[k.strip()] = v
    rec = {"script": script, "script_sha256": script_sha,
           "runner_sha256": sha(os.path.abspath(__file__)), "pari_version": version, "seconds": secs,
           "error": err, "data_sha256_before_run": data_sha_before,
           "data_sha256_after_run": {f: sha(os.path.join(data, f)) for f in sorted(os.listdir(data))},
           "results": vals}
    json.dump(rec, open(os.path.join(HERE, "out", base + ".json"), "w"), indent=1)
    print(json.dumps({"script": script, "seconds": secs, "error": err, "lines": len(vals)}))
    if err:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
