"""Fetch the current coordinates of selected entries of the live table (https://spherical-codes.org/data/D/N).

  nice -n 10 python3 -I fetch_live.py SELECT.tsv OUTDIR

One request at a time, at most one every 1.5 seconds; a file already present is not fetched again; five failures in a
row stop the run. OUTDIR/manifest.tsv gets d, n, SHA-256, bytes, the UTC time and the HTTP status of every fetch.
The files are third-party data: they are kept locally and are not redistributed.
"""
import hashlib, os, sys, time, urllib.request
from datetime import datetime, timezone

UA = "quantavane-census-currency-check/7 (research; sequential, rate-limited)"


def main():
    sel, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    man = os.path.join(out, "manifest.tsv")
    todo = []
    for line in open(sel):
        if line.startswith("#"):
            continue
        p = line.split("\t")
        todo.append((int(p[0]), int(p[1])))
    fails = 0
    for d, n in todo:
        path = os.path.join(out, f"d{d}_n{n}.txt")
        if os.path.exists(path):
            continue
        url = f"https://spherical-codes.org/data/{d}/{n}"
        t0 = time.time()
        status, body = None, None
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as r:
                status, body = r.status, r.read()
        except Exception as e:
            status = f"error {type(e).__name__}: {e}"
        now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        if status == 200 and body:
            with open(path + ".part", "wb") as f:
                f.write(body)
            os.replace(path + ".part", path)
            sha = hashlib.sha256(body).hexdigest()
            fails = 0
            with open(man, "a") as f:
                f.write(f"{d}\t{n}\t{sha}\t{len(body)}\t{now}\t200\n")
        else:
            fails += 1
            with open(man, "a") as f:
                f.write(f"{d}\t{n}\t-\t0\t{now}\t{status}\n")
            if fails >= 5:
                print("five failures in a row: stopping", flush=True)
                break
        time.sleep(max(0.0, 1.5 - (time.time() - t0)))
    print("done", flush=True)


if __name__ == "__main__":
    main()
