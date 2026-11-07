"""Writes the two lists of codes that run_batch.py certifies (sorted by N * d).

  python3 -I certify/make_lists.py
lists/uploaded.txt: every code in codes/ (the last code we uploaded for each entry).
lists/converged.txt: the codes that observations/structure.py measures (runs of the active-set polish finish.py whose
summary line says polish=converged, finish lane after 1 Nov 01:35 and fast lane after 1 Nov 00:30; 433 entries), with
the same selection rule, restricted to the 377 whose output files are in this repository (the other 56 are the table's
codes polished in the last digits, which observations/README.md explains are not redistributed).
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "observations"))
from common import lane_runs, output, FAST_SPLIT  # noqa: E402

FIN_SPLIT = 1793482500  # 2026-11-01 01:35 +04, as in observations/structure.py


def main():
    os.makedirs(os.path.join(HERE, "lists"), exist_ok=True)
    up = []
    for f in os.listdir(os.path.join(ROOT, "codes")):
        m = re.fullmatch(r"d(\d+)_n(\d+)\.txt", f)
        if m:
            up.append((int(m.group(1)) * int(m.group(2)), f"codes/{f}", int(m.group(1)), int(m.group(2))))
    codes = {}
    for lane, split in (("finish", FIN_SPLIT), ("fast", FAST_SPLIT)):
        for r in lane_runs(lane):
            if r["polish"] != "converged" or r["mtime"] < split or (lane == "finish" and r["secs"] >= 1700):
                continue
            k = (r["d"], r["n"])
            if k not in codes or r["checked"] < codes[k]["checked"]:
                r["lane"] = lane
                codes[k] = r
    conv = []
    for (d, n), r in codes.items():
        if output(f"runs/{r['lane']}/{r['file']}") is not None:
            conv.append((d * n, f"observations/data/runs/{r['lane']}/{r['file']}", d, n))
    for name, rows in (("uploaded", up), ("converged", conv)):
        with open(os.path.join(HERE, "lists", name + ".txt"), "w") as f:
            for _, p, d, n in sorted(rows):
                f.write(f"{p} {d} {n}\n")
        print(f"lists/{name}.txt: {len(rows)} codes")


if __name__ == "__main__":
    main()
