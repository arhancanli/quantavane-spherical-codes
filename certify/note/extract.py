# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""Write the stored coordinates of the entry (22,116) of the archived table to data/stored_d22_n116.txt.

  python3 -I certify/note/extract.py SNAPDIR

SNAPDIR holds data.txt of the archived table (DSpace@MIT handle 1721.1/153543, data as of 31 July 2026). The entry is
found by its header lines ("Dimension: 22", "Number of points: 116"); the 116 lines after "Points:" are copied verbatim.
Standard library only; nothing here is a check. The SHA-256 of data.txt and of the output are printed.
"""
import hashlib, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    snapdir = sys.argv[1]
    path = os.path.join(snapdir, "data.txt")
    d, n = 22, 116
    pts, state = [], 0
    with open(path, encoding="utf-8") as f:
        prev = ""
        for line in f:
            s = line.strip()
            if state == 0:
                if s == f"Number of points: {n}" and prev == f"Dimension: {d}":
                    state = 1
                prev = s
            elif state == 1:
                if s == "Points:":
                    state = 2
            else:
                parts = s.split(",")
                if len(parts) != d:
                    break
                pts.append(s)
                if len(pts) == n:
                    break
    if len(pts) != n:
        raise SystemExit(f"read {len(pts)} points")
    out = os.path.join(HERE, "data", "stored_d22_n116.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out, "w").write("\n".join(pts) + "\n")
    sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
    print("data.txt sha256", sha(path))
    print("stored_d22_n116.txt sha256", sha(out))


if __name__ == "__main__":
    main()
