"""Write the stored decimal coordinates of one entry of the archived table to a file, one point per line.
  python3 -I certify/v5/extract_code.py SNAPDIR D N OUT.txt"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "census"))
import snapshot
snapdir, d, n, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
idx = {(r["d"], r["n"]): r for r in snapshot.load_index() if r["listed"]}
pts = snapshot.read_points(snapdir, idx[(d, n)])
open(out, "w").write("".join(",".join(p) + "\n" for p in pts))
