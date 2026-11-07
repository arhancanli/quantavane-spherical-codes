# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""The Galois relabellings of the second code X' (the analogue of Theorem 3.3(e) for X'), found exactly.

  python3 -I certify/note/galois_xprime.py

Input: data_xprime/class_matrix.gp (xprime.gp). The three automorphisms of L over Q(sqrt145) permute the classes of X'
as they do for X* (s1 -> -s1 exchanges the classes 4, 11 and 7, 9; s2 -> -s2 exchanges 3, 10; the new values t0', t1'
lie in Q(sqrt145) and are fixed). For each, an exhaustive individualisation-refinement search (aut116.isomorphism,
integers only) finds a relabelling pi of the 116 points with C'[pi(a), pi(b)] = map(C'[a, b]) for all a, b; the equality
is checked entry by entry here and again, independently, by check_xprime.gp, which also checks that each class map is
the conjugation of the values. Writes data_xprime/galois.gp (GALPERM, the format of data/aut_generators.gp) and
out/galois_xprime.json.
"""
import hashlib, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [p for p in os.environ.get("QV_LIBS", "").split(os.pathsep) if p]   # optional library folders
sys.path.insert(0, HERE)
import numpy as np
from aut116 import isomorphism

MAPS = (("s1->-s1", {4: 11, 11: 4, 7: 9, 9: 7}), ("s2->-s2", {3: 10, 10: 3}),
        ("s1->-s1,s2->-s2", {4: 11, 11: 4, 7: 9, 9: 7, 3: 10, 10: 3}))


def main():
    path = os.path.join(HERE, "data_xprime", "class_matrix.gp")
    C = np.array([json.loads(line) for line in open(path) if line.strip()], dtype=np.int64)
    assert C.shape == (116, 116) and (C == C.T).all() and (np.diag(C) == -1).all()
    gal = []
    for label, sw in MAPS:
        C2 = C.copy()
        for a, b in sw.items():
            C2[C == a] = b
        g = isomorphism(C2, C)               # C[g(a), g(b)] = C2[a, b] = map(C[a, b])
        ok = g is not None and bool((C[np.ix_(g, g)] == C2).all())
        gal.append({"conjugation": label, "class_map": {str(k): v for k, v in sw.items()},
                    "relabelling": None if g is None else [int(x) for x in g], "verified": ok})
    if not all(x["verified"] for x in gal):
        raise SystemExit("a conjugation of X' is not realised by a relabelling")
    with open(os.path.join(HERE, "data_xprime", "galois.gp"), "w") as f:
        f.write("\\\\ The Galois relabellings of X' (galois_xprime.py), checked by check_xprime.gp.\n")
        f.write("\\\\ GALPERM[j] = [class map as [k, image], relabelling pi with C[pi(a), pi(b)] = map(C[a, b])].\n")
        f.write("GALPERM = [" + ", ".join("[[" + "; ".join("%s, %d" % (k, v) for k, v in x["class_map"].items()) + "], ["
                                         + ",".join(map(str, x["relabelling"])) + "]]" for x in gal) + "];\n")
    out = {"program": "galois_xprime.py",
           "code_sha256": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
           "class_matrix_sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
           "relabellings": gal, "status": "3 of 3 conjugations of X' realised by relabellings"}
    json.dump(out, open(os.path.join(HERE, "out", "galois_xprime.json"), "w"), indent=1)
    print(out["status"])


if __name__ == "__main__":
    main()
