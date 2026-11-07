# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""The certificate data for (22,116) that do not depend on the symmetry computation (UNTRUSTED: every fact used is
re-checked by check.gp).

  python3 -I certify/note/classes.py RECOGNITION.json

Inputs: data/stored_d22_n116.txt (extract.py) and the recognition file of version 5
(certify/v5/out/recognition_d22_n116.json), from which only the thirteen exact values (integer relations found by PSLQ),
the field and the basis B are taken. Writes, in PARI/GP syntax:
  data/values.gp        the field L = Q(r, s1, s2), r^2 = 145, s1^2 = D1, s2^2 = D2, and the values t_0..t_12 as
                        coefficient vectors in the basis 1, r, s1, r s1, s2, r s2, s1 s2, r s1 s2;
  data/class_matrix.gp  116 lines, line a = [k(a,0), ..., k(a,115)], the class of each pair (-1 on the diagonal);
  data/basis.gp         the 22 points of B (0-based).
The class of a pair is the index of the value nearest to the stored inner product (stored coordinates read exactly as
decimals, normalised and multiplied at 60 digits). Printed: the largest distance of a stored inner product from its
value, the smallest distance between two values, and whether the classes agree with those of the recognition file.
"""
import json, os, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [p for p in os.environ.get("QV_LIBS", "").split(os.pathsep) if p]   # optional library folders
import mpmath
from mpmath import mp, mpf


def main():
    rec = json.load(open(sys.argv[1]))
    n, d = 116, 22
    # ---------------------------------------------------------------- the field and the values
    def k0(rel):                     # a Delta + b + c r = 0
        a, b, c = rel
        return (Fraction(-b, a), Fraction(-c, a))
    D1, D2 = k0(rec["Delta1_rel"]), k0(rec["Delta2_rel"])
    assert rec["sigma1"] == -1 and rec["sigma2"] == -1
    vals = []
    for rel in rec["class_rel"]:     # r0 t + sum c_i basis_i = 0
        r0 = rel[0]
        vals.append([Fraction(-x, r0) for x in rel[1:]])
    mp.dps = 60
    r = mpmath.sqrt(145)
    s1 = -mpmath.sqrt(D1[0].numerator / mpf(D1[0].denominator) + D1[1].numerator / mpf(D1[1].denominator) * r)
    s2 = -mpmath.sqrt(D2[0].numerator / mpf(D2[0].denominator) + D2[1].numerator / mpf(D2[1].denominator) * r)
    basis = [mpf(1), r, s1, r * s1, s2, r * s2, s1 * s2, r * s1 * s2]
    tv = [mpmath.fsum(mpf(c.numerator) / c.denominator * b for c, b in zip(v, basis)) for v in vals]
    order = sorted(range(len(tv)), key=lambda k: tv[k])
    assert order == list(range(len(tv))), "values not increasing"
    gap = min(tv[k + 1] - tv[k] for k in range(len(tv) - 1))
    # ---------------------------------------------------------------- classes from the stored code
    rows = [line.split(",") for line in open(os.path.join(HERE, "data", "stored_d22_n116.txt")).read().split()]
    assert len(rows) == n and all(len(x) == d for x in rows)
    Xs = []
    for row in rows:
        v = [mpf(t) for t in row]
        nr = mpmath.sqrt(mpmath.fsum(x * x for x in v))
        Xs.append([x / nr for x in v])
    cls = [[-1] * n for _ in range(n)]
    dev = mpf(0)
    for a in range(n):
        for b in range(a + 1, n):
            ip = mpmath.fsum(x * y for x, y in zip(Xs[a], Xs[b]))
            k = min(range(len(tv)), key=lambda j: abs(ip - tv[j]))
            dev = max(dev, abs(ip - tv[k]))
            cls[a][b] = cls[b][a] = k
    pc = rec["pair_class"]
    same = [cls[a][b] for a in range(n) for b in range(a + 1, n)] == pc
    # ---------------------------------------------------------------- write
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    fr = lambda q: str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
    with open(os.path.join(HERE, "data", "values.gp"), "w") as f:
        f.write("\\\\ The field L = Q(r, s1, s2): r^2 = 145, s1^2 = D1, s2^2 = D2 with D = D[1] + D[2]*r; real embedding\n"
                "\\\\ r > 0, s1 < 0, s2 < 0. Row k+1 of the matrix TC is the value t_k as a coefficient vector in the basis\n"
                "\\\\ 1, r, s1, r*s1, s2, r*s2, s1*s2, r*s1*s2. Written by classes.py from the version 5 recognition.\n")
        f.write(f"D1 = [{fr(D1[0])}, {fr(D1[1])}];\n")
        f.write(f"D2 = [{fr(D2[0])}, {fr(D2[1])}];\n")
        f.write("TC = [" + "; ".join(", ".join(fr(c) for c in v) for v in vals) + "];\n")
    with open(os.path.join(HERE, "data", "class_matrix.gp"), "w") as f:
        for a in range(n):
            f.write("[" + ",".join(str(x) for x in cls[a]) + "]\n")
    with open(os.path.join(HERE, "data", "basis.gp"), "w") as f:
        f.write("[" + ",".join(str(x) for x in rec["basis"]) + "]\n")
    print(json.dumps({"largest_distance_stored_ip_to_value": mpmath.nstr(dev, 3),
                      "smallest_gap_between_values": mpmath.nstr(gap, 6),
                      "classes_agree_with_recognition": same,
                      "pairs_per_class": [sum(1 for a in range(n) for b in range(a + 1, n) if cls[a][b] == k)
                                          for k in range(len(tv))]}))


if __name__ == "__main__":
    main()
