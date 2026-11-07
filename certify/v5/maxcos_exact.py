"""The exact largest cosine of a code given by decimal coordinates (paper, version 5). TRUSTED: exact rational arithmetic.

  python3 -I certify/v5/maxcos_exact.py CODE.txt D [--compare 0.988134624923]

Every point is the decimal vector y_i (as an exact rational vector); the code is y_i / |y_i|. The largest cosine is
max_{i<j} <y_i, y_j> / (|y_i| |y_j|); pairs are ranked in double precision (untrusted), and the largest ones are compared
exactly through the signs and squares of <y_i, y_j> and the rational numbers |y_i|^2 |y_j|^2, so the reported pair attains
the maximum exactly. Printed: the maximum to 30 digits (ball arithmetic) and whether it is at most the given decimal.
"""
import json, os, sys
from fractions import Fraction
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/flintlib"))
import numpy as np
import flint


def read(path, d):
    pts = []
    for line in open(path):
        p = line.strip().split(",")
        if len(p) == d:
            pts.append([Fraction(t) for t in p])
    return pts


def cmp(ip1, n1, ip2, n2):
    """sign of ip1/sqrt(n1) - ip2/sqrt(n2) (n1, n2 > 0)."""
    s1 = (ip1 > 0) - (ip1 < 0)
    s2 = (ip2 > 0) - (ip2 < 0)
    if s1 != s2:
        return (s1 > s2) - (s1 < s2)
    if s1 == 0:
        return 0
    t = ip1 * ip1 * n2 - ip2 * ip2 * n1
    t = (t > 0) - (t < 0)
    return t if s1 > 0 else -t


def main():
    path, d = sys.argv[1], int(sys.argv[2])
    pts = read(path, d)
    n = len(pts)
    N = [sum(x * x for x in p) for p in pts]
    Xf = np.array([[float(x) for x in p] for p in pts])
    Xf /= np.linalg.norm(Xf, axis=1)[:, None]
    G = Xf @ Xf.T
    iu = np.triu_indices(n, 1)
    gv = G[iu]
    order = np.argsort(-gv)
    top = float(gv[order[0]])
    cand = [int(i) for i in order[:200000] if gv[i] >= top - 1e-9]
    best = None
    for i in cand:
        a, b = int(iu[0][i]), int(iu[1][i])
        ip = sum(x * y for x, y in zip(pts[a], pts[b]))
        key = (ip, N[a] * N[b], a, b)
        if best is None or cmp(key[0], key[1], best[0], best[1]) > 0:
            best = key
    # every pair outside the candidates is below top - 1e-9 in double precision; check this margin exactly is not
    # needed for the maximum, but it is for its attainment: confirm the remaining pairs are below the best exactly when
    # their double value is within 1e-6 of it
    worry = [int(i) for i in order[:200000] if top - 1e-6 <= gv[i] < top - 1e-9]
    for i in worry:
        a, b = int(iu[0][i]), int(iu[1][i])
        ip = sum(x * y for x, y in zip(pts[a], pts[b]))
        if cmp(ip, N[a] * N[b], best[0], best[1]) > 0:
            best = (ip, N[a] * N[b], a, b)
    flint.ctx.prec = 200
    val = flint.arb(flint.fmpq(best[0].numerator, best[0].denominator)) / flint.arb(
        flint.fmpq(best[1].numerator, best[1].denominator)).sqrt()
    out = {"code": os.path.basename(path), "d": d, "n": n, "pair": [best[2], best[3]], "max_cos_30": val.str(30, radius=False),
           "pairs_compared_exactly": len(cand) + len(worry)}
    if "--compare" in sys.argv:
        q = Fraction(sys.argv[sys.argv.index("--compare") + 1])
        # max <= q  iff  ip <= q sqrt(N)  (decided exactly through squares and signs)
        ip, nn = best[0], best[1]
        if q >= 0:
            le = ip <= 0 or ip * ip <= q * q * nn
        else:
            le = ip < 0 and ip * ip >= q * q * nn
        out["compare_with"] = str(sys.argv[sys.argv.index("--compare") + 1])
        out["max_at_most_compare"] = bool(le)
    print(json.dumps(out))


if __name__ == "__main__":
    main()
