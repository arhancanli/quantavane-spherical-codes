"""Independent check of a spherical code (written apart from engine/spherical/sphopt.py; standard library only).
Reads one point per line, coordinates separated by commas (other lines are ignored, as on spherical-codes.org),
normalises every point in 50-digit decimal arithmetic and reports the largest cosine between two distinct points.

  python3 -I independent/sph_check.py CODE.txt D N [TABLE_COS]
Prints the point count, the dimension, the largest |norm - 1| before normalising, the largest cosine to 20 digits and,
given TABLE_COS, whether it is smaller than the table's value rounded to 12 decimals (the table's precision)."""
import sys
from decimal import Decimal, getcontext

getcontext().prec = 50


def main():
    path, d, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    pts = []
    for line in open(path):
        parts = line.strip().split(",")
        try:
            vals = [Decimal(x.strip()) for x in parts]
        except Exception:
            continue
        if len(vals) == d:
            pts.append(vals)
    assert len(pts) == n, f"{len(pts)} points of dimension {d}, expected {n}"
    worst_norm = Decimal(0)
    unit = []
    for p in pts:
        norm = sum(x * x for x in p).sqrt()
        worst_norm = max(worst_norm, abs(norm - 1))
        unit.append([x / norm for x in p])
    best = Decimal(-2)
    for i in range(n):
        a = unit[i]
        for j in range(i + 1, n):
            c = sum(x * y for x, y in zip(a, unit[j]))
            if c > best:
                best = c
    print(f"points {n} dimension {d} max|norm-1| {worst_norm:.3e} max_cos {best:.20f}")
    if len(sys.argv) > 4:
        table = Decimal(sys.argv[4])
        better = best.quantize(Decimal("1e-12")) < table
        print(f"table {table} ours rounded {best.quantize(Decimal('1e-12'))} improves at 12 decimals: {'YES' if better else 'NO'}")


if __name__ == "__main__":
    main()
