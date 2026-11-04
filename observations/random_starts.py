"""The 25 Oct from-scratch tests (data/runs/tests/*.log: d N best table improved secs): how far behind the table.

  python3 -I random_starts.py > random_starts.out
"""
import os
from decimal import Decimal
T = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "runs", "tests")
for f in sorted(os.listdir(T)):
    if f.endswith(".log"):
        p = open(os.path.join(T, f)).read().split()
        print(f"{f:22s} d{p[0]} N{p[1]}: ours {p[2]} table {p[3]} -> behind by {float(Decimal(p[2]) - Decimal(p[3])):.2e} ({p[5]} s)")
