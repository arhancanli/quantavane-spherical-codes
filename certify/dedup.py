"""Which uploaded codes are, number for number, the same code as a converged code (the same file content apart from
comment lines). Writes lists/uploaded_same.txt ("uploaded path  converged path") and lists/uploaded_rest.txt (the
uploaded codes that are not), so that each distinct code is certified once; analysis.py uses the converged code's
certificate for the first kind.

  python3 -I certify/dedup.py
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def body(p):
    return [l.strip() for l in open(os.path.join(ROOT, p)) if l.strip() and not l.startswith("#")]


def read(name):
    return [l.split() for l in open(os.path.join(HERE, "lists", name + ".txt")) if l.strip()]


conv = {(int(d), int(n)): p for p, d, n in read("converged")}
same, rest = [], []
for p, d, n in read("uploaded"):
    q = conv.get((int(d), int(n)))
    if q is not None and body(p) == body(q):
        same.append(f"{p} {q}\n")
    else:
        rest.append(f"{p} {d} {n}\n")
open(os.path.join(HERE, "lists", "uploaded_same.txt"), "w").writelines(same)
open(os.path.join(HERE, "lists", "uploaded_rest.txt"), "w").writelines(rest)
print(f"uploaded codes identical to a converged code: {len(same)}; others: {len(rest)}")
