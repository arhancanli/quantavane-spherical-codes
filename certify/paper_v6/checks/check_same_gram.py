"""Closes the identification step of the certified distances for the exact-only confirmations (MCOM version).

For each of the 18 entries of certify/v5/distances.py, the configuration confirmed by certify/census/exact_strict.py
(recognised exact coordinates y in Q(sqrt D), or a recognised rational Gram matrix) and the configuration whose distance
to the stored code certify/v5/core.py bounds (the core Gram matrix recognised in Q(sqrt disc) and verified exactly)
are compared EXACTLY: every normalised inner product <y_a, y_b> / (|y_a| |y_b|) of the first must equal the entry of
the second (checked through signs and squares in Q(sqrt D)), and the two must have the same core. Equal Gram matrices
mean the two configurations differ by an orthogonal map, so the certified distance applies to the confirmed one.

  nice -n 10 python3 -I check_same_gram.py SNAPDIR OUT.json

The recognitions are re-run with the programs' own (untrusted) functions; the comparison is exact rational arithmetic.
"""
import hashlib, json, os, sys, time
from fractions import Fraction

CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
sys.path.insert(0, os.path.join(CERT, "v5"))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
import numpy as np
import core            # certify/v5/core.py
import exact_strict    # certify/census/exact_strict.py
import snapshot        # certify/census/snapshot.py

ENTRIES = [(3, 12, 5), (4, 120, 5), (6, 27, 3), (6, 72, 3), (4, 10, 1), (14, 42, 1), (14, 64, 1), (14, 70, 1),
           (20, 56, 1), (20, 81, 1), (20, 105, 1), (20, 120, 1), (21, 50, 1), (21, 77, 1), (21, 112, 1), (22, 100, 1),
           (21, 162, 1), (22, 275, 1)]


def fr(x):
    """A FLINT rational (flint.fmpq) or a Fraction as a Fraction."""
    return Fraction(int(x.p), int(x.q)) if hasattr(x, "p") else Fraction(x)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    snap, out = sys.argv[1], sys.argv[2]
    Qc = exact_strict.Q
    res = {"program": "check_same_gram.py",
           "code_sha256": {"check_same_gram.py": sha(os.path.abspath(__file__)),
                           "core.py": sha(os.path.join(CERT, "v5", "core.py")),
                           "exact_strict.py": sha(os.path.join(CERT, "census", "exact_strict.py"))},
           "entries": []}
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    for d, n, disc in ENTRIES:
        t0 = time.process_time()
        rec = {"d": d, "n": n, "disc": disc}
        C = core.build(snap, d, n, disc)
        rec["core_status"] = C["status"]
        V = C["V"]
        G5 = C["G"]
        pos = {a: i for i, a in enumerate(V)}
        pts = snapshot.read_points(snap, idx[(d, n)])
        dec = [[t.strip() for t in p] for p in pts]
        y, D = exact_strict.recognise_coords(dec)
        mism, checked, outside_bad = 0, 0, 0
        c5 = C["c"]
        if y is None:
            X = np.array([[float(t) for t in row] for row in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            G4, dev = exact_strict.recognise_gram(X)
            rec["v4_object"] = "rational Gram matrix" if G4 is not None else "not recognised"
            if G4 is None:
                rec["status"] = "FAIL: exact_strict recognition failed"
                res["entries"].append(rec)
                continue
            assert fr(c5.b) == 0, "core cosine irrational while the v4 Gram matrix is rational"
            for a in range(n):
                for b in range(a + 1, n):
                    if a in pos and b in pos:
                        g = G5[pos[a]][pos[b]]
                        checked += 1
                        if not (fr(g.b) == 0 and fr(g.a) == G4[a][b]):
                            mism += 1
                    else:
                        if not G4[a][b] < fr(c5.a):
                            outside_bad += 1
        else:
            rec["v4_object"] = f"coordinates in Q(sqrt {D})"
            Dm = max(D, disc)
            if D != 1 and disc != 1 and D != disc:
                rec["status"] = f"FAIL: fields differ (exact_strict D = {D}, core disc = {disc})"
                res["entries"].append(rec)
                continue
            Qc.D = Dm
            yy = [[Qc(x.a, x.b) for x in row] for row in y]
            N = [exact_strict.dot(r, r) for r in yy]
            cq = Qc(fr(c5.a), fr(c5.b))
            one = Qc(1)
            for a in range(n):
                for b in range(a + 1, n):
                    ip = exact_strict.dot(yy[a], yy[b])
                    nn = N[a] * N[b]
                    if a in pos and b in pos:
                        g = G5[pos[a]][pos[b]]
                        gq = Qc(fr(g.a), fr(g.b))
                        checked += 1
                        # ip / sqrt(nn) == g  <=>  same sign and ip^2 == g^2 nn
                        if ip.sign() != gq.sign() or not (ip * ip - gq * gq * nn).iszero():
                            mism += 1
                    else:
                        if exact_strict.cmp_cos(ip, nn, cq, one) >= 0:
                            outside_bad += 1
        rec.update({"core_size": len(V), "pairs_compared": checked, "mismatches": mism,
                    "pairs_outside_core_not_below_c": outside_bad,
                    "cpu_seconds": round(time.process_time() - t0, 2)})
        rec["status"] = "SAME GRAM MATRIX" if (mism == 0 and outside_bad == 0 and C["status"] == "ok") else "FAIL"
        print(d, n, rec["status"], rec["v4_object"], checked, mism, outside_bad, rec["cpu_seconds"], flush=True)
        res["entries"].append(rec)
    json.dump(res, open(out, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
