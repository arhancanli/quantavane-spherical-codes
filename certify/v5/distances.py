"""Certified distances for the exact-arithmetic-only confirmations of version 4 (paper, version 5, Section 4).

  python3 -I certify/v5/distances.py SNAPDIR OUT.json D,N,DISC [D,N,DISC ...]

exact_strict.py (version 4) proved, with no floating-point arithmetic, that the exact configuration of each of these
entries is a strict local minimum modulo rotations and rattlers; the exact configuration was defined by recognised
coordinates (in Q(sqrt D)) or by a recognised rational Gram matrix, and its distance to the stored code was not
bounded. Here core.py recognises the inner products of the core again, verifies the Gram matrix exactly (rank d,
contacts, rattlers below c) and bounds, in ball arithmetic with an exactly orthogonal rational rotation, the distance
between the configuration with this Gram matrix and the exactly normalised stored code. Two configurations with the
same Gram matrix differ by an orthogonal map, so the bound applies to the configuration of exact_strict.py whenever the
two recognised matrices coincide, which the program checks against exact_strict.py's recorded cosine and counts.
"""
import hashlib, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import core

SHA = {}
for _f in ("distances.py", "core.py", "exact2.py", "kfield.py", "dist.py"):
    SHA[_f] = hashlib.sha256(open(os.path.join(HERE, _f), "rb").read()).hexdigest()


def main():
    snapdir, out = sys.argv[1], sys.argv[2]
    ex4 = os.path.join(os.path.dirname(HERE), "census", "exact_v4", "strict")
    res = {"program": "distances.py", "code_sha256": SHA, "entries": []}
    for spec in sys.argv[3:]:
        d, n, disc = map(int, spec.split(","))
        t0 = time.time()
        C = core.build(snapdir, d, n, disc)
        rec = {"d": d, "n": n, "disc": disc, "status": C["status"]}
        rec.update(C.get("log", {}))
        try:
            v4 = json.load(open(os.path.join(ex4, f"d{d}_n{n}.json")))
            rec["v4_status"] = v4.get("status")
            rec["v4_contacts"] = v4.get("contacts")
            rec["v4_rattlers"] = v4.get("rattlers")
            rec["same_contacts_as_v4"] = (v4.get("contacts") == rec.get("contacts")
                                          and v4.get("rattlers") == rec.get("rattlers"))
        except FileNotFoundError:
            rec["v4_status"] = None
        rec["seconds"] = round(time.time() - t0, 1)
        res["entries"].append(rec)
        print(json.dumps(rec, default=str), flush=True)
    json.dump(res, open(out, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
