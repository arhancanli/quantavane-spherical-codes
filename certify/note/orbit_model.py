# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""Each orbit of Aut(X*) on pairs of points, described in the model of Theorem 3.3 (for Tables 1 and 4 of the note).

  python3 -I certify/note/orbit_model.py

Inputs: out/aut116.json (the 30 orbits on pairs, with all their pairs) and out/construct.json (the bijection from the
points of X* to the model: pentads 'A:...', pairs 'B:z:i,j' of a point z of PG(1,9) and two pencils at z, and P+-,
Q+-), both written by the exact programs. The circle geometry is rebuilt with the helpers of construct.py. For every
pair of every orbit the program computes the invariants below and checks that they are constant on the orbit, so each
orbit is described by them (exact, integers only):
  A x A  |O n O'|, and whether O, O' lie in the same A6-orbit;
  A x B  whether the point z of beta lies on O, and how many pencils of beta are of type I for O;
  B x B  same point or not; at the same point, whether the two pairs of pencils share a pencil; at different points,
         how many pencils the pair of beta, moved to z', shares with the pair of beta';
  P, Q   the signs (A_eta, B_eta).
Writes out/orbit_model.json.
"""
import itertools, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [p for p in os.environ.get("QV_LIBS", "").split(os.pathsep) if p]   # optional library folders
sys.path.insert(0, HERE)
import construct as M


def main():
    aut = json.load(open(os.path.join(HERE, "out", "aut116.json")))
    con = json.load(open(os.path.join(HERE, "out", "construct.json")))
    names = con["relabelling_Xstar_point_to_model_point"]
    signs = con["signs_of_Xstar_points"]
    GEN = [M.semilinear(M.ONE, M.ONE, M.ZERO, M.ONE, 0), M.semilinear(M.W, M.ZERO, M.ZERO, M.ONE, 0),
           M.semilinear(M.ZERO, M.ONE, M.ONE, M.ZERO, 0), M.semilinear(M.ONE, M.ZERO, M.ZERO, M.ONE, 1)]
    GAMMA = M.closure(GEN)
    base_circle = frozenset(M.IDX[x] for x in (M.INF, (0, 0), (1, 0), (2, 0)))
    circles = sorted({frozenset(g[i] for i in base_circle) for g in GAMMA}, key=sorted)
    pencils = {}
    for p in range(10):               # as in construct.py (same order, so the pencil indices of the names agree)
        through = [c for c in circles if p in c]
        groups = []
        for c in through:
            for grp in groups:
                if all(len(c & c2) == 1 for c2 in grp):
                    grp.append(c)
                    break
            else:
                groups.append([c])
        pencils[p] = [frozenset(gp) for gp in groups]

    def parse(nm):
        if nm.startswith("A:"):
            return ("A", frozenset(int(x) for x in nm[2:].split(",")))
        if nm.startswith("B:"):
            _, z, pp = nm.split(":")
            z = int(z)
            return ("B", (z, frozenset(pencils[z][int(i)] for i in pp.split(","))))
        return (nm[0], 1 if nm[1] == "+" else -1)

    def type_I(O, pen):
        return sum(1 for c in pen if len(c & O) == 3) == 1

    def transport(z, z2, pen):
        c = next(c for c in pen if z2 in c)
        return next(x for x in pencils[z2] if c in x)

    def invariant(a, b):
        (ka, va), (kb, vb) = parse(names[a]), parse(names[b])
        sa, sb = signs[a], signs[b]
        if (ka, kb) in (("B", "A"), ("P", "A"), ("Q", "A"), ("P", "B"), ("Q", "B"), ("Q", "P")):
            (ka, va, sa), (kb, vb, sb) = (kb, vb, sb), (ka, va, sa)
        if (ka, kb) == ("A", "A"):
            return ("AA", len(va & vb), "same" if sa == sb else "different")
        if (ka, kb) == ("A", "B"):
            z, pp = vb
            return ("AB", "z in O" if z in va else "z not in O", sum(type_I(va, x) for x in pp))
        if (ka, kb) == ("B", "B"):
            (z, e), (z2, f) = va, vb
            if z == z2:
                return ("BB", "same point", len(e & f))
            return ("BB", "different points", len(frozenset(transport(z, z2, x) for x in e) & f))
        if (ka, kb) in (("A", "P"), ("A", "Q"), ("B", "P"), ("B", "Q")):
            return (ka + kb, "same sign" if sa == sb else "opposite signs")
        return (ka + kb, "same sign" if sa == sb else "opposite signs")

    out = []
    for o in aut["orbitals"]:
        invs = {invariant(a, b) for a, b, _ in o["pairs"]}
        assert len(invs) == 1, (o["class"], o["orbits"], invs)
        inv = invs.pop()
        out.append({"class": o["class"], "orbits": o["orbits"], "size": o["size"], "split": o["split"],
                    "rep": o["rep"], "invariant": list(inv)})
    # the invariants separate the 30 orbits
    keys = [(x["class"], tuple(x["invariant"])) for x in out]
    rep = {"program": "orbit_model.py", "orbits": out, "invariants_distinct": len(set(keys)) == len(keys),
           "invariant_alone_distinct": len({tuple(x["invariant"]) for x in out}) == len(out)}
    # a worked example for the figure of the note: beta0 = (inf, {directions 1, w}) and a pentad O through inf with
    # both pencils of beta0 of type I (so O is in A_+ and the class of (O, beta0) is 3)
    INFI = M.IDX[M.INF]
    def pencil_dir(dv):
        return next(pn for pn in pencils[INFI] if frozenset([INFI] + [M.IDX[M.add(M.ZERO, M.mul((t, 0), dv))]
                                                                   for t in range(3)]) in pn)
    beta0 = frozenset([pencil_dir(M.ONE), pencil_dir(M.W)])
    pent = sorted({frozenset(parse(nm)[1]) for nm in names if nm.startswith("A:")}, key=sorted)
    Os = [O for O in pent if INFI in O and all(type_I(O, pn) for pn in beta0)]
    O = Os[0]
    elt = lambda i: M.PI[i]
    lines = {}
    for dname, dv in (("1", M.ONE), ("w", M.W)):
        pn = pencil_dir(dv)
        lines[dname] = [sorted(elt(i) for i in c if i != INFI) for c in pn if len(c & O) == 3]
    i0 = names.index("B:%d:%s" % (INFI, ",".join(str(pencils[INFI].index(pn))
                                                     for pn in sorted(beta0, key=pencils[INFI].index))))
    iO = names.index("A:" + ",".join(map(str, sorted(O))))
    rep["example"] = {"pentad_finite_points_(a,b)=a+bw": sorted(elt(i) for i in O if i != INFI),
                      "lines_meeting_O_in_three_points": lines, "beta0_point_index": i0, "pentad_point_index": iO,
                      "sign_of_pentad": signs[iO], "sign_of_beta0": signs[i0],
                      "class_of_the_pair": invariant(iO, i0)}
    json.dump(rep, open(os.path.join(HERE, "out", "orbit_model.json"), "w"), indent=1)
    for x in out:
        print(x["class"], x["orbits"], x["size"], "split" if x["split"] else "", x["invariant"])
    print({k: v for k, v in rep.items() if k != "orbits"})


if __name__ == "__main__":
    main()
