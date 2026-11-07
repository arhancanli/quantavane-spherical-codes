# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""The symmetry group of the code X* of (22,116): the permutations of the 116 points that preserve the class matrix.

  python3 -I certify/note/aut116.py

Input: data/class_matrix.gp (classes.py). Since the thirteen values are distinct and X* spans R^22, these permutations
are exactly the orthogonal maps of R^22 that permute X* (each permutation preserves the Gram matrix, hence extends to a
unique orthogonal map, and conversely).

Method (exact; integers only): individualisation and refinement on the complete graph whose edges are coloured by the
classes. The refinement splits cells by the number of neighbours of each colour in each cell, and orders the new cells
by these counts, so it commutes with relabelling the points; the search that maps one individualised point to another
is exhaustive. The order is the product of the orbit lengths along a stabiliser chain (base points chosen on the first
path to a discrete partition), and every generator found is checked to preserve the class matrix entry by entry. The
whole computation is repeated after a random relabelling of the points as a cross-check of the order.
Finally, for each conjugation of L over Q(sqrt145) (which permutes the classes), a relabelling of the points
carrying the class matrix to its conjugate is searched for in the same way.

Then, by enumerating the group: element orders; the derived subgroup G' (order, perfect, simple: the normal closure in
G' of every non-identity element is G'); the centraliser of G' in G (trivial, so G embeds in Aut(G')); the three
subgroups of index 2 containing G'; the orbits on points, on unordered pairs (with their classes) and on contacts; and
for the index-2 subgroup M fixing the four points of the two smallest orbits, which orbits split under M (with the
sign +1/-1 of the two halves). Writes out/aut116.json, data/aut_generators.gp and data/orbitals.gp.
"""
import collections, hashlib, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [p for p in os.environ.get("QV_LIBS", "").split(os.pathsep) if p]   # optional library folders
import numpy as np


def load_classes():
    rows = [json.loads(line) for line in open(os.path.join(HERE, "data", "class_matrix.gp")) if line.strip()]
    C = np.array(rows, dtype=np.int64)
    assert C.shape[0] == C.shape[1] and (C == C.T).all() and (np.diag(C) == -1).all()
    return C


class Refiner:
    def __init__(self, C):
        self.C = C
        self.n = C.shape[0]
        k = int(C.max()) + 1
        D = C.copy()
        D[np.diag_indices(self.n)] = k          # the diagonal is its own colour
        self.OH = np.stack([(D == j).astype(np.int64) for j in range(k + 1)])

    def refine(self, cells):
        """Equitable refinement; returns the cells and a trace (the sorted signatures), both relabelling-invariant."""
        n, trace = self.n, []
        while True:
            P = np.zeros((n, len(cells)), dtype=np.int64)
            for j, cell in enumerate(cells):
                P[cell, j] = 1
            sig = np.einsum("kab,bj->akj", self.OH, P).reshape(n, -1)
            new, changed = [], False
            for cell in cells:
                if len(cell) == 1:
                    new.append(cell)
                    continue
                groups = {}
                for v in cell:
                    groups.setdefault(sig[v].tobytes(), []).append(v)
                keys = sorted(groups)
                trace.append(tuple((key, len(groups[key])) for key in keys))
                changed |= len(keys) > 1
                new.extend(groups[key] for key in keys)
            cells = new
            if not changed:
                return cells, trace


def indiv(cells, i, u):
    return cells[:i] + [[u], [x for x in cells[i] if x != u]] + cells[i + 1:]


def first_open(cells):
    return next((i for i, c in enumerate(cells) if len(c) > 1), None)


def is_aut(C, g):
    return bool((C[np.ix_(g, g)] == C).all())


def extend(R, P1, P2):
    """An automorphism mapping the cells of P1 onto those of P2 in order, or None (exhaustive)."""
    i = first_open(P1)
    if i is None:
        g = np.empty(R.n, dtype=np.int64)
        for c1, c2 in zip(P1, P2):
            g[c1[0]] = c2[0]
        return g if is_aut(R.C, g) else None
    u = P1[i][0]
    Q1, t1 = R.refine(indiv(P1, i, u))
    for v in P2[i]:
        Q2, t2 = R.refine(indiv(P2, i, v))
        if t1 != t2:
            continue
        g = extend(R, Q1, Q2)
        if g is not None:
            return g
    return None


def isomorphism(C1, C2):
    """A permutation g with C2[g(a), g(b)] = C1[a, b] for all a, b, or None (exhaustive, as in extend)."""
    R1, R2 = Refiner(C1), Refiner(C2)
    P1, t1 = R1.refine([list(range(R1.n))])
    P2, t2 = R2.refine([list(range(R2.n))])
    if t1 != t2:
        return None

    def ext(Q1, Q2):
        i = first_open(Q1)
        if i is None:
            g = np.empty(R1.n, dtype=np.int64)
            for c1, c2 in zip(Q1, Q2):
                g[c1[0]] = c2[0]
            return g if bool((C2[np.ix_(g, g)] == C1).all()) else None
        u = Q1[i][0]
        S1, s1 = R1.refine(indiv(Q1, i, u))
        for v in Q2[i]:
            S2, s2 = R2.refine(indiv(Q2, i, v))
            if s1 != s2:
                continue
            g = ext(S1, S2)
            if g is not None:
                return g
        return None
    return ext(P1, P2)


def orbit(b, gens):
    orb, todo = {b}, [b]
    while todo:
        x = todo.pop()
        for g in gens:
            y = int(g[x])
            if y not in orb:
                orb.add(y)
                todo.append(y)
    return orb


def automorphisms(C):
    R = Refiner(C)
    P0, _ = R.refine([list(range(R.n))])
    parts, base, P = [P0], [], P0
    while first_open(P) is not None:
        i = first_open(P)
        base.append(P[i][0])
        P, _ = R.refine(indiv(P, i, P[i][0]))
        parts.append(P)
    gens, lengths = [], []
    for lev in reversed(range(len(base))):
        Pl = parts[lev]
        i = first_open(Pl)
        b = base[lev]
        sg = [g for g in gens if all(g[x] == x for x in base[:lev])]
        orb = orbit(b, sg)
        for v in Pl[i]:
            if v in orb:
                continue
            Q, _ = R.refine(indiv(Pl, i, v))
            g = extend(R, parts[lev + 1], Q)
            if g is not None:
                assert is_aut(C, g)
                gens.append(g)
                sg.append(g)
                orb = orbit(b, sg)
        lengths.append(len(orb))
    lengths.reverse()
    order = 1
    for x in lengths:
        order *= x
    return order, base, lengths, gens, [len(c) for c in P0]


def main():
    t0 = time.time()
    C = load_classes()
    n = C.shape[0]
    order, base, lengths, gens, eq = automorphisms(C)
    # cross-check after a random relabelling
    rng = np.random.default_rng(1440)
    perm = rng.permutation(n)
    C2 = C[np.ix_(perm, perm)]
    order2, base2, lengths2, gens2, _ = automorphisms(C2)
    assert order2 == order, ("the order depends on the labelling", order, order2)
    assert all(is_aut(C2, g) for g in gens2)
    # ---------------------------------------------------------------- enumerate the group
    e = tuple(range(n))
    T = [tuple(int(x) for x in g) for g in gens]
    mul = lambda g, h: tuple(g[x] for x in h)            # (g h)(x) = g(h(x))

    def inv(g):
        r = [0] * n
        for x in range(n):
            r[g[x]] = x
        return tuple(r)

    def closure(S):
        H, todo = {e}, [e]
        while todo:
            x = todo.pop()
            for s in S:
                y = mul(s, x)
                if y not in H:
                    H.add(y)
                    todo.append(y)
        return H
    G = closure(T)
    assert len(G) == order

    def eorder(g):
        k, h = 1, g
        while h != e:
            h, k = mul(g, h), k + 1
        return k
    ords = {g: eorder(g) for g in G}

    def normal_closure(S, within_gens):
        """The subgroup generated by the conjugates of S under the group generated by within_gens."""
        H = closure(S)
        while True:
            new = {mul(mul(w, h), inv(w)) for h in H for w in within_gens} - H
            if not new:
                return H
            H = closure(list(H | new))
    comm = [mul(mul(inv(a), inv(b)), mul(a, b)) for a in T for b in T]
    D = normal_closure(comm, T)                           # G' = normal closure of the commutators of generators
    Dg = list(D)
    Dgens = []                                            # a small generating set of G'
    Hc = {e}
    for x in Dg:
        if x not in Hc:
            Dgens.append(x)
            Hc = closure(Dgens)
            if len(Hc) == len(D):
                break
    comm2 = [mul(mul(inv(a), inv(b)), mul(a, b)) for a in Dgens for b in Dgens]
    perfect = len(normal_closure(comm2, Dgens)) == len(D)
    classes_D = []
    seen = set()
    for x in Dg:
        if x in seen:
            continue
        cl = {mul(mul(w, x), inv(w)) for w in D}
        seen |= cl
        classes_D.append((eorder(x), len(cl), x))
    simple = all(len(normal_closure([x], Dgens)) == len(D) for (o, s, x) in classes_D if x != e)
    centraliser = [g for g in G if all(mul(g, x) == mul(x, g) for x in Dgens)]
    # the index-2 subgroups containing G': G/G' has order 4; a coset is determined by g modulo G'
    cosets = {}
    for g in G:
        if g in D:
            continue
        key = min(mul(g, x) for x in D)
        cosets.setdefault(key, set()).add(g)
    cos = list(cosets.values())
    sub2 = []
    for cset in cos:
        H2 = D | cset
        st = collections.Counter(ords[g] for g in H2)
        sub2.append({"order": len(H2), "element_orders": dict(sorted(st.items())), "set": H2})
    # ---------------------------------------------------------------- orbits
    def orbits_points(grp):
        seen, out = set(), []
        for v in range(n):
            if v not in seen:
                o = sorted({g[v] for g in grp})
                seen |= set(o)
                out.append(o)
        return out
    pts = orbits_points(G)
    small = sorted(v for o in pts if len(o) == 2 for v in o)
    M = [s for s in sub2 if all(g[v] == v for g in s["set"] for v in small)]
    assert len(M) == 1
    Mset = M[0]["set"]
    eps = lambda g: 1 if g in Mset else -1
    name = {}
    for o in pts:
        for v in o:
            name[v] = {72: "A", 40: "B"}.get(len(o))
    # P: the 2-point orbit whose points meet A in classes 9 and 12; Q: in classes 7 and 12
    for o in pts:
        if len(o) == 2:
            ks = {int(C[o[0], a]) for a in range(n) if name.get(a) == "A"}
            name[o[0]] = name[o[1]] = "P" if ks == {9, 12} else "Q" if ks == {7, 12} else "?"
    pts_M = orbits_points(Mset)
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    norm = lambda a, b: (a, b) if a < b else (b, a)
    seen, orbitals = set(), []
    for p in pairs:
        if p in seen:
            continue
        og = sorted({norm(g[p[0]], g[p[1]]) for g in G})
        seen |= set(og)
        om = sorted({norm(g[p[0]], g[p[1]]) for g in Mset})
        split = len(om) < len(og)
        sign = {q: (1 if q in set(om) else -1) for q in og} if split else None
        orbitals.append({"class": int(C[p]), "size": len(og), "orbits": "".join(sorted(name[p[0]] + name[p[1]])),
                         "split": split, "rep": list(p),
                         "pairs": [[q[0], q[1], (sign[q] if split else 0)] for q in og]})
    orbitals.sort(key=lambda o: (o["class"], o["orbits"], o["size"], o["rep"]))
    types = [o for o in orbitals if o["class"] == int(C.max())]   # the contacts: the largest class
    def name2(sg):
        eo = sg["element_orders"]
        return "M10" if (8 in eo and 6 not in eo and 10 not in eo) else "S6" if 6 in eo else "PGL(2,9)" if 10 in eo else "?"
    names2 = [name2(sg) for sg in sub2]
    assert sorted(names2) == ["M10", "PGL(2,9)", "S6"]
    char_idx = [names2.index(k) for k in ("M10", "S6", "PGL(2,9)")]
    chargen = [[1 if g in sub2[i]["set"] else -1 for i in char_idx] for g in T]
    # the conjugations of L over Q(sqrt145) permute the classes: s1 -> -s1 swaps 4<->11 and 7<->9, s2 -> -s2 swaps
    # 3<->10 (read off the values in data/values.gp); is the conjugate class matrix a relabelling of C?
    galois = []
    for label, sw in (("s1->-s1", {4: 11, 11: 4, 7: 9, 9: 7}), ("s2->-s2", {3: 10, 10: 3}),
                      ("s1->-s1,s2->-s2", {4: 11, 11: 4, 7: 9, 9: 7, 3: 10, 10: 3})):
        C2 = C.copy()
        for a_, b_ in sw.items():
            C2[C == a_] = b_
        g = isomorphism(C, C2)
        galois.append({"conjugation": label, "class_map": {str(k): v for k, v in sw.items()},
                       "relabelling": None if g is None else [int(x) for x in g],
                       "verified": g is not None and bool((C2[np.ix_(g, g)] == C).all())})
    # each point of A (B) determines the Sylow 5-subgroup (3-subgroup) of its stabiliser: the fibres of these maps
    fibres = {}
    for nm_, pr in (("A", 5), ("B", 3)):
        fib = collections.defaultdict(list)
        for a in [v for v in range(n) if name[v] == nm_]:
            st_ = [g for g in G if g[a] == a]
            sub = frozenset(g for g in st_ if ords[g] in ((1, 5) if pr == 5 else (1, 3, 9)))
            fib[sub].append(a)
        allsyl = len([g for g in G if ords[g] == pr]) // (4 if pr == 5 else 8)
        fibres[nm_] = {"sylow_order": sorted({len(k) for k in fib}), "images": len(fib), "sylow_subgroups_in_G": allsyl,
                       "fibre_sizes": sorted(collections.Counter(len(v) for v in fib.values()).items())}
    out = {
        "program": "aut116.py",
        "sylow_fibres": fibres,
        "code_sha256": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
        "class_matrix_sha256": hashlib.sha256(open(os.path.join(HERE, "data", "class_matrix.gp"), "rb").read()).hexdigest(),
        "equitable_partition": eq,
        "order": order, "base": base, "orbit_lengths_along_chain": lengths, "generators": [list(map(int, g)) for g in gens],
        "order_after_random_relabelling": order2, "base_after_relabelling": base2,
        "orbit_lengths_along_chain_after_relabelling": lengths2,
        "element_orders": dict(sorted(collections.Counter(ords.values()).items())),
        "derived_order": len(D), "derived_perfect": perfect, "derived_simple": simple,
        "derived_classes": sorted([(o, s) for (o, s, x) in classes_D]),
        "centraliser_of_derived_order": len(centraliser),
        "index2_subgroups": [{"order": s["order"], "element_orders": s["element_orders"],
                              "fixes_small_points": s["set"] is Mset} for s in sub2],
        "point_orbits": pts, "point_orbit_names": {k: sorted(v for v in range(n) if name[v] == k) for k in "ABPQ"},
        "point_orbits_under_M": [len(o) for o in pts_M],
        "stabiliser_orders": {name[o[0]]: order // len(o) for o in pts},
        "eps_of_generators": [eps(g) for g in T],
        "index2_names": names2, "characters_on_generators_[M10,S6,PGL]": chargen,
        "galois_relabellings": galois,
        "orbitals": orbitals,
        "orbitals_count": len(orbitals),
        "contact_orbits": len(types), "contact_orbits_split": sum(1 for o in types if o["split"]),
        "seconds": round(time.time() - t0, 1)}
    os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, "out", "aut116.json"), "w"))
    with open(os.path.join(HERE, "data", "aut_generators.gp"), "w") as f:
        f.write("\\\\ Generators of Aut(X*) as permutations of the points 0..115 (entry i+1 = image of point i).\n")
        f.write("\\\\ Written by aut116.py. EPSGEN: their signs under the character with kernel M (= M10).\n")
        f.write("AUTGEN = [" + ", ".join("[" + ",".join(map(str, g)) + "]" for g in gens) + "];\n")
        f.write("EPSGEN = [" + ",".join(str(eps(g)) for g in T) + "];\n")
        f.write("\\\\ CHARGEN[i] = values on generator i of the characters with kernels M10, S6, PGL(2,9).\n")
        f.write("CHARGEN = [" + "; ".join(", ".join(map(str, cg)) for cg in chargen) + "];\n")
        f.write("\\\\ GALPERM[j] = [class map as [k, image], relabelling pi with C[pi(a), pi(b)] = map(C[a, b])].\n")
        f.write("GALPERM = [" + ", ".join("[[" + "; ".join("%s, %d" % (k, v) for k, v in gl["class_map"].items()) + "], ["
                                         + ",".join(map(str, gl["relabelling"])) + "]]" for gl in galois) + "];\n")
    with open(os.path.join(HERE, "data", "orbitals.gp"), "w") as f:
        f.write("\\\\ The orbits of Aut(X*) on unordered pairs of points (0-based): [class, size, split under M, names,\n"
                "\\\\ [[a, b, sign], ...]]; sign = +1 / -1 on the two halves of a split orbit, 0 otherwise. aut116.py.\n")
        f.write("ORBS = List();\n")
        for o in orbitals:
            f.write("listput(ORBS, [%d, %d, %d, \"%s\", [%s]]);\n" % (
                o["class"], o["size"], int(o["split"]), o["orbits"],
                ", ".join("[%d,%d,%d]" % tuple(q) for q in o["pairs"])))
    brief = {k: v for k, v in out.items() if k not in ("orbitals", "generators", "point_orbits", "point_orbit_names")}
    brief["orbitals_summary"] = [(o["class"], o["size"], o["orbits"], o["split"]) for o in orbitals]
    print(json.dumps(brief))


if __name__ == "__main__":
    main()
