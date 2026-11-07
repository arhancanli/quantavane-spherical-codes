# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""A second, separately written check that the model of Definition 3.1 gives the deposited class matrices of X* and X'.

  python3 -I certify/note/model_check.py OUT.json

Written from the text of the note alone (Section 3, Table 2, Figure 1 and Section 7); it shares no code with
construct.py, aut116.py or the GP scripts, and uses only Python integers and lists (exact).

 1. F9 = F3[z]/(z^2 - z - 1), the element a + b z encoded as a + 3b; P^1 = F9 u {inf}. PGammaL(2,9) is built as the set
    of maps z -> (a1 phi(z) + a2)/(a3 phi(z) + a4), a1 a4 != a2 a3, phi = identity or Frobenius (1440 distinct maps);
    PSL(2,9) = the maps with phi = identity and a square determinant, M10 = PSL(2,9) together with the Frobenius maps
    of non-square determinant.
 2. The 30 circles (images of F3 u {inf}), the 4 pencils of 3 pairwise tangent circles at each point, the
    perpendicularity of pencils at inf ((d'/d)^2 = -1 for the directions d, d'), transported to every point (checked
    to be well defined), the 72 pentads (orbits of the 144 elements of order 5), the 40 pairs beta, the A6-orbits on A
    and on B (two each, preserved by M10), the type-I relation, and A+, B+ fixed by Figure 1: the pentad
    O = {inf, 0, z, 1+2z, 2+2z} lies in A+ and beta0 = (inf, {pencils of directions 1 and z}) in B+.
 3. The class matrix K of Table 2 on the 116 objects (A, then B, then p+, q+, p-, q-), checked symmetric, preserved by
    all 1440 maps, with the class counts of Table 1; and the 360 pairs (O, beta) with both pencils of type I for O
    split as 180 + 180 between A+ x B+ and A- x B-.
 4. An exhaustive individualisation-refinement search for a bijection pi with C[pi(a), pi(b)] = K[a, b] on all pairs,
    where C is data/class_matrix.gp (the matrix that both certificates verify); the equality is then checked on all
    6,670 pairs. The images of p+, q+, p-, q- are reported.
 5. Section 7: K' = K with the classes of beta and q+, q- exchanged (12 if beta in B_{-eta}, 11 otherwise) equals
    data_xprime/class_matrix.gp under the same bijection pi, on all pairs.
 6. The number of permutations preserving K, by the same exhaustive search (all leaves).
"""
import itertools, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))


def f9(a, b):
    return (a % 3) + 3 * (b % 3)


def comps(x):
    return x % 3, x // 3


def add(x, y):
    a, b = comps(x); c, d = comps(y)
    return f9(a + c, b + d)


def neg(x):
    a, b = comps(x)
    return f9(-a, -b)


def mul(x, y):
    a, b = comps(x); c, d = comps(y)
    return f9(a * c + b * d, a * d + b * c + b * d)        # z^2 = z + 1


INV = {x: y for x in range(1, 9) for y in range(1, 9) if mul(x, y) == 1}
ZETA, INF, PTS = f9(0, 1), 9, list(range(10))
SQUARES = {mul(x, x) for x in range(1, 9)}
MINUS1 = f9(2, 0)


def frob(x):
    return mul(x, mul(x, x))


def mobius(a1, a2, a3, a4, phi):
    out = []
    for z in PTS:
        if z == INF:
            out.append(INF if a3 == 0 else mul(a1, INV[a3]))
        else:
            w = frob(z) if phi else z
            num, den = add(mul(a1, w), a2), add(mul(a3, w), a4)
            out.append(INF if den == 0 else mul(num, INV[den]))
    return tuple(out)


def need(cond, msg):
    if not cond:
        raise SystemExit("FAILED: " + msg)


def build_model():
    assert mul(mul(ZETA, ZETA), mul(ZETA, ZETA)) == MINUS1        # z^4 = -1
    G = {}
    for a1, a2, a3, a4 in itertools.product(range(9), repeat=4):
        det = add(mul(a1, a4), neg(mul(a2, a3)))
        if det == 0:
            continue
        for phi in (0, 1):
            p = mobius(a1, a2, a3, a4, phi)
            info = (phi, det in SQUARES)
            need(G.get(p, info) == info, "a map with two determinant classes")
            G[p] = info
    GROUP = sorted(G)
    need(len(GROUP) == 1440, "|PGammaL(2,9)| != 1440")
    PSL = [g for g in GROUP if G[g] == (0, True)]
    M10 = [g for g in GROUP if G[g] in [(0, True), (1, False)]]
    inM10 = {g: (G[g] in [(0, True), (1, False)]) for g in GROUP}
    need(len(PSL) == 360 and len(M10) == 720, "PSL(2,9) or M10 of the wrong order")
    base = frozenset([f9(0, 0), f9(1, 0), f9(2, 0), INF])
    CIRCLES = sorted({frozenset(g[x] for x in base) for g in GROUP}, key=sorted)
    need(len(CIRCLES) == 30 and all(len(C) == 4 for C in CIRCLES), "not 30 circles of 4 points")
    PENCILS = {}
    for z in PTS:
        through = [C for C in CIRCLES if z in C]
        need(len(through) == 12, "not 12 circles through a point")
        pens = set()
        for C in through:
            pens.add(frozenset([C] + [D for D in through if D != C and C & D == {z}]))
        pens = sorted(pens, key=lambda P: sorted(sorted(C) for C in P))
        need(len(pens) == 4 and all(len(P) == 3 for P in pens), "not 4 pencils of 3 circles")
        for P in pens:
            for C, D in itertools.combinations(P, 2):
                need(C & D == {z}, "circles of a pencil not tangent")
        PENCILS[z] = pens

    def direction(P):
        C = next(iter(P))
        aff = sorted(x for x in C if x != INF)
        d = add(aff[1], neg(aff[0]))
        return min(d, neg(d))

    def perp_inf(P, Q):
        r = mul(direction(Q), INV[direction(P)])
        return mul(r, r) == MINUS1

    def gpencil(g, P):
        return frozenset(frozenset(g[x] for x in C) for C in P)

    PERP = {}
    for z in PTS:
        pairs = None
        for g in GROUP:
            if g[INF] != z:
                continue
            s = set()
            for P, Q in itertools.combinations(PENCILS[INF], 2):
                if perp_inf(P, Q):
                    s.add(frozenset([gpencil(g, P), gpencil(g, Q)]))
            if pairs is None:
                pairs = s
            need(pairs == s, "perpendicularity not well defined")
        need(len(pairs) == 2, "not two perpendicular pairs at a point")
        PERP[z] = pairs
    for g in GROUP:
        if g[INF] == INF:
            for P, Q in itertools.combinations(PENCILS[INF], 2):
                need(perp_inf(P, Q) == perp_inf(gpencil(g, P), gpencil(g, Q)), "stabiliser of inf breaks perpendicularity")
    Bset = []
    for z in PTS:
        for P, Q in itertools.combinations(PENCILS[z], 2):
            if frozenset([P, Q]) not in PERP[z]:
                Bset.append((z, frozenset([P, Q])))
    need(len(Bset) == 40, "not 40 pairs beta")
    ident = tuple(PTS)

    def order(g):
        k, h = 1, g
        while h != ident:
            h = tuple(g[h[i]] for i in PTS)
            k += 1
        return k
    pentads, n5 = set(), 0
    for g in GROUP:
        if order(g) == 5:
            n5 += 1
            seen = set()
            for z in PTS:
                if z in seen:
                    continue
                orb, w = [], z
                while w not in orb:
                    orb.append(w)
                    w = g[w]
                seen |= set(orb)
                need(len(orb) == 5, "an element of order 5 with a fixed point")
                pentads.add(frozenset(orb))
    need(n5 == 144, "not 144 elements of order 5")
    Aset = sorted(pentads, key=sorted)
    need(len(Aset) == 72, "not 72 pentads")

    def gB(g, b):
        z, pr = b
        return (g[z], frozenset(gpencil(g, P) for P in pr))

    def gA(g, O):
        return frozenset(g[x] for x in O)

    def orbits(objs, act, grp):
        idx = {o: i for i, o in enumerate(objs)}
        lab, k = [-1] * len(objs), 0
        for i, o in enumerate(objs):
            if lab[i] >= 0:
                continue
            for g in grp:
                lab[idx[act(g, o)]] = k
            k += 1
        return lab
    A6orbA, A6orbB = orbits(Aset, gA, PSL), orbits(Bset, gB, PSL)
    need(sorted(set(A6orbA)) == [0, 1] and sorted(set(A6orbB)) == [0, 1], "A6 does not have two orbits on A and on B")
    need(orbits(Aset, gA, M10) == A6orbA and orbits(Bset, gB, M10) == A6orbB, "M10 does not preserve the A6-orbits")

    def typeI(P, O):
        return sum(1 for C in P if len(C & O) == 3) == 1
    Ofig = frozenset([INF, f9(0, 0), f9(0, 1), f9(1, 2), f9(2, 2)])
    need(Ofig in pentads, "the pentad of Figure 1 is not a pentad")

    def pencil_dir(d):
        return [P for P in PENCILS[INF] if direction(P) == min(d, neg(d))][0]
    b0 = (INF, frozenset([pencil_dir(f9(1, 0)), pencil_dir(ZETA)]))
    need(b0 in Bset and all(typeI(P, Ofig) for P in b0[1]), "beta0 of Figure 1")
    Aplus, Bplus = A6orbA[Aset.index(Ofig)], A6orbB[Bset.index(b0)]
    Asign = [(1 if l == Aplus else -1) for l in A6orbA]
    Bsign = [(1 if l == Bplus else -1) for l in A6orbB]
    pairs360 = [(i, j) for i, O in enumerate(Aset) for j, b in enumerate(Bset)
                if b[0] in O and all(typeI(P, O) for P in b[1])]
    split = (len(pairs360), sum(1 for i, j in pairs360 if Asign[i] == Bsign[j] == 1),
             sum(1 for i, j in pairs360 if Asign[i] == Bsign[j] == -1))
    need(split == (360, 180, 180), "the 360 pairs do not split as 180 + 180 in A+ x B+ and A- x B-")
    return GROUP, inM10, PENCILS, Aset, Bset, A6orbA, Asign, Bsign, gA, gB, typeI


def class_matrix(model, xprime=False):
    GROUP, inM10, PENCILS, Aset, Bset, A6orbA, Asign, Bsign, gA, gB, typeI = model
    NA, NB, N = 72, 40, 116
    iP, iQ = {1: 112, -1: 114}, {1: 113, -1: 115}

    def move(P, z2):
        C = [C for C in P if z2 in C]
        need(len(C) == 1, "moving a pencil")
        Q = [Q for Q in PENCILS[z2] if C[0] in Q]
        need(len(Q) == 1, "moving a pencil")
        return Q[0]
    K = [[-1] * N for _ in range(N)]
    for i, O in enumerate(Aset):
        for j, O2 in enumerate(Aset):
            if i != j:
                m = len(O & O2)
                need(m in (0, 1, 2, 3, 4), "pentads meeting in 5 points")
                K[i][j] = 12 if m in (0, 2) else (5 if m in (1, 4) else (12 if A6orbA[i] == A6orbA[j] else 5))
    for i, (z, pr) in enumerate(Bset):
        for j, (z2, pr2) in enumerate(Bset):
            if i != j:
                if z == z2:
                    need(len(pr & pr2) < 2, "equal pairs of pencils")
                    K[NA + i][NA + j] = 8 if not (pr & pr2) else 2
                else:
                    K[NA + i][NA + j] = 6 if frozenset(move(P, z2) for P in pr) == pr2 else 12
    for i, O in enumerate(Aset):
        for j, b in enumerate(Bset):
            k = {2: 3, 0: 10, 1: 12}[sum(1 for P in b[1] if typeI(P, O))]
            K[i][NA + j] = K[NA + j][i] = k
    for eta in (1, -1):
        for i in range(NA):
            K[i][iP[eta]] = K[iP[eta]][i] = 9 if Asign[i] == eta else 12
            K[i][iQ[eta]] = K[iQ[eta]][i] = 7 if Asign[i] == eta else 12
        for j in range(NB):
            K[NA + j][iP[eta]] = K[iP[eta]][NA + j] = 12 if Bsign[j] == eta else 4
            qs = -eta if xprime else eta                      # Section 7: the classes with q+ and q- exchanged
            K[NA + j][iQ[eta]] = K[iQ[eta]][NA + j] = 12 if Bsign[j] == qs else 11
        for eta2 in (1, -1):
            K[iP[eta]][iQ[eta2]] = K[iQ[eta2]][iP[eta]] = 1 if eta == eta2 else 0
    K[iP[1]][iP[-1]] = K[iP[-1]][iP[1]] = 12
    K[iQ[1]][iQ[-1]] = K[iQ[-1]][iQ[1]] = 12
    need(all(K[a][b] == K[b][a] and (K[a][b] >= 0) == (a != b) for a in range(N) for b in range(N)), "K not symmetric")
    return K


def group_perms(model):
    GROUP, inM10, PENCILS, Aset, Bset, A6orbA, Asign, Bsign, gA, gB, typeI = model
    Aidx = {O: i for i, O in enumerate(Aset)}
    Bidx = {b: i for i, b in enumerate(Bset)}
    out = []
    for g in GROUP:
        p = [0] * 116
        for i, O in enumerate(Aset):
            p[i] = Aidx[gA(g, O)]
        for j, b in enumerate(Bset):
            p[72 + j] = 72 + Bidx[gB(g, b)]
        sw = not inM10[g]
        for x, y in ((112, 114), (113, 115)):
            p[x], p[y] = (y, x) if sw else (x, y)
        out.append(p)
    return out


def refine(K1, c1, K2, c2, N):
    while True:
        n1 = len(set(c1))
        sig1 = [(c1[v], tuple(sorted((K1[v][u], c1[u]) for u in range(N) if u != v))) for v in range(N)]
        sig2 = [(c2[v], tuple(sorted((K2[v][u], c2[u]) for u in range(N) if u != v))) for v in range(N)]
        if sorted(sig1) != sorted(sig2):
            return None, None
        keys = {s: i for i, s in enumerate(sorted(set(sig1)))}
        c1, c2 = [keys[s] for s in sig1], [keys[s] for s in sig2]
        if len(set(c1)) == n1:
            return c1, c2


def search(K1, c1, K2, c2, N, count_all=False):
    """Bijections perm with K2[perm[a]][perm[b]] = K1[a][b] (the first one, or all of them)."""
    c1, c2 = refine(K1, c1, K2, c2, N)
    if c1 is None:
        return []
    if len(set(c1)) == N:
        pos2 = {c: v for v, c in enumerate(c2)}
        perm = [pos2[c1[v]] for v in range(N)]
        ok = all(K2[perm[a]][perm[b]] == K1[a][b] for a in range(N) for b in range(N) if a != b)
        return [perm] if ok else []
    sizes = {}
    for c in c1:
        sizes[c] = sizes.get(c, 0) + 1
    cell = min((s, c) for c, s in sizes.items() if s > 1)[1]
    v = min(u for u in range(N) if c1[u] == cell)
    out = []
    for w in [u for u in range(N) if c2[u] == cell]:
        new = max(c1) + 1
        d1, d2 = list(c1), list(c2)
        d1[v], d2[w] = new, new
        res = search(K1, d1, K2, d2, N, count_all)
        out += res
        if res and not count_all:
            return out
    return out


def read_gp_matrix(path):
    rows = [json.loads(line) for line in open(path) if line.strip()]
    need(len(rows) == 116 and all(len(r) == 116 for r in rows), "class matrix of the wrong size: " + path)
    return rows


def main():
    out = sys.argv[1]
    t0 = time.time()
    model = build_model()
    K = class_matrix(model)
    N = 116
    counts = [sum(1 for a in range(N) for b in range(a + 1, N) if K[a][b] == k) for k in range(13)]
    perms = group_perms(model)
    need(all(K[p[a]][p[b]] == K[a][b] for p in perms for a in range(N) for b in range(N) if a != b),
         "the group does not preserve K")
    C = read_gp_matrix(os.path.join(HERE, "data", "class_matrix.gp"))
    Cp = read_gp_matrix(os.path.join(HERE, "data_xprime", "class_matrix.gp"))
    iso = search(K, [0] * N, C, [0] * N, N)
    need(len(iso) == 1, "no bijection between the model and data/class_matrix.gp")
    pi = iso[0]
    pairs_equal = sum(1 for a in range(N) for b in range(a + 1, N) if C[pi[a]][pi[b]] == K[a][b])
    need(pairs_equal == N * (N - 1) // 2, "the bijection fails on some pair")
    Kp = class_matrix(model, xprime=True)
    pairs_equal_xp = sum(1 for a in range(N) for b in range(a + 1, N) if Cp[pi[a]][pi[b]] == Kp[a][b])
    need(pairs_equal_xp == N * (N - 1) // 2, "the model of Section 7 differs from data_xprime/class_matrix.gp")
    t1 = time.time()
    auts = search(K, [0] * N, K, [0] * N, N, count_all=True)
    res = {"program": "model_check.py", "group_order": len(perms), "class_counts": counts,
           "split_360": [360, 180, 180], "bijection_found": True, "pairs_checked": pairs_equal,
           "PQ_images_in_stored_order": [pi[i] for i in (112, 113, 114, 115)],
           "xprime_pairs_checked": pairs_equal_xp, "class_preserving_permutations": len(auts),
           "seconds_model_and_bijection": round(t1 - t0, 1), "seconds_count": round(time.time() - t1, 1),
           "status": "MODEL = DEPOSITED CLASS MATRICES (X* and X')"}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
