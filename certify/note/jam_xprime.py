# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""The second code X' of the note (Theorem 7.1), checked by the Python/Arb path: an adapted copy of
certify/v5/jam116.py (the first verification of the note, Section 6), independent of the PARI/GP programs.

  python3 -I certify/note/jam_xprime.py OUT.json [--red gram|ratio]

X' is X* with the two points of Q reflected in the hyperplane spanned by the orbit A. Nothing computed by PARI/GP
enters a trusted step: X* is read from the version 5 recognition (certify/v5/out/recognition_d22_n116.json, the data
that jam116.py certifies), the orbits from out/aut116.json (Python), and X' is computed here.

TRUSTED and exact (arithmetic in L = Q(sqrt145, s1, s2) of degree 8 with certify/v5/tower.py; signs at the real
embedding in FLINT ball arithmetic; ranks over a prime field):
 (L) [L:Q] = 8 and the embedding is real.
 (X) the Gram matrix G' of X' from that of X*: 21 points LA of A with G[LA,LA] invertible (exact Gauss-Jordan) and
     G[a,a'] = G[a,LA] G[LA,LA]^-1 G[LA,a'] for all a, a' in A, so A spans a hyperplane n_A^perp of R^22 and
     <u,n_A><v,n_A> = G[u,v] - G[u,LA] G[LA,LA]^-1 G[LA,v]; the reflection changes only the entries between Q and the
     other points: G'[u,q] = G[u,q] - 2<u,n_A><q,n_A>. Every entry of G' is one of thirteen values; the two values on
     the pairs between P and Q are not values of X* (so X' is not congruent to X*); and G' agrees entry by entry with
     data_xprime/ (written by xprime.gp), a cross-check of the two programs.
 (C) c = t_12 satisfies 30c^2 + 5c - 1 = 0, every other value is below c, every point lies in at least 22 contacts.
 (G) G' is the Gram matrix of 116 unit vectors spanning R^22: G'_BB positive definite (exact pivots) and
     G' = G'_{.B} G'_BB^{-1} G'_{B.} exactly (B the basis of the note, which contains no point of Q).
 (R) the bar rigidity matrix of (X', E) has rank at least 116*22 - 231 - 1 over a prime field: modulo rotations the
     infinitesimal flexes form a space of dimension at most one.
 (F) a nonzero first-order Gram change S of a flex, zero on the diagonal and the contacts, solved exactly from sampled
     equations and checked entry by entry: X' is not infinitesimally jammed.
 (S) an equilibrium stress, constant on the orbits of contacts, checked exactly at every point, every entry positive.
 (H) H_omega(F) > 0 exactly.
By the note's Theorem 2.5, (G), (C), (R), (F), (S), (H) make X' a strict local minimum of the largest cosine modulo
rotations, so jammed, and not infinitesimally jammed. There is no stored code for X', so jam116.py's step (D) is left
out. UNTRUSTED (floating point, proposals only): the choice of LA, the sign pattern of a numerical flex, the stress
nearest to a linear-programming stress.

Red controls (each must FAIL): --red gram checks the configuration with the classes of X' but the values t0, t1 of X* on
the pairs between P and Q (not a Gram matrix of rank 22: (G) must fail); --red ratio proposes a strictly positive stress
with omega(q+q-)/omega(p+p-) = 0.2, below the threshold of Theorem 7.1 ((H) must fail).
"""
import argparse, hashlib, json, math, os, sys, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
V5 = os.path.join(os.path.dirname(HERE), "v5")
sys.path[:0] = [p for p in os.environ.get("QV_LIBS", "").split(os.pathsep) if p]   # optional library folders
sys.path.insert(0, V5)
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.seterr(all="ignore")
import flint
from kfield import K, fq
from tower import Tower
import exact2 as X2

SHA = {"jam_xprime.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()}
for _f in ("tower.py", "kfield.py", "exact2.py"):
    SHA[_f] = hashlib.sha256(open(os.path.join(V5, _f), "rb").read()).hexdigest()


def fail(log, out, msg):
    log["status"] = "FAILED: " + msg
    json.dump(log, open(out, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in log.items() if k != "code_sha256"}, default=str)[:3000])
    raise SystemExit(msg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--red", choices=["gram", "ratio"], default=None)
    args = ap.parse_args()
    out = args.out
    t0 = time.time()
    recf = os.path.join(V5, "out", "recognition_d22_n116.json")
    rec = json.load(open(recf))
    aut = json.load(open(os.path.join(HERE, "out", "aut116.json")))
    d, n = rec["d"], rec["n"]
    g = d * (d - 1) // 2
    log = {"d": d, "n": n, "program": "jam_xprime.py", "code_sha256": SHA, "red_control": args.red,
           "recognition_sha256": hashlib.sha256(open(recf, "rb").read()).hexdigest()}
    # ------------------------------------------------------------------ (L) the field (as jam116.py)
    def k0(r):
        a, b, c = r
        return K(fq(-b, a), fq(-c, a))
    K.D = 145
    T = Tower(145, k0(rec["Delta1_rel"]), k0(rec["Delta2_rel"]), rec["sigma1"], rec["sigma2"])
    facts = T.check_degree()
    log["field"] = {k: (v if isinstance(v, bool) else str(v)) for k, v in facts.items()}
    if not (facts["degree_8"] and facts["Delta1_positive"] and facts["Delta2_positive"]):
        fail(log, out, "field check failed")

    def from_rel(r):
        r0 = r[0]
        co = [fq(-x, r0) for x in r[1:]]
        return ((K(co[0], co[1]), K(co[2], co[3])), (K(co[4], co[5]), K(co[6], co[7])))
    tstar = [from_rel(r) for r in rec["class_rel"]]
    nc = len(tstar)
    ONE, ZERO = T.const(1), T.const(0)
    ci = max(range(nc), key=lambda k: rec["classes"][k])
    pc = rec["pair_class"]
    cls_star = {}
    idx = 0
    for a in range(n):
        for b in range(a + 1, n):
            cls_star[(a, b)] = pc[idx]
            idx += 1

    def Gs(a, b):
        if a == b:
            return ONE
        return tstar[cls_star[(a, b)] if a < b else cls_star[(b, a)]]
    # ------------------------------------------------------------------ (X) the Gram matrix of X'
    A = aut["point_orbit_names"]["A"]
    Pp = aut["point_orbit_names"]["P"]
    Qq = aut["point_orbit_names"]["Q"]
    # untrusted: 21 points of A spanning span(A), greedily in floating point
    Gf = np.array([[T.to_float(Gs(a, b)) for b in range(n)] for a in range(n)])
    LA, dprev = [], 1.0
    for a in A:
        if len(LA) == 21:
            break
        tst = LA + [a]
        dn = np.linalg.det(Gf[np.ix_(tst, tst)])
        if dn / dprev > 1e-6:
            LA, dprev = tst, dn
    if len(LA) != 21:
        fail(log, out, "no 21 spanning points of A")
    m = len(LA)
    Mx = [[Gs(LA[i], LA[j]) for j in range(m)] + [ONE if i == j else ZERO for j in range(m)] for i in range(m)]
    for k in range(m):
        pk = Mx[k][k]
        if T.iszero(pk):
            fail(log, out, "G[LA,LA] singular")
        inv = T.inv(pk)
        Mx[k] = [T.mul(x, inv) for x in Mx[k]]
        for i in range(m):
            if i != k and not T.iszero(Mx[i][k]):
                fct = Mx[i][k]
                Mx[i] = [T.sub(x, T.mul(fct, y)) for x, y in zip(Mx[i], Mx[k])]
    GLAi = [row[m:] for row in Mx]
    VLA = [[None] * n for _ in range(m)]            # G[LA,LA]^-1 G[LA, v]
    for i in range(m):
        for v in range(n):
            acc = ZERO
            for j in range(m):
                gj = Gs(LA[j], v)
                if not T.iszero(GLAi[i][j]) and not T.iszero(gj):
                    acc = T.add(acc, T.mul(GLAi[i][j], gj))
            VLA[i][v] = acc

    def proj(u, v):                                  # <u, Pi_A v>
        acc = ZERO
        for i in range(m):
            gu = Gs(u, LA[i])
            if not T.iszero(gu) and not T.iszero(VLA[i][v]):
                acc = T.add(acc, T.mul(gu, VLA[i][v]))
        return acc
    badA = sum(1 for i, a in enumerate(A) for a2 in A[i:] if not T.eq(proj(a, a2), Gs(a, a2)))
    log["A_in_span_of_LA_failures"] = badA
    if badA:
        fail(log, out, "A does not lie in the span of LA")
    Gp = {}
    for u in range(n):
        if u in Qq:
            continue
        for q in Qq:
            mn = T.sub(Gs(u, q), proj(u, q))            # <u, n_A><q, n_A>
            Gp[(min(u, q), max(u, q))] = T.sub(Gs(u, q), T.add(mn, mn))

    def Gx(a, b):                                     # the Gram matrix of X'
        if a == b:
            return ONE
        key = (a, b) if a < b else (b, a)
        return Gp[key] if key in Gp else Gs(a, b)
    # the values of X': those of X*, with the two new values on the pairs between P and Q
    pq_vals = []
    for p_ in Pp:
        for q in Qq:
            v = Gx(p_, q)
            if not any(T.eq(v, w) for w in pq_vals):
                pq_vals.append(v)
    if len(pq_vals) != 2:
        fail(log, out, "expected two values on the pairs between P and Q")
    pq_vals.sort(key=lambda v: T.to_float(v))
    new_vals = all(not any(T.eq(v, w) for w in tstar) for v in pq_vals)
    log["new_values_not_values_of_Xstar"] = new_vals
    if not new_vals:
        fail(log, out, "a value between P and Q is a value of X*")
    tvals = [pq_vals[0], pq_vals[1]] + tstar[2:]
    if not all(T.lt(tvals[k], tvals[k + 1]) for k in range(nc - 1)):
        fail(log, out, "the values of X' are not increasing")
    cls_of = {}
    for a in range(n):
        for b in range(a + 1, n):
            v = Gx(a, b)
            k = next((k for k in range(nc) if T.eq(v, tvals[k])), None)
            if k is None:
                fail(log, out, "an entry of G' is not one of the thirteen values")
            cls_of[(a, b)] = k
    changed = sum(1 for key in cls_of if cls_of[key] != cls_star[key])
    log.update(LA=LA, X_entries_checked=len(cls_of), classes_changed_from_Xstar=changed,
               new_values_float=[T.to_float(v) for v in pq_vals],
               new_values_tower=[[str(x) for y in v for z in y for x in (z.a, z.b)] for v in pq_vals])
    # cross-check with the data of xprime.gp (PARI): the same values and the same class matrix
    import re
    vt = open(os.path.join(HERE, "data_xprime", "values.gp")).read()
    TC = [[Fraction(x.strip()) for x in row.split(",")] for row in re.search(r"^TC = \[(.*?)\];", vt, re.M).group(1).split(";")]
    def tw(c):
        return ((K(fq(c[0].numerator, c[0].denominator), fq(c[1].numerator, c[1].denominator)),
                 K(fq(c[2].numerator, c[2].denominator), fq(c[3].numerator, c[3].denominator))),
                (K(fq(c[4].numerator, c[4].denominator), fq(c[5].numerator, c[5].denominator)),
                 K(fq(c[6].numerator, c[6].denominator), fq(c[7].numerator, c[7].denominator))))
    vals_gp = [tw(c) for c in TC]
    cm_gp = [json.loads(line) for line in open(os.path.join(HERE, "data_xprime", "class_matrix.gp")) if line.strip()]
    agree = (all(T.eq(x, y) for x, y in zip(vals_gp, tvals)) and
             all(cm_gp[a][b] == cls_of[(a, b)] for (a, b) in cls_of))
    log["agrees_with_xprime_gp"] = agree
    if not agree:
        fail(log, out, "X' differs from data_xprime (xprime.gp)")
    log["seconds_X"] = round(time.time() - t0, 1)
    if args.red == "gram":
        # RED: the classes of X' with the values t0, t1 of X* on the pairs between P and Q
        tvals = tstar[:2] + tvals[2:]
        log["red_gram"] = "values t0, t1 of X* on the pairs between P and Q"

    def Gv(a, b):
        if a == b:
            return ONE
        return tvals[cls_of[(a, b)] if a < b else cls_of[(b, a)]]
    # ------------------------------------------------------------------ (C) the cosine (as jam116.py)
    c = tvals[ci]
    poly = T.add(T.add(T.mul(T.const(30), T.mul(c, c)), T.mul(T.const(5), c)), T.const(-1))
    below = [k for k in range(nc) if k != ci and T.lt(tvals[k], c)]
    log.update(cos_class=ci, cos_value_float=T.to_float(c), cos_satisfies_30x2_5x_1=T.iszero(poly),
               classes_below_c=len(below), classes=nc)
    if not (T.iszero(poly) and len(below) == nc - 1):
        fail(log, out, "cosine check failed")
    E = [(a, b) for (a, b), k in cls_of.items() if k == ci]
    deg = [0] * n
    for a, b in E:
        deg[a] += 1
        deg[b] += 1
    log.update(contacts=len(E), min_contacts_per_point=min(deg), max_contacts_per_point=max(deg),
               rattlers=sum(1 for x in deg if x == 0))
    if min(deg) < d:
        fail(log, out, "a point has fewer than d contacts")
    # ------------------------------------------------------------------ (G) the Gram matrix (as jam116.py)
    B = rec["basis"]
    if any(b in Qq for b in B):
        fail(log, out, "the basis contains a point of Q")
    GBB = [[Gv(a, b) for b in B] for a in B]
    Mx = [row[:] + [ONE if i == j else ZERO for j in range(d)] for i, row in enumerate(GBB)]
    pivots = []
    for k in range(d):
        pk = Mx[k][k]
        if T.iszero(pk):
            fail(log, out, "zero pivot")
        pivots.append(pk)
        inv = T.inv(pk)
        Mx[k] = [T.mul(x, inv) for x in Mx[k]]
        for i in range(d):
            if i != k and not T.iszero(Mx[i][k]):
                fct = Mx[i][k]
                Mx[i] = [T.sub(x, T.mul(fct, y)) for x, y in zip(Mx[i], Mx[k])]
    GBBinv = [row[d:] for row in Mx]
    pos_piv = all(T.sign(p) > 0 for p in pivots)
    log["G_BB_positive_definite"] = pos_piv
    if not pos_piv:
        fail(log, out, "G_BB not positive definite")
    W = [[None] * n for _ in range(d)]
    for k in range(d):
        for v in range(n):
            acc = ZERO
            for l in range(d):
                acc = T.add(acc, T.mul(GBBinv[k][l], Gv(B[l], v)))
            W[k][v] = acc
    bad = 0
    for a in range(n):
        ga = [Gv(a, b) for b in B]
        for v in range(a, n):
            acc = ZERO
            for l in range(d):
                if not T.iszero(ga[l]) and not T.iszero(W[l][v]):
                    acc = T.add(acc, T.mul(ga[l], W[l][v]))
            if not T.eq(acc, Gv(a, v)):
                bad += 1
    log["gram_rank_d_check_failures"] = bad
    if bad:
        fail(log, out, "G does not have rank d")
    log["seconds_gram"] = round(time.time() - t0, 1)
    # ------------------------------------------------------------------ (R) rank over a prime field (as jam116.py)
    p = (1 << 31) - 1
    setup = None
    while setup is None:
        p -= 2
        if flint.fmpz(p).is_prime():
            setup = T.modp_setup(p)
    vmod = {}

    def gm(a, b):
        key = (a, b) if a <= b else (b, a)
        if key not in vmod:
            vmod[key] = T.modp(Gv(a, b), setup)
        return vmod[key]
    rows = []
    for a in range(n):
        rows.append({a * d + k: gm(a, B[k]) for k in range(d)})
    for a, b in E:
        row = {}
        for k in range(d):
            row[b * d + k] = (row.get(b * d + k, 0) + gm(a, B[k])) % p
            row[a * d + k] = (row.get(a * d + k, 0) + gm(b, B[k])) % p
        rows.append(row)
    rk = X2.rank_mod(rows, n * d, p, sketch_seed=2026)
    log.update(prime=p, rank_mod_p=int(rk), rank_needed=n * d - g - 1, cols=n * d)
    if rk < n * d - g - 1:
        fail(log, out, "rank bound too small")
    log["seconds_rank"] = round(time.time() - t0, 1)
    # ------------------------------------------------------------------ untrusted proposals in floating point
    Gfl = np.array([[T.to_float(Gv(a, b)) for b in range(n)] for a in range(n)])
    lam, U = np.linalg.eigh(Gfl)
    Xf = U[:, -d:] * np.sqrt(np.maximum(lam[-d:], 0))
    Rm = np.zeros((len(E) + n, n * d))
    for r_, (a, b) in enumerate(E):
        Rm[r_, a * d:(a + 1) * d] = Xf[b]
        Rm[r_, b * d:(b + 1) * d] = Xf[a]
    for a in range(n):
        Rm[len(E) + a, a * d:(a + 1) * d] = Xf[a]
    _, sv, Vt = np.linalg.svd(Rm, full_matrices=False)
    rkf = int((sv > 1e-8 * sv.max()).sum())
    Nul = Vt[rkf:].T
    rots = []
    for i in range(d):
        for j in range(i + 1, d):
            Am = np.zeros((d, d)); Am[i, j], Am[j, i] = 1, -1
            rots.append((Xf @ Am.T).ravel())
    Qr, _ = np.linalg.qr(np.array(rots).T)
    Fl = Nul - Qr @ (Qr.T @ Nul)
    U2, s2, _ = np.linalg.svd(Fl, full_matrices=False)
    yfl = U2[:, 0].reshape(n, d)
    Sf = Xf @ yfl.T + yfl @ Xf.T
    smax = np.abs(Sf).max()
    fcls = sorted({cls_of[(a, b)] for (a, b) in cls_of if abs(Sf[a, b]) > 1e-6 * smax})
    sgn = {(a, b): int(np.sign(Sf[a, b])) for (a, b) in cls_of if cls_of[(a, b)] in fcls and abs(Sf[a, b]) > 1e-6 * smax}
    log.update(proposal_flex_numerical_nullity_mod_rotations=int(sum(s2 > 1e-6)), flex_classes=fcls)
    # the orbits of contacts (aut116.py): edge types
    orb_of = {}
    for oi, o in enumerate(aut["orbitals"]):
        for a, b, _ in o["pairs"]:
            orb_of[(min(a, b), max(a, b))] = oi
    ctypes = sorted({orb_of[e] for e in E})
    et = [ctypes.index(orb_of[e]) for e in E]
    ntype = len(ctypes)
    tPP = next(j for j, oi in enumerate(ctypes) if aut["orbitals"][oi]["orbits"] == "PP")
    tQQ = next(j for j, oi in enumerate(ctypes) if aut["orbitals"][oi]["orbits"] == "QQ")
    # a linear-programming stress constant on the edge types (floating point)
    from scipy.optimize import linprog
    cf = float(T.to_float(c))
    EQ = np.zeros((n * d, ntype))
    for (a, b), j in zip(E, et):
        EQ[a * d:(a + 1) * d, j] += Xf[b] - cf * Xf[a]
        EQ[b * d:(b + 1) * d, j] += Xf[a] - cf * Xf[b]
    _, sv3, V3 = np.linalg.svd(EQ, full_matrices=True)
    Ns = V3[int((sv3 > 1e-8 * sv3.max()).sum()):].T
    cnt = np.bincount(et, minlength=ntype).astype(float)
    nz = Ns.shape[1]
    A_ub = np.hstack([-Ns, np.ones((ntype, 1))]); b_ub = np.zeros(ntype)
    if args.red == "ratio":
        # RED: omega_QQ <= 0.2 omega_PP (a strictly positive stress below the threshold)
        A_ub = np.vstack([A_ub, np.concatenate([Ns[tQQ] - 0.2 * Ns[tPP], [0]])[None, :]]); b_ub = np.concatenate([b_ub, [0]])
        log["red_ratio"] = "stress with omega(q+q-)/omega(p+p-) <= 0.2"
    cvec = np.zeros(nz + 1); cvec[-1] = -1
    lp = linprog(cvec, A_ub=A_ub, b_ub=b_ub, A_eq=np.concatenate([cnt @ Ns, [0]])[None, :], b_eq=[1.0],
                 bounds=[(None, None)] * (nz + 1), method="highs")
    if lp.status != 0:
        fail(log, out, "no linear-programming stress (proposal)")
    wlp_types = Ns @ lp.x[:nz]
    log["proposal_lp_min_times_E"] = float(lp.x[-1] * len(E))
    # ------------------------------------------------------------------ (F) the flex (as jam116.py)
    def S_basis(j, a, b):
        if a == b:
            return 0
        key = (a, b) if a < b else (b, a)
        if cls_of[key] != fcls[j]:
            return 0
        return sgn.get(key, 0)

    def M_of(Sfun):
        Mm = [[None] * d for _ in range(n)]
        SBB = [[Sfun(B[k], B[l]) for l in range(d)] for k in range(d)]
        for v in range(n):
            for k in range(d):
                acc = Sfun(B[k], v)
                acc = T.const(acc) if not isinstance(acc, tuple) else acc
                for l in range(d):
                    s = SBB[k][l]
                    if isinstance(s, tuple):
                        if T.iszero(s):
                            continue
                        acc = T.sub(acc, T.mul(T.mul(T.const(fq(1, 2)), s), W[l][v]))
                    elif s:
                        acc = T.sub(acc, T.mul(T.const(fq(s, 2)), W[l][v]))
                Mm[v][k] = acc
        return Mm

    def resid(Sfun, Mm, a, v):
        s = Sfun(a, v)
        acc = T.const(s) if not isinstance(s, tuple) else s
        for l in range(d):
            acc = T.sub(acc, T.add(T.mul(W[l][a], Mm[v][l]), T.mul(Mm[a][l], W[l][v])))
        return acc
    rng = np.random.default_rng(116)
    nonB = [v for v in range(n) if v not in set(B)]
    Ms = [M_of(lambda a, b, j=j: S_basis(j, a, b)) for j in range(len(fcls))]
    small = [v for v in nonB if deg[v] == min(deg)]
    order = [(v, v) for v in nonB] + [(s, v) for s in small for v in nonB if v != s]
    rest = [(a, b) for i, a in enumerate(nonB) for b in nonB[i + 1:] if (a, b) not in set(order)]
    rng.shuffle(rest)
    order += rest
    used = 0

    def nullspace(rows_in):
        rowsA_ = [r[:] for r in rows_in]
        piv_cols, r_ = [], 0
        for col in range(len(fcls)):
            pr = next((i for i in range(r_, len(rowsA_)) if not T.iszero(rowsA_[i][col])), None)
            if pr is None:
                continue
            rowsA_[r_], rowsA_[pr] = rowsA_[pr], rowsA_[r_]
            inv = T.inv(rowsA_[r_][col])
            rowsA_[r_] = [T.mul(x, inv) for x in rowsA_[r_]]
            for i in range(len(rowsA_)):
                if i != r_ and not T.iszero(rowsA_[i][col]):
                    fct = rowsA_[i][col]
                    rowsA_[i] = [T.sub(x, T.mul(fct, y)) for x, y in zip(rowsA_[i], rowsA_[r_])]
            piv_cols.append(col)
            r_ += 1
        return piv_cols, rowsA_[:r_]
    piv_cols, red = [], []
    while used < len(order):
        chunk = order[used:used + 200]
        used += len(chunk)
        new_rows = [[resid(lambda a, b, j=j: S_basis(j, a, b), Ms[j], a, v) for j in range(len(fcls))] for a, v in chunk]
        piv_cols, red = nullspace(red + new_rows)
        if len(fcls) - len(piv_cols) <= 1:
            break
    rowsA = red
    free = [col for col in range(len(fcls)) if col not in piv_cols]
    log["flex_sample_pairs"] = used
    log["flex_sample_nullity"] = len(free)
    if len(free) < 1:
        fail(log, out, "no flex of the given pattern")
    fcol = free[0]
    delta = [ZERO] * len(fcls)
    delta[fcol] = ONE
    for i, pcol in enumerate(piv_cols):
        delta[pcol] = T.neg(rowsA[i][fcol])

    def S_full(a, b):
        acc = ZERO
        for j in range(len(fcls)):
            s = S_basis(j, a, b)
            if s:
                acc = T.add(acc, T.mul(T.const(s), delta[j]))
        return acc
    Mfull = [[ZERO] * d for _ in range(n)]
    for j in range(len(fcls)):
        if T.iszero(delta[j]):
            continue
        for v in range(n):
            for k in range(d):
                Mfull[v][k] = T.add(Mfull[v][k], T.mul(delta[j], Ms[j][v][k]))
    badF = 0
    for a in range(n):
        for v in range(a, n):
            if not T.iszero(resid(S_full, Mfull, a, v)):
                badF += 1
    nonzero = sum(1 for (a, b) in cls_of if not T.iszero(S_full(a, b)))
    zero_on_contacts = all(T.iszero(S_full(a, b)) for a, b in E)
    log.update(flex_check_failures=badF, flex_S_nonzero_pairs=nonzero, flex_S_zero_on_contacts=zero_on_contacts,
               flex_delta_float=[T.to_float(x) for x in delta])
    if badF or nonzero == 0 or not zero_on_contacts:
        fail(log, out, "flex check failed")
    log["seconds_flex"] = round(time.time() - t0, 1)
    # ------------------------------------------------------------------ (S) the stress (as jam116.py)
    C = {}
    for (a, b), j in zip(E, et):
        for (x, y) in ((a, b), (b, a)):
            for k in range(d):
                term = T.sub(W[k][y], T.mul(c, W[k][x]))
                key = (x, k)
                if key not in C:
                    C[key] = [ZERO] * ntype
                C[key][j] = T.add(C[key][j], term)
    keys = sorted(C)
    comb = []
    for r in range(2 * ntype + 4):
        coef = rng.integers(-1, 2, size=len(keys))
        row = [ZERO] * ntype
        for kk, cf_ in zip(keys, coef):
            if cf_ == 0:
                continue
            for j in range(ntype):
                if not T.iszero(C[kk][j]):
                    row[j] = T.add(row[j], C[kk][j]) if cf_ > 0 else T.sub(row[j], C[kk][j])
        comb.append(row)
    piv_cols, r_ = [], 0
    for col in range(ntype):
        pr = next((i for i in range(r_, len(comb)) if not T.iszero(comb[i][col])), None)
        if pr is None:
            continue
        comb[r_], comb[pr] = comb[pr], comb[r_]
        inv = T.inv(comb[r_][col])
        comb[r_] = [T.mul(x, inv) for x in comb[r_]]
        for i in range(len(comb)):
            if i != r_ and not T.iszero(comb[i][col]):
                fct = comb[i][col]
                comb[i] = [T.sub(x, T.mul(fct, y)) for x, y in zip(comb[i], comb[r_])]
        piv_cols.append(col)
        r_ += 1
    freeS = [col for col in range(ntype) if col not in piv_cols]
    basisS = []
    for fcol_ in freeS:
        v = [ZERO] * ntype
        v[fcol_] = ONE
        for i, pcol in enumerate(piv_cols):
            v[pcol] = T.neg(comb[i][fcol_])
        basisS.append(v)
    Bf = np.array([[T.to_float(x) for x in v] for v in basisS]).T
    z = np.linalg.lstsq(Bf, wlp_types, rcond=None)[0]
    zq = [fq(int(round(x * 2 ** 40)), 2 ** 40) for x in z]
    omega_t = [ZERO] * ntype
    for zz, v in zip(zq, basisS):
        for j in range(ntype):
            omega_t[j] = T.add(omega_t[j], T.mul(T.const(zz), v[j]))
    pos = all(T.sign(x) > 0 for x in omega_t)
    badS = 0
    for kk in keys:
        acc = ZERO
        for j in range(ntype):
            if not T.iszero(C[kk][j]):
                acc = T.add(acc, T.mul(omega_t[j], C[kk][j]))
        if not T.iszero(acc):
            badS += 1
    tot = ZERO
    for j in range(ntype):
        tot = T.add(tot, T.mul(T.const(int(cnt[j])), omega_t[j]))
    wmin = min(T.to_float(x) for x in omega_t) / T.to_float(tot)
    log.update(stress_types=ntype, stress_space_dim_on_types=len(basisS), stress_positive=pos,
               stress_equilibrium_failures=badS, stress_min_relative=wmin * len(E),
               stress_ratio_QQ_over_PP=T.to_float(omega_t[tQQ]) / T.to_float(omega_t[tPP]))
    if not pos or badS:
        fail(log, out, "stress check failed")
    log["seconds_stress"] = round(time.time() - t0, 1)
    # ------------------------------------------------------------------ (H) the second-order form (as jam116.py)
    Y = [[None] * d for _ in range(n)]
    for a in range(n):
        for k in range(d):
            acc = ZERO
            for l in range(d):
                if not T.iszero(Mfull[a][l]):
                    acc = T.add(acc, T.mul(Mfull[a][l], GBBinv[l][k]))
            Y[a][k] = acc

    def FF(a, b):
        acc = ZERO
        for k in range(d):
            if not T.iszero(Y[a][k]) and not T.iszero(Mfull[b][k]):
                acc = T.add(acc, T.mul(Y[a][k], Mfull[b][k]))
        return acc
    diag = [FF(a, a) for a in range(n)]
    H = ZERO
    for (a, b), j in zip(E, et):
        term = T.sub(T.mul(T.const(2), FF(a, b)), T.mul(c, T.add(diag[a], diag[b])))
        H = T.add(H, T.mul(omega_t[j], term))
    flexnorm2 = ZERO
    for a in range(n):
        flexnorm2 = T.add(flexnorm2, diag[a])
    Hs = T.sign(H)
    log.update(H_sign=Hs, H_float=T.to_float(H), flex_norm2_float=T.to_float(flexnorm2),
               H_over_norm2_total_stress=T.to_float(H) / T.to_float(flexnorm2) / T.to_float(tot))
    log["seconds_H"] = round(time.time() - t0, 1)
    if Hs <= 0:
        fail(log, out, "H not positive")
    log["status"] = "JAMMED, NOT INFINITESIMALLY JAMMED (prestress stable: strict local minimum modulo rotations)"
    log["seconds"] = round(time.time() - t0, 1)
    json.dump(log, open(out, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in log.items() if k not in ("code_sha256", "LA", "new_values_tower")}, default=str))


if __name__ == "__main__":
    main()
