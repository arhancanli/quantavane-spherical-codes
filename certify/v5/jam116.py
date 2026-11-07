"""A jammed spherical code that is not infinitesimally jammed: the exact certificate for the entry (22,116) of the
archived table (paper, version 5, Theorem 5.4).

  python3 -I certify/v5/jam116.py SNAPDIR RECOGNITION.json OUT.json

RECOGNITION.json is written by recognise116.py (UNTRUSTED): the classes of pairs, the tower of quadratic fields and the
class values in it, the basis B, the sign pattern of a numerical flex, the edge types and a numerical stress.
Everything below is TRUSTED and exact (arithmetic in the field L = Q(sqrt145, s1, s2) of degree 8, tower.py; signs at
the fixed real embedding in FLINT ball arithmetic; ranks over a prime field), except the distance bound (D), which is
in ball arithmetic with rigorous radii.

 (L) [L:Q] = 8 and the embedding is real (tower.check_degree).
 (G) G = I + sum_k t_k A_k (A_k the 0/1 matrix of class k) is the Gram matrix of a configuration X* of 116 unit vectors
     spanning R^22: G_BB is positive definite (exact LDL^T pivots, positive) and G = G_{.B} G_BB^{-1} G_{B.} exactly.
 (C) c = t_12 satisfies 30c^2 + 5c - 1 = 0; every other class value is below c; every point lies in at least 22 contacts.
     So M(X*) = c, the contacts E are the pairs of class 12, and there are no rattlers.
 (R) the bar rigidity matrix of (X*, E) (coordinates in the basis x_b, b in B) has rank at least 116*22 - 231 - 1 over a
     prime field (a lower bound for its rank over L): modulo rotations, the infinitesimal flexes form a space of
     dimension at most one.
 (F) a nonzero symmetric matrix S over L, zero on the diagonal and on the contacts, and a matrix M with
     S = W^T M^T + M W (W = G_BB^{-1} G_{B.}) checked entry by entry: F = M G_BB^{-1} X_B is an infinitesimal flex with
     first-order Gram change S != 0, hence not an infinitesimal rotation. So X* is not infinitesimally jammed.
 (S) an equilibrium stress omega, constant on the given edge types, checked exactly at every point, with every entry
     positive.
 (H) H_omega(F) = sum_E omega_ab (2<f_a,f_b> - c|f_a|^2 - c|f_b|^2) > 0, with <f_a,f_b> = (M G_BB^{-1} M^T)_ab exactly.
 By the paper's second-order theorem, (G), (C), (R), (F), (S) and (H) make X* a strict local minimum of the largest
 cosine modulo rotations (with quadratic but not linear growth); in particular X* is jammed.
 (D) the distance between X* (suitably rotated) and the exactly normalised stored code.
"""
import hashlib, json, math, os, sys, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/flintlib"))
import numpy as np
import flint
from kfield import K, fq, prime_for
from tower import Tower
import exact2 as X2
import dist

SHA = {}
for _f in ("jam116.py", "tower.py", "kfield.py", "exact2.py", "dist.py"):
    SHA[_f] = hashlib.sha256(open(os.path.join(HERE, _f), "rb").read()).hexdigest()


def fail(log, out, msg):
    log["status"] = "FAILED: " + msg
    json.dump(log, open(out, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in log.items() if k != "code_sha256"}, default=str))
    raise SystemExit(msg)


def main():
    snapdir, recf, out = sys.argv[1], sys.argv[2], sys.argv[3]
    t0 = time.time()
    rec = json.load(open(recf))
    d, n = rec["d"], rec["n"]
    g = d * (d - 1) // 2
    log = {"d": d, "n": n, "program": "jam116.py", "code_sha256": SHA,
           "recognition_sha256": hashlib.sha256(open(recf, "rb").read()).hexdigest()}
    # ------------------------------------------------------------------ (L) the field
    def k0(r):   # Delta from a PSLQ relation [a, b, c]: a Delta + b + c sqrt145 = 0
        a, b, c = r
        return K(fq(-b, a), fq(-c, a))
    K.D = 145
    T = Tower(145, k0(rec["Delta1_rel"]), k0(rec["Delta2_rel"]), rec["sigma1"], rec["sigma2"])
    facts = T.check_degree()
    log["field"] = {k: (v if isinstance(v, bool) else str(v)) for k, v in facts.items()}
    log["field"]["Delta1"] = repr(T.Delta1)
    log["field"]["Delta2"] = repr(T.Delta2)
    if not (facts["degree_8"] and facts["Delta1_positive"] and facts["Delta2_positive"]):
        fail(log, out, "field check failed")

    def from_rel(r):
        r0 = r[0]
        co = [fq(-x, r0) for x in r[1:]]
        return ((K(co[0], co[1]), K(co[2], co[3])), (K(co[4], co[5]), K(co[6], co[7])))
    tvals = [from_rel(r) for r in rec["class_rel"]]
    nc = len(tvals)
    ONE, ZERO = T.const(1), T.const(0)
    # untrusted sanity: the recognised values agree with the double-precision classes
    log["class_values_max_dev"] = max(abs(T.to_float(x) - v) for x, v in zip(tvals, rec["classes"]))
    # ------------------------------------------------------------------ (C) the cosine
    ci = max(range(nc), key=lambda k: rec["classes"][k])
    c = tvals[ci]
    poly = T.add(T.add(T.mul(T.const(30), T.mul(c, c)), T.mul(T.const(5), c)), T.const(-1))
    below = [k for k in range(nc) if k != ci and T.lt(tvals[k], c)]
    log.update(cos_class=ci, cos_value_float=T.to_float(c), cos_satisfies_30x2_5x_1=T.iszero(poly),
               classes_below_c=len(below), classes=nc)
    if not (T.iszero(poly) and len(below) == nc - 1):
        fail(log, out, "cosine check failed")
    pc = rec["pair_class"]
    cls_of = {}
    idx = 0
    for a in range(n):
        for b in range(a + 1, n):
            cls_of[(a, b)] = pc[idx]
            idx += 1

    def Gv(a, b):
        if a == b:
            return ONE
        return tvals[cls_of[(a, b)] if a < b else cls_of[(b, a)]]
    E = [(a, b) for (a, b), k in cls_of.items() if k == ci]
    deg = [0] * n
    for a, b in E:
        deg[a] += 1
        deg[b] += 1
    log.update(contacts=len(E), min_contacts_per_point=min(deg), rattlers=sum(1 for x in deg if x == 0))
    if min(deg) < d:
        fail(log, out, "a point has fewer than d contacts")
    # ------------------------------------------------------------------ (G) the Gram matrix
    B = rec["basis"]
    GBB = [[Gv(a, b) for b in B] for a in B]
    # LDL^T pivots and the inverse (Gauss-Jordan), exactly
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
    # Gauss-Jordan without row exchanges: the pivots are those of LDL^T (leading principal minors ratios)
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
    # ------------------------------------------------------------------ (R) rank over a prime field
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
    # ------------------------------------------------------------------ (F) the flex
    fcls = rec["flex_classes"]
    sgn = {tuple(map(int, k.split(","))): v for k, v in rec["flex_sign"].items()}

    def S_basis(j, a, b):
        """S for delta = e_j (class fcls[j]): sign * [class(a,b) == fcls[j]]."""
        if a == b:
            return 0
        key = (a, b) if a < b else (b, a)
        if cls_of[key] != fcls[j]:
            return 0
        return sgn.get(key, 0)

    def M_of(Sfun):
        # M_vk = S(b_k, v) - 1/2 sum_l S(b_k, b_l) W_lv
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
    # untrusted choice of sample pairs for solving delta; exact solve; exact check on all pairs below
    rng = np.random.default_rng(116)
    nonB = [v for v in range(n) if v not in set(B)]   # pairs with a point of B give identically zero equations
    Ms = [M_of(lambda a, b, j=j: S_basis(j, a, b)) for j in range(len(fcls))]
    small = [v for v in nonB if deg[v] == min(deg)]     # the points of the smallest orbits meet few pairs: include them
    order = [(v, v) for v in nonB] + [(s, v) for s in small for v in nonB if v != s]   # diagonal: tangency
    rest = [(a, b) for i, a in enumerate(nonB) for b in nonB[i + 1:] if (a, b) not in set(order)]
    rng.shuffle(rest)
    order += rest
    rowsA, used = [], 0

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
        # keep only the reduced nonzero rows (at most len(fcls)) for the next round
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
    firstbad = []
    for a in range(n):
        for v in range(a, n):
            r_ = resid(S_full, Mfull, a, v)
            if not T.iszero(r_):
                badF += 1
                if len(firstbad) < 5:
                    firstbad.append((a, v, T.to_float(r_)))
    log["flex_first_failures"] = firstbad
    nonzero = sum(1 for (a, b) in cls_of if not T.iszero(S_full(a, b)))
    zero_on_contacts = all(T.iszero(S_full(a, b)) for a, b in E)
    log.update(flex_check_failures=badF, flex_S_nonzero_pairs=nonzero, flex_S_zero_on_contacts=zero_on_contacts,
               flex_delta_float=[T.to_float(x) for x in delta])
    if badF or nonzero == 0 or not zero_on_contacts:
        fail(log, out, "flex check failed")
    log["seconds_flex"] = round(time.time() - t0, 1)
    # ------------------------------------------------------------------ (S) the stress
    et = rec["edge_type"]
    ntype = len(rec["edge_types"])
    # coefficient C[(a,k)][j] = sum over edges ab of type j of (W_kb - c W_ka)
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
    # untrusted: a few random combinations of the equations, null space over L
    comb = []
    for r in range(2 * ntype + 4):
        coef = rng.integers(-1, 2, size=len(keys))
        row = [ZERO] * ntype
        for kk, cf in zip(keys, coef):
            if cf == 0:
                continue
            for j in range(ntype):
                if not T.iszero(C[kk][j]):
                    row[j] = T.add(row[j], C[kk][j]) if cf > 0 else T.sub(row[j], C[kk][j])
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
    # choose the combination closest (least squares, floating point, untrusted) to the LP stress, rounded to rationals
    Bf = np.array([[T.to_float(x) for x in v] for v in basisS]).T
    wlp = np.array(rec["lp_stress_by_type"])
    z = np.linalg.lstsq(Bf, wlp, rcond=None)[0]
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
    cnt = np.bincount(et, minlength=ntype)
    for j in range(ntype):
        tot = T.add(tot, T.mul(T.const(int(cnt[j])), omega_t[j]))
    wmin = min(T.to_float(x) for x in omega_t) / T.to_float(tot)
    log.update(stress_types=ntype, stress_space_dim_on_types=len(basisS), stress_positive=pos,
               stress_equilibrium_failures=badS, stress_min_relative=wmin * len(E))
    if not pos or badS:
        fail(log, out, "stress check failed")
    log["seconds_stress"] = round(time.time() - t0, 1)
    # ------------------------------------------------------------------ (H) the second-order form on the flex
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
    # ------------------------------------------------------------------ (D) the distance to the stored code
    entry, dec, Xf = X2.read_code(snapdir, d, n)
    old = flint.ctx.prec
    flint.ctx.prec = 300
    Ga = [[T.arb(Gv(a, b), 300) for b in B] for a in B]
    L = [[flint.arb(0)] * d for _ in range(d)]
    for j in range(d):
        s = Ga[j][j] - sum((L[j][k] * L[j][k] for k in range(j)), flint.arb(0))
        if not s > 0:
            fail(log, out, "Cholesky pivot not provably positive")
        L[j][j] = s.sqrt()
        for i in range(j + 1, d):
            L[i][j] = (Ga[i][j] - sum((L[i][k] * L[j][k] for k in range(j)), flint.arb(0))) / L[j][j]
    rowsP = []
    for a in range(n):
        gvec = [T.arb(Gv(a, b), 300) for b in B]
        y = []
        for i in range(d):
            y.append((gvec[i] - sum((L[i][k] * y[k] for k in range(i)), flint.arb(0))) / L[i][i])
        rowsP.append(y)
    P = np.array([[float(x.mid()) for x in r] for r in rowsP])
    Xn = Xf / np.linalg.norm(Xf, axis=1)[:, None]
    U_, _, Vt_ = np.linalg.svd(P.T @ Xn)
    Q = dist.rational_orthogonal(U_ @ Vt_, d)
    orth = Q.transpose() * Q == flint.fmpq_mat(d, d, [int(i == j) for i in range(d) for j in range(d)])
    Qa = [[flint.arb(Q[j, k]) for k in range(d)] for j in range(d)]
    tot2 = flint.arb(0)
    mx = 0.0
    for a in range(n):
        q = [Fraction(t_) for t_ in dec[a]]
        s2 = sum(x * x for x in q)
        nrm = flint.arb(flint.fmpq(s2.numerator, s2.denominator)).sqrt()
        for k in range(d):
            z_ = sum((rowsP[a][j] * Qa[j][k] for j in range(d)), flint.arb(0))
            df = z_ - flint.arb(flint.fmpq(q[k].numerator, q[k].denominator)) / nrm
            tot2 = tot2 + df * df
            mx = max(mx, math.nextafter(float(abs(df).upper()), math.inf))
    flint.ctx.prec = old
    log.update(rotation_exactly_orthogonal=bool(orth), dist_ls_upper=math.nextafter(float(tot2.sqrt().upper()), math.inf),
               dist_coord_upper=mx)
    ok = (Hs > 0)
    log["status"] = ("JAMMED, NOT INFINITESIMALLY JAMMED (prestress stable: strict local minimum modulo rotations)"
                     if ok else "H not positive: undecided")
    log["seconds"] = round(time.time() - t0, 1)
    json.dump(log, open(out, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in log.items() if k not in ("code_sha256",)}, default=str))


if __name__ == "__main__":
    main()
