"""Census version 9, item A2: second-order (prestress-stability) certificates for exact codes whose contact framework
is not bar-rigid but carries an exact strictly positive equilibrium stress (Theorem 2.9 of paper_v9).

  nice -n 10 python3 -I prestress9.py OUTDIR [--only d_n,...] [--max-mb 1500]

The entries are those of paper_v8/above (run_above.py) with an exact strictly positive stress whose rank bound fell
short: (11,604), (22,166), (30,656); --only restricts them. For each:
 (C) TRUSTED, certify/v5/core.py (released, unchanged): the exact core over K = Q or Q(sqrt D) (the field recorded by
     run_above.py): the Gram matrix G recognised and checked exactly (unit diagonal, rank d: G_BB positive definite and
     G = G_VB G_BB^-1 G_BV), the contacts E exactly at c, every other core pair and every stored rattler below c.
 (S) TRUSTED: a stress omega constant on classes of contacts (run_above.class_stress, untrusted coefficients), checked
     exactly: equilibrium at every core point (v5/rigidity.py check_stress) and every entry positive (exact signs).
 (F) the infinitesimal flexes. UNTRUSTED: K_p = M^T M mod p for the rigidity matrix M of Proposition 3.19(d) (FLINT
     nmod_mat, primes below 2^26 in which D is a square), its nullspace in reduced echelon form, a set J of k pivot
     rows complementing the g infinitesimal rotations, and the rows of J lifted to K by the Chinese remainder theorem
     and rational reconstruction (two embeddings sqrt D -> +-r per prime when D > 1), primes being added until the
     lift is stable. TRUSTED: every lifted vector f is an exact flex (M f = 0 over K, v5/rigidity.py apply_R); the g
     rotations (exact, rigidity.rotations) and the k lifted vectors are independent (rank g + k of their images mod a
     prime, a lower bound for the rank over K); and the rank of K_p over F_p (FLINT) is |V| d - g - k, so the flexes
     form a space of dimension at most g + k over K (Proposition 3.19(c)). Hence the flexes are exactly the rotations
     plus the span of the lifted vectors.
 (H) TRUSTED: the matrix H_omega(f_i, f_j) = sum_E omega_ab (<f_ia, f_jb> + <f_ib, f_ja> - c <f_ia, f_ja>
     - c <f_ib, f_jb>) over K (v5/exact2.py H_matrix) and its exact LDL^T pivots (kfield.is_pos_def).
If every pivot is positive, Theorem 2.9 (second-order criterion) makes the exact core, with the stored rattlers, a
strict local minimum of the largest cosine modulo rotations and rattlers (prestress stable; quadratic growth along the
flexes). The output records the pivots, the heights of the lifted vectors and the timings.
"""
import hashlib, json, math, os, sys, time, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
sys.path.insert(0, os.path.join(CERT, "paper_v8", "above"))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, os.path.join(CERT, "v5"))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import flint
import run_above as RA
import core
import rigidity as RG
import exact2 as X2
from kfield import K, fq, prime_for, is_pos_def

SHA = {"prestress9.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
       "run_above.py": RA.SHA}
ENTRIES = [(22, 166), (11, 604), (30, 656)]


# ============================================================================================ UNTRUSTED: lifting
def ratrec(a, m):
    """Rational reconstruction: p/q = a mod m with |p|, q <= sqrt(m/2), or None."""
    a %= m
    bound = math.isqrt(m // 2)
    r0, r1, s0, s1 = m, a, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    if s1 < 0:
        r1, s1 = -r1, -s1
    if math.gcd(r1, s1) != 1:
        return None
    return fq(r1, s1)


def kernel_rref_modp(G, B, E, d, nV, p, rt, max_mb):
    """K_p = M^T M mod p; its rank (FLINT) and its nullspace in reduced echelon form (rows) with the pivot columns."""
    n = nV * d
    if n * n * 8 / 2 ** 20 > max_mb:
        raise MemoryError("dense K_p above the memory limit")
    M = RA.modp_csr(G, B, E, d, nV, p, rt)
    assert M.data.min() >= 0 and M.data.max() < p < 2 ** 26
    assert int(np.diff(M.tocsc().indptr).max()) <= 2 ** 11      # sums of at most 2^11 products below 2^52
    Kc = (M.T.tocsr() @ M.tocsc()).tocoo()
    del M
    A = flint.nmod_mat(n, n, p)
    vals = np.mod(Kc.data, p)
    for i, j, v in zip(Kc.row.tolist(), Kc.col.tolist(), vals.tolist()):
        if v:
            A[i, j] = v
    del Kc, vals
    X, nul = A.nullspace()
    del A
    Bm = flint.nmod_mat(nul, n, p)
    for c in range(nul):
        for i in range(n):
            x = int(X[i, c])
            if x:
                Bm[c, i] = x
    del X
    R, rk = Bm.rref()
    rows = [[int(R[c, i]) for i in range(n)] for c in range(rk)]
    piv = [next(i for i in range(n) if row[i]) for row in rows]
    return n - nul, rows, piv


def lift(G, B, E, d, nV, D, rot, max_mb, log, max_primes=40):
    """The pivot rows J of the echelon nullspace that complement the rotations, lifted to K (untrusted)."""
    n = nV * d
    g = len(rot)
    start = 1 << 26
    acc = None                     # (modulus, {pivot: (list a-parts, list b-parts)})
    J = None
    prev = None
    ranks = []
    for it in range(max_primes):
        p, rt = prime_for(D, start=start)
        start = p
        try:
            emb = []
            for r_ in ([rt] if D == 1 else [rt, p - rt]):
                rk, rows, piv = kernel_rref_modp(G, B, E, d, nV, p, r_, max_mb)
                emb.append((rk, rows, piv))
        except ValueError:
            continue                                     # a denominator divisible by p
        if any(e[2] != emb[0][2] for e in emb):
            continue
        rk, rows, piv = emb[0]
        ranks.append({"prime": p, "rank": rk})
        if J is None:
            # J: pivots of the nullspace not covered by the rotations (rref of the rotations restricted to the pivots)
            Rm = flint.nmod_mat(g, len(piv), p)
            for a, u in enumerate(rot):
                for c, col in enumerate(piv):
                    x = u[col]
                    if not x.iszero():
                        Rm[a, c] = x.modp(p, rt)
            Rr, rr = Rm.rref()
            if rr != g:
                raise ValueError("rotations not independent on the pivot columns")
            used = set(next(c for c in range(len(piv)) if int(Rr[a, c])) for a in range(g))
            J = [piv[c] for c in range(len(piv)) if c not in used]
            log["pivots_J"] = J
            pivset = list(piv)
        elif piv != pivset:
            continue
        rowof = {pv: i for i, pv in enumerate(piv)}
        inv2 = pow(2, -1, p)
        cur = {}
        for pv in J:
            if D == 1:
                a_part = emb[0][1][rowof[pv]]
                b_part = [0] * n
            else:
                x1, x2 = emb[0][1][rowof[pv]], emb[1][1][rowof[pv]]
                i2r = pow(2 * rt % p, -1, p)
                a_part = [(u + v) * inv2 % p for u, v in zip(x1, x2)]
                b_part = [(u - v) * i2r % p for u, v in zip(x1, x2)]
            cur[pv] = (a_part, b_part)
        if acc is None:
            acc = (p, cur)
        else:
            m0, old = acc
            m1 = m0 * p
            inv = pow(m0, -1, p)
            new = {}
            for pv in J:
                na, nb = [], []
                for (oa, ob), (ca, cb) in [((old[pv][0][i], old[pv][1][i]), (cur[pv][0][i], cur[pv][1][i])) for i in range(n)]:
                    na.append(oa + m0 * ((ca - oa) * inv % p))
                    nb.append(ob + m0 * ((cb - ob) * inv % p))
                new[pv] = (na, nb)
            acc = (m1, new)
        m, vals = acc
        rec = {}
        ok = True
        for pv in J:
            va, vb = [], []
            for x, y in zip(*vals[pv]):
                ra_ = ratrec(x, m) if x else fq(0)
                rb_ = ratrec(y, m) if y else fq(0)
                if ra_ is None or rb_ is None:
                    ok = False
                    break
                va.append(ra_)
                vb.append(rb_)
            if not ok:
                break
            rec[pv] = [K(a_, b_) for a_, b_ in zip(va, vb)]
        log["primes_used"] = len(ranks)
        log["modulus_bits"] = m.bit_length()
        if ok and prev is not None and all(rec[pv] == prev[pv] for pv in J):
            log["rank_trials"] = ranks
            return [rec[pv] for pv in J], ranks
        prev = rec if ok else None
    raise ValueError("lift did not stabilise")


# ============================================================================================ one entry
def one(d, n, rec8, max_mb, log):
    disc = int(rec8.get("disc", 1))
    K.D = disc
    t0 = time.time()
    C = core.build(SNAP, d, n, disc)
    if C["status"] != "ok":
        log["status"] = "core: " + C["status"]
        return
    log.update(C["log"])
    G, V, E, c, r, B, W, GBB, GBBinv = (C[k] for k in ("G", "V", "E", "c", "r", "B", "W", "GBB", "GBBinv"))
    nV = len(V)
    g = d * (d - 1) // 2
    log.update(field_disc=disc, core_rank_ok=(r == d), seconds_core=round(time.time() - t0, 1))
    if r != d:
        log["status"] = "core does not span"
        return
    # (S) the stress
    om = None
    unif = [K(1)] * len(E)
    if RG.check_stress(W, E, d, nV, unif):
        om, how = unif, "uniform"
    elif RA.class_stress(W, E, d, nV, log, RA.neighbour_classes(E, nV), "neighbours") is not None:
        om, how = log.pop("_class_stress"), "constant on classes of contacts (common-neighbour counts)"
    if om is None:
        log["status"] = "no exact strictly positive stress"
        return
    ok_eq = RG.check_stress(W, E, d, nV, om)
    ok_pos = all(x.sign() > 0 for x in om)
    tot = om[0]
    for x in om[1:]:
        tot = tot + x
    log.update(stress=how, stress_equilibrium_exact=ok_eq, stress_positive_exact=ok_pos,
               stress_min_relative=float(min(om, key=lambda x: float(x))) * len(E) / float(tot))
    if not (ok_eq and ok_pos):
        log["status"] = "stress check failed"
        return
    # (F) the flexes
    rot = RG.rotations(nV, W, GBBinv, d)
    for u in rot:
        if any(not x.iszero() for x in RG.apply_R(G, B, E, u, d, nV)):
            log["status"] = "a rotation is not a flex (impossible)"
            return
    t1 = time.time()
    flex, ranks = lift(G, B, E, d, nV, disc, rot, max_mb, log)
    log["seconds_lift"] = round(time.time() - t1, 1)
    k = len(flex)
    exact = [all(x.iszero() for x in RG.apply_R(G, B, E, f, d, nV)) for f in flex]
    log["lifted_flexes_exact"] = exact
    if not all(exact):
        log["status"] = "a lifted vector is not an exact flex"
        return
    rk_lo = ranks[-1]["rank"]
    pchk, rtchk = prime_for(disc, start=ranks[-1]["prime"])
    Mi = flint.nmod_mat(g + k, nV * d, pchk)
    for a, u in enumerate(rot + flex):
        for i, x in enumerate(u):
            if not x.iszero():
                Mi[a, i] = x.modp(pchk, rtchk)
    ind = Mi.rank()
    log.update(flexes_beyond_rotations=k, independent_rank=int(ind), rank_K_lower=int(rk_lo),
               flex_space_exact=bool(ind == g + k and nV * d - rk_lo == g + k))
    if not log["flex_space_exact"]:
        log["status"] = "flex space not determined exactly"
        return
    hts = [max(max(x.a.p.bit_length() if False else int(abs(x.a.p)).bit_length(), int(x.a.q).bit_length(),
                   int(abs(x.b.p)).bit_length(), int(x.b.q).bit_length()) for x in f) for f in flex]
    log["lift_height_bits"] = hts
    # (H)
    t2 = time.time()
    H = X2.H_matrix(GBB, list(range(nV)), E, c, om, flex, d)
    pd, piv = is_pos_def(H)
    log.update(H=[[repr(x) for x in row] for row in H], H_float=[[float(x) for x in row] for row in H],
               H_pivots=[repr(x) for x in piv], H_pos_def=bool(pd), seconds_H=round(time.time() - t2, 1))
    log["status"] = ("STRICT (prestress stable: positive stress, H positive definite on the flexes modulo rotations; "
                     "Theorem 2.9)") if pd else "undecided: H not positive definite for this stress"


def main():
    out = sys.argv[1]
    args = sys.argv[2:]
    opt = lambda k_, dflt: args[args.index(k_) + 1] if k_ in args else dflt
    only = set(opt("--only", "").split(",")) if "--only" in args else None
    max_mb = int(opt("--max-mb", "1500"))
    os.makedirs(out, exist_ok=True)
    for d, n in ENTRIES:
        if only and f"{d}_{n}" not in only:
            continue
        js = os.path.join(out, f"d{d}_n{n}.json")
        if os.path.exists(js):
            continue
        rec8 = json.load(open(os.path.join(CERT, "paper_v8", "above", "out", f"d{d}_n{n}.json")))
        log = {"d": d, "n": n, "program_sha256": SHA}
        t0, c0 = time.time(), time.process_time()
        try:
            one(d, n, rec8, max_mb, log)
        except MemoryError as e:
            log["status"] = "skipped: " + str(e)
        except Exception:
            log["status"] = "error"
            log["error"] = traceback.format_exc()[-1500:]
        log["wall"] = round(time.time() - t0, 1)
        log["cpu_seconds"] = round(time.process_time() - c0, 1)
        with open(js + ".part", "w") as f:
            json.dump(log, f, default=str)
        os.replace(js + ".part", js)
        print(d, n, log["status"], log.get("H_float"), log["wall"], flush=True)


if __name__ == "__main__":
    main()
