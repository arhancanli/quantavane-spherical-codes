"""Census version 9, item A1: strictness of hyperstatic entries from an exact configuration recognised at higher degree
(Theorem 3.14 with Proposition 3.15 of the paper), decided in exact arithmetic. The trusted part is self-contained:
no function of the census is replaced at run time; every exact test is called explicitly.

  nice -n 10 python3 -I strict9.py OUTDIR [--only d_n,...] [--bits 2000] [--maxdeg 30] [--data DATADIR] [--class local]

UNTRUSTED (recog9.survey): the census certificate re-run; P' refined to 2^-BITS; the minimal polynomial p of a
primitive element theta (an integer polynomial of degree <= MAXDEG found by LLL); the common cosine and every distinct
inner product of P' as elements of Q(theta) = Q[x]/(p); the choice of B, of E' and of the coefficients lambda_f; the
approximate inverse; the rattler positions.

TRUSTED:
 (F)    p is irreducible over Q (FLINT); theta* is the unique real root of p (FLINT's isolation, PREC bits) within 2^-200
        of the refined value; c* = c(theta*) for the recognised element c of Q[x]/(p).
 (E)    the recognised symmetric matrix G over Q[x]/(p), unit diagonal: G_ab = c exactly on every contact; for every
        other pair of V, G_ab != c and c* - G_ab(theta*) > 0 in Arb; G_BB positive definite (ball Cholesky with
        positive pivots); G = G_VB G_BB^-1 G_BV exactly in Q[x]/(p). By Proposition 3.15 G is the Gram matrix of a
        configuration P* in (S^{d-1})^V spanning R^d, with E(P*) = E and M(P*) = c*.
 (S')   on the box Pi around the Cholesky-gauge enclosure of P*: DF_E'(p) nonsingular and the combined stress of E
        strictly positive for every p in Pi (exactcheck.square_test, mode 'box', exact rationals; Lemma 3.4).
 (span) d points of the box: ||B^-1||_inf (max row sum of the radii) < 1, B^-1 exact.
 (G')   every rattler, placed by untrusted code, below c* against every point (strict.rattler_gaps, exact rationals).
 (D)    the distance to the exactly normalised stored code (strict.distance_to_stored, Arb).
Theorem 3.14 then makes P* with the rattlers a strict local minimum modulo rotations and rattlers, with value c*.
With --data, the exact description (field, element of every distinct inner product, index matrix) is written to
DATADIR/d{d}_n{n}.json.gz.
"""
import gzip, hashlib, json, math, os, sys, time, traceback
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rec9lib as L                     # library paths
import recog9 as R9
CERT = L.CERT
sys.path.insert(0, os.path.join(CERT, "paper_v8", "verify"))
import exactcheck as XC                 # the exact tests of record (paper_v8)
import numpy as np
import flint
import cq
import strict                           # census/strict.py: rattler_gaps, distance_to_stored (both trusted)
import sphcert as S
import snapshot

PREC = 1024
BLK = 64
SHA = dict(R9.SHA)
SHA.update({"strict9.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
            "exactcheck.py": XC.SHA, "strict.py": strict.SHA["strict.py"]})


# ============================================================================================ TRUSTED: the field
class Field9:
    """Q[x]/(p) for an irreducible integer polynomial p, with the real embedding x -> theta*, the unique real root of p
    whose isolating ball lies within 2^-200 of the given rational approximation."""

    def __init__(self, coeffs, theta):
        self.p = flint.fmpq_poly([int(c) for c in coeffs])
        fac = self.p.factor()[1]
        if len(fac) != 1 or fac[0][1] != 1:
            raise ValueError("polynomial not irreducible")
        self.D = self.p.degree()
        num, den_bits = theta
        old = flint.ctx.prec
        flint.ctx.prec = PREC
        th = flint.arb(flint.fmpq(int(num), 1 << int(den_bits)))
        tol = flint.arb(flint.fmpq(1, 1 << 200))
        roots = self.p.numer().complex_roots()
        real = [r.real for r, mlt in roots if r.imag == 0 and bool(abs(r.real - th) < tol)]
        flint.ctx.prec = old
        if len(real) != 1:
            raise ValueError("no unique real root of p near the refined value")
        self.theta = real[0]
        self.theta_rad = float(self.theta.rad())

    def elem(self, co):
        return flint.fmpq_poly([flint.fmpq(int(a), int(b)) for a, b in co]) % self.p

    def value(self, a):
        old = flint.ctx.prec
        flint.ctx.prec = PREC
        acc = flint.arb(0)
        for c in reversed(a.coeffs()):
            acc = acc * self.theta + flint.arb(c)
        flint.ctx.prec = old
        return acc

    def mul(self, a, b):
        return (a * b) % self.p

    def inv(self, a):
        g, s, _ = a.xgcd(self.p)
        if g.degree() != 0:
            raise ZeroDivisionError
        return (s / g) % self.p


def coeff_mats(F, rows, ncols):
    """TRUSTED. A matrix of field elements (list of rows) as D rational matrices of the power basis."""
    m = len(rows)
    ent = [[flint.fmpq(0)] * (m * ncols) for _ in range(F.D)]
    for i, row in enumerate(rows):
        for j, el in enumerate(row):
            for t, q in enumerate(el.coeffs()):
                ent[t][i * ncols + j] = q
    return [flint.fmpq_mat(m, ncols, e) for e in ent]


def rank_d_exact(F, els, cid, B, nV, d):
    """TRUSTED. G = G_VB G_BB^-1 G_BV exactly over Q[x]/(p), for the matrix G with entries els[cid[i][j]]: the inverse of
    G_BB by Gauss-Jordan over the field, W = G_BB^-1 G_BV, and every row block of G_VB W, computed through the power
    basis (products of rational matrices, x^t reduced by p), compared with G. Returns (ok, message)."""
    D = F.D
    zero = F.elem([(0, 1)])
    one = F.elem([(1, 1)])
    M_ = [[els[cid[i][j]] for j in B] for i in B]
    Ainv = [[one if i == j else zero for j in range(d)] for i in range(d)]
    for c_ in range(d):
        pv = next((r for r in range(c_, d) if M_[r][c_] != 0), None)
        if pv is None:
            return False, "G_BB singular"
        M_[c_], M_[pv] = M_[pv], M_[c_]
        Ainv[c_], Ainv[pv] = Ainv[pv], Ainv[c_]
        iv = F.inv(M_[c_][c_])
        M_[c_] = [F.mul(iv, e) for e in M_[c_]]
        Ainv[c_] = [F.mul(iv, e) for e in Ainv[c_]]
        for r in range(d):
            if r != c_ and M_[r][c_] != 0:
                f = M_[r][c_]
                M_[r] = [(e - F.mul(f, g)) % F.p for e, g in zip(M_[r], M_[c_])]
                Ainv[r] = [(e - F.mul(f, g)) % F.p for e, g in zip(Ainv[r], Ainv[c_])]
    # W = Ainv G_BV (d x nV), over the field
    W = [[zero] * nV for _ in range(d)]
    for k in range(d):
        for j in range(nV):
            acc = zero
            for l in range(d):
                acc = acc + Ainv[k][l] * els[cid[B[l]][j]]
            W[k][j] = acc % F.p
    Wm = coeff_mats(F, W, nV)
    pc = F.p.coeffs()
    lead = pc[D]
    redc = [-pc[k] / lead for k in range(D)]
    for i0 in range(0, nV, BLK):
        i1 = min(nV, i0 + BLK)
        Am = coeff_mats(F, [[els[cid[i][b]] for b in B] for i in range(i0, i1)], d)
        Sm = [None] * (2 * D - 1)
        for r in range(D):
            if not any(Am[r].entries()):
                continue
            for s in range(D):
                P = Am[r] * Wm[s]
                Sm[r + s] = P if Sm[r + s] is None else Sm[r + s] + P
        for t in range(2 * D - 2, D - 1, -1):          # x^t = x^(t-D) x^D, x^D = sum_k redc_k x^k
            T = Sm[t]
            if T is None:
                continue
            for k in range(D):
                if redc[k] != 0:
                    tgt = t - D + k
                    Sm[tgt] = T * redc[k] if Sm[tgt] is None else Sm[tgt] + T * redc[k]
            Sm[t] = None
        Gm = coeff_mats(F, [[els[cid[i][j]] for j in range(nV)] for i in range(i0, i1)], nV)
        for k in range(D):
            lhs = Sm[k] if Sm[k] is not None else flint.fmpq_mat(i1 - i0, nV)
            if not lhs == Gm[k]:
                return False, "the recognised Gram matrix does not have rank d"
    return True, ""


def span_exact_box(Xf, rad, d):
    """TRUSTED. d points of the box midpoint (pivoted QR, untrusted choice) with ||B^-1||_inf times the largest row sum
    of the radii below 1, B^-1 exact: every matrix within the box of B is nonsingular, so the points span R^d."""
    import scipy.linalg as sl
    _, _, piv = sl.qr(np.asarray(Xf).T, pivoting=True, mode="economic")
    idx = sorted(int(i) for i in piv[:d])
    Bm = flint.fmpq_mat(d, d, [flint.fmpq(*float(Xf[i][k]).as_integer_ratio()) for i in idx for k in range(d)])
    try:
        Bi = Bm.inv()
    except Exception:
        return False, float("inf")
    nrm = max(sum(abs(Fraction(int(Bi[i, k].p), int(Bi[i, k].q))) for k in range(d)) for i in range(d))
    rs = max(sum(Fraction(float(rad[i][k])) for k in range(d)) for i in idx)
    val = nrm * rs
    return bool(val < 1), XC.up_float(val)


# ============================================================================================ one entry
def strict_entry9(d, n, X, dec, tol, bits, maxdeg, log, datadir=None, probe=False):
    t0 = time.time()
    out = {}
    R9.survey(d, n, X, dec, tol, bits, maxdeg, 20000, log, out=out, probe=probe)
    if not str(log.get("status", "")).startswith("recognised"):
        return log
    log["recognition"] = log.pop("status")
    V, E = out["V"], out["E"]
    nV = len(V)
    Z, C, K = out["Z"], out["C"], out["K"]
    # ---- (F)
    try:
        F = Field9(out["field_poly"], out["theta"])
    except ValueError as e:
        log["status"] = "field: " + str(e)
        return log
    cel = F.elem(out["c_elem"])
    log["field_degree"] = F.D
    # ---- the recognised matrix: the elements of the clusters and the index matrix (cluster keys recomputed)
    shift, half = out["cluster_shift"], out["cluster_half"]
    els = [F.elem([(1, 1)])]
    key_id = {}
    for a, b, rc in out["clusters"]:
        g = sum(x * y for x, y in zip(Z[a], Z[b]))
        key_id[(g + half) >> shift] = len(els)
        els.append(F.elem(rc))
    cid = [[0] * nV for _ in range(nV)]
    for a in range(nV):
        za = Z[a]
        for b in range(a + 1, nV):
            k = key_id.get((sum(x * y for x, y in zip(za, Z[b])) + half) >> shift)
            if k is None:
                log["status"] = "a pair without a recognised inner product"
                return log
            cid[a][b] = cid[b][a] = k
    # ---- (E): contacts exactly c, every other pair provably below c*
    Eset = set((min(a, b), max(a, b)) for a, b in E)
    for (a, b) in E:
        if els[cid[a][b]] != cel:
            log["status"] = "a contact is not exactly c* in the recognised Gram matrix"
            return log
    other = set()
    for a in range(nV):
        for b in range(a + 1, nV):
            if (a, b) not in Eset:
                other.add(cid[a][b])
    worst = None
    for k in other:
        if els[k] == cel:
            log["status"] = "a pair outside E has inner product c*"
            return log
        diff = F.value(cel - els[k])
        if not diff > 0:
            log["status"] = "a pair outside E not provably below c*"
            return log
        lo = float(diff.lower())
        worst = lo if worst is None else min(worst, lo)
    log["gap_lower"] = worst
    cstar = F.value(cel)
    log["c_star"] = cstar.str(40, radius=False)
    # B: d well-conditioned points (untrusted choice)
    import scipy.linalg as sl
    Xf0 = np.array([[z / (1 << K) for z in row] for row in Z])
    _, _, piv = sl.qr(Xf0.T, pivoting=True, mode="economic")
    B = sorted(int(i) for i in piv[:d])
    # G_BB positive definite: ball Cholesky with positive pivots
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    val = [F.value(e) for e in els]
    Lc = [[flint.arb(0)] * d for _ in range(d)]
    for i in range(d):
        for j in range(i + 1):
            s_ = val[cid[B[i]][B[j]]] if i != j else flint.arb(1)
            for k in range(j):
                s_ = s_ - Lc[i][k] * Lc[j][k]
            if i == j:
                if not s_ > 0:
                    flint.ctx.prec = old
                    log["status"] = "G_BB not provably positive definite"
                    return log
                Lc[i][i] = s_.sqrt()
            else:
                Lc[i][j] = s_ / Lc[j][j]
    flint.ctx.prec = old
    t1 = time.time()
    ok, msg = rank_d_exact(F, els, cid, B, nV, d)
    log["seconds_gram_check"] = round(time.time() - t1, 1)
    if not ok:
        log["status"] = msg
        return log
    log["exact_gram"] = True
    # ---- enclosure of P* in the Cholesky gauge (balls), the box Pi
    flint.ctx.prec = PREC
    Yb = []
    for i in range(nV):
        y = [flint.arb(0)] * d
        for k in range(d):
            s_ = val[cid[i][B[k]]] if i != B[k] else flint.arb(1)
            for j in range(k):
                s_ = s_ - Lc[k][j] * y[j]
            y[k] = s_ / Lc[k][k]
        Yb.append(y)
    Xf = np.empty((nV, d))
    rad = np.empty((nV, d))
    for i in range(nV):
        for k in range(d):
            f = float(Yb[i][k].mid())
            Xf[i, k] = f
            rad[i, k] = math.nextafter(float(abs(Yb[i][k] - flint.arb(f)).upper()), math.inf)
    flint.ctx.prec = old
    rad = np.nextafter(rad * (1 + 2.0 ** -40), np.inf)
    log["enclosure_radius"] = float(rad.max())
    # ---- (S'): untrusted E', inverse and coefficients; the exact box test
    m = len(E)
    iso = nV * (d - 1) - d * (d - 1) // 2 + 1
    full = S.System(Xf, E, d)
    sub = cq.choose_square_subset_qr(full, Xf, d, nV, iso)
    Esq = [E[e] for e in sub]
    ss = set(sub)
    extra = [E[e] for e in range(m) if e not in ss]
    sysm = S.System(Xf, Esq, d)
    N = sysm.nu
    try:
        Rinv = cq.approx_inverse(sysm.jac(Xf), np.float64)
    except (np.linalg.LinAlgError, ValueError):
        log["status"] = "square Jacobian numerically singular at P*"
        return log
    beta_f = cq.choose_beta(Esq, extra, Xf, Rinv, N, nV, m, d, log)
    if beta_f is None or not np.all(np.asarray(beta_f) > 0):
        log["status"] = "no positive stress found at P*"
        return log
    Lx, zi = XC.dyadic_exp(Xf)
    Zc = [zi[a * d:(a + 1) * d] for a in range(nV)]
    sq = XC.Square(nV, d, Esq, Xf, Zc, 0, Lx)
    inv = XC.Inverse(np.asarray(Rinv, dtype=np.float64))
    t2 = time.time()
    res = XC.square_test(sq, inv, "box", rad=rad, bspec=(list(extra), np.asarray(beta_f, dtype=float)))
    XC._ABS.clear()
    log["seconds_box_test"] = round(time.time() - t2, 1)
    log["box_beta_K"] = XC.up_float(res["beta_K"])
    if not res.get("ok_nonsingular"):
        log["status"] = "square Jacobian not proved nonsingular on the box"
        return log
    lo = min(res["stress_lo"])
    wmin = min(lo, min(Fraction(float(x)) for x in beta_f))
    log.update(stress_min_lower=XC.dn_float(wmin), stress_min_rel=XC.dn_float(wmin * m))
    if not wmin > 0:
        log["status"] = "stress not provably positive on the box"
        return log
    ok_span, span_val = span_exact_box(Xf, rad, d)
    log["span_value"] = span_val
    if not ok_span:
        log["status"] = "span failed on the box"
        return log
    # ---- (G') rattlers in the frame of P*, and (D) the distance
    rat = np.setdiff1d(np.arange(n), V)
    U_, _, Vt_ = np.linalg.svd(X[V].T @ Xf)
    O = U_ @ Vt_
    Ycfg = np.empty((n, d))
    Ycfg[V] = Xf
    Ycfg[rat] = X[rat] @ O
    c_lo = Fraction(math.nextafter(float(cstar.lower()), -math.inf))
    if len(rat):
        Yr, rinfo = S.place_rattlers(Ycfg, V, Xf, float(c_lo), d)
        log["rattler_placement"] = rinfo
        ok_r, rbad = strict.rattler_gaps(Yr, V, rat, Xf, rad, c_lo, d)
        log["rattler_gaps_ok"] = ok_r
        if not ok_r:
            log["status"] = "rattler gaps failed"
            return log
    try:
        log.update(strict.distance_to_stored(dec, V, Yb, O, d))
    except Exception:
        log["distance_error"] = traceback.format_exc()[-300:]
    log["status"] = "STRICT"
    log["seconds_strict9"] = round(time.time() - t0, 1)
    if datadir:
        os.makedirs(datadir, exist_ok=True)
        desc = {"d": d, "n": n, "field_poly": [int(c) for c in out["field_poly"]],
                "theta_ball": F.theta.str(60, radius=True), "c_elem": [str(q) for q in cel.coeffs()],
                "c_star": log["c_star"], "nonrattlers": [int(v) for v in V], "B": B,
                "elements": [[str(q) for q in e.coeffs()] for e in els],
                "index_matrix": cid, "contacts": [[int(a), int(b)] for a, b in E]}
        with gzip.open(os.path.join(datadir, f"d{d}_n{n}.json.gz"), "wt") as fo:
            json.dump(desc, fo)
    return log


def main():
    outdir = sys.argv[1]
    args = sys.argv[2:]
    opt = lambda k, dflt: args[args.index(k) + 1] if k in args else dflt
    bits, maxdeg = int(opt("--bits", "2000")), int(opt("--maxdeg", "30"))
    datadir = opt("--data", None)
    cls = opt("--class", "local")
    probe = "--probe" in args
    only = set(opt("--only", "").split(",")) if "--only" in args else None
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    todo = []
    for d, n in R9.open_entries(cls):
        if only and f"{d}_{n}" not in only:
            continue
        tol, U = R9.census_tol(d, n)
        todo.append((U, d, n, tol))
    todo.sort()
    os.makedirs(outdir, exist_ok=True)
    for U, d, n, tol in todo:
        js = os.path.join(outdir, f"d{d}_n{n}.json")
        if os.path.exists(js):
            continue
        log = {"d": d, "n": n, "tol": tol, "unknowns": U, "bits": bits, "maxdeg": maxdeg, "probe": probe,
               "program_sha256": SHA}
        t0, c0 = time.time(), time.process_time()
        try:
            pts = snapshot.read_points(R9.SNAP, idx[(d, n)])
            dec = [[t.strip() for t in p] for p in pts]
            X = np.array([[float(x) for x in p] for p in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            strict_entry9(d, n, X, dec, tol, bits, maxdeg, log, datadir, probe=probe)
        except Exception:
            log["status"] = "error"
            log["error"] = traceback.format_exc()[-1500:]
        log["wall"] = round(time.time() - t0, 1)
        log["cpu_seconds"] = round(time.process_time() - c0, 1)
        with open(js + ".part", "w") as fo:
            json.dump(log, fo, default=str)
        os.replace(js + ".part", js)
        print(f"d{d}_n{n} U={U} {log['status']} deg={log.get('field_degree')} {log['wall']}s", flush=True)


if __name__ == "__main__":
    main()
