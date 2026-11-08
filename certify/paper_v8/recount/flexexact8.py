"""The negative certificate (Theorem 4.1) decided in exact arithmetic for every entry, with sparse exact products
(census version 8). This is the decision of record for the negative verdict.

  nice -n 10 python3 -I flexexact8.py OUTDIR [--only d_n,...]

Several processes may run at once on the same OUTDIR: an entry is claimed by creating OUTDIR/dD_nN.lock exclusively.

UNTRUSTED (proposals only): E, the pairs within tau of the largest cosine of the double-precision code X (exactly as
cq.certify_flex), with the SHA-256 of its pair list; approximate right inverses of L (normal equations, and a QR
factorisation when these give beta >= 1/2 and L has at most 1.2e7 entries, as the census); approximate inner products.

TRUSTED, exact (Python and NumPy integers whose bounds are proved below, FLINT integers, Arb balls for square roots):
  - the reference configuration X0 = X rounded to the grid 2^-52, so L0 = L(X0) has integer entries L_int / 2^52 with
    |L_int| <= 2^53 and at most 2d nonzeros per row;
  - C' = C_int / 2^t with |C_int| < 2^52 (the approximate right inverse, rounded; any matrix will do);
  - L_int C_int exactly: L_int = L1 2^26 + L0, C_int = C1 2^26 + C0 with 0 <= L0, C0 < 2^26, |L1| <= 2^27, |C1| <= 2^26;
    every product of limbs is at most 2^53 in absolute value and every row of L has at most 2d <= 64 nonzeros, so the
    sparse integer products L_a C_b (int64) are exact (|sums| <= 2^59, two of them added <= 2^60); their combination is
    normalised into base-2^26 digits, and |I - L0 C'| is summed row by row in exact integers;
  - beta = the largest row sum of |I - L0 C'| and ||C'||_inf exactly; h0 = (1 - beta) / (2 ||C'||_inf) when beta < 1;
  - h = h0 - max_i |xhat_i - x0_i|_1, xhat the exactly normalised stored decimals (Arb, 256 bits);
  - tau' = M_lo - U: M_lo a lower bound (Arb) of the inner product of xhat at the best pair of E; U a rational with
    <xhat_i, xhat_j> <= U for every pair outside E, checked pair by pair in exact integers from the integer Gram matrix
    of the stored decimals (G_ij <= 0, or G_ij^2 Ud^2 <= Un^2 N_i N_j);
  - delta = min(h / s, tau' / 4) with a rational s >= sqrt(d); FLEX when h > 0 and tau' > 0.
"""
import hashlib, json, math, os, sys, time, traceback
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
sys.path.insert(0, os.path.join(CERT, "census"))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import scipy.sparse as sps
import flint
import snapshot

SHA = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
APREC = 256
GRID = 52
M26 = (1 << 26) - 1
QR_MAX = 12_000_000


def F2(num, e):
    return Fraction(num, 1 << e) if e >= 0 else Fraction(num * (1 << -e))


def up_float(q):
    f = float(q)
    return f if Fraction(f) >= q else math.nextafter(f, math.inf)


def dn_float(q):
    f = float(q)
    return f if Fraction(f) <= q else math.nextafter(f, -math.inf)


def arb_up(a):
    m, e = a.upper().mid().man_exp()
    return Fraction(int(m)) * (Fraction(2) ** int(e))


def arb_dn(a):
    m, e = a.lower().mid().man_exp()
    return Fraction(int(m)) * (Fraction(2) ** int(e))


def sqrt_up(q):
    old = flint.ctx.prec
    flint.ctx.prec = APREC
    r = arb_up(flint.arb(flint.fmpq(q.numerator, q.denominator)).sqrt())
    flint.ctx.prec = old
    assert r * r >= q
    return r


def dec_ints(dec):
    """The stored decimals as integers over a common 10^k, exactly."""
    fr = [[Fraction(t) for t in row] for row in dec]
    k = 0
    for row in fr:
        for q in row:
            den = q.denominator
            a = (den & -den).bit_length() - 1
            rest, b = den >> a, 0
            while rest % 5 == 0:
                rest //= 5
                b += 1
            assert rest == 1, "not a decimal"
            k = max(k, a, b)
    sc = 10 ** k
    out = []
    for row in fr:
        r_ = []
        for q in row:
            v = q * sc
            assert v.denominator == 1
            r_.append(int(v))
        out.append(r_)
    return out


def proposal(X, n, tol):
    G = np.empty((n, n))
    for i0 in range(0, n, 256):
        G[i0:i0 + 256] = X[i0:i0 + 256] @ X.T
    np.fill_diagonal(G, -2.0)
    M = float(G.max())
    I_, J_ = np.nonzero(np.triu(G >= M - tol, 1))
    return I_, J_, G


def build_L(Xi, Xf, pts, loc, I_, J_, d):
    """L(X0) as int64 sparse (entries L_int, |L_int| <= 2^53) and as float (L_int / 2^52), rows: points of W, then E."""
    nE, m = len(pts), len(I_)
    ar = np.arange(d)
    rN = np.repeat(np.arange(nE), d)
    cN = (loc[pts][:, None] * d + ar).ravel()
    vN = 2 * Xi[pts].ravel()
    rC = np.repeat(nE + np.arange(m), 2 * d)
    cC = np.concatenate([loc[I_][:, None] * d + ar, loc[J_][:, None] * d + ar], 1).ravel()
    vC = np.concatenate([Xi[J_], Xi[I_]], 1).ravel()
    rows, cols = np.concatenate([rN, rC]), np.concatenate([cN, cC])
    vals = np.concatenate([vN, vC]).astype(np.int64)
    shape = (nE + m, nE * d)
    Li = sps.csr_matrix((vals, (rows, cols)), shape=shape, dtype=np.int64)
    Li.sum_duplicates()
    assert int(np.max(np.abs(Li.data))) <= (1 << 53)
    Lf = sps.csr_matrix((Li.data.astype(np.float64) / 2.0 ** GRID, Li.indices, Li.indptr), shape=shape)
    return Li, Lf


def exact_beta(Li, L1, L0, Cblock, ncols_C, UL, UC):
    """beta and ||C'||_inf exactly for C' = round(C 2^t) / 2^t, C given by column blocks Cblock(j0, j1) (floats).
    Returns (beta, cnorm, t) or (None, None, None)."""
    BLK = 256
    mx = 0.0
    for j0 in range(0, UL, BLK):
        mx = max(mx, float(np.max(np.abs(Cblock(j0, min(UL, j0 + BLK))))))
    if not (mx > 0 and math.isfinite(mx)):
        return None, None, None
    t = 52 - math.frexp(mx)[1]                 # |C_int| < 2^52
    if not (0 <= t <= 60):
        return None, None, None
    rowA = [0] * UL
    rowBh = np.zeros(UL, dtype=np.int64)
    rowBl = np.zeros(UL, dtype=np.int64)
    rowA_np = np.zeros(UL, dtype=np.int64)
    cs_hi = np.zeros(UC, dtype=np.int64)
    cs_lo = np.zeros(UC, dtype=np.int64)
    two52 = np.int64(1) << np.int64(52)
    for j0 in range(0, UL, BLK):
        j1 = min(UL, j0 + BLK)
        Ci = np.rint(np.ldexp(Cblock(j0, j1), t)).astype(np.int64)
        A_ = np.abs(Ci)
        cs_hi += (A_ >> 26).sum(1)
        cs_lo += (A_ & M26).sum(1)
        del A_
        C1 = Ci >> 26
        C0 = Ci & M26
        del Ci
        Q2 = np.asarray(L1 @ C1)
        Q1 = np.asarray(L1 @ C0) + np.asarray(L0 @ C1)
        Q0 = np.asarray(L0 @ C0)
        del C1, C0
        d0 = Q0 & M26
        k0 = Q0 >> 26
        c1 = Q1 + k0
        d1 = c1 & M26
        c2 = Q2 + (c1 >> 26)
        del Q0, Q1, Q2, k0, c1
        # E = 2^(52+t) delta - V, V = c2 2^52 + d1 2^26 + d0 (0 <= d0, d1 < 2^26)
        e2 = -c2
        w = np.arange(j1 - j0)
        rr = j0 + w
        e2[rr, w] += np.int64(1) << np.int64(t)
        low = (d1 << 26) + d0                  # 0 <= low < 2^52
        pos = e2 > 0
        A_ = np.where(pos, e2 - 1, -e2)
        B_ = np.where(pos, two52 - low, low)
        carry = B_ == two52
        A_ = A_ + carry
        B_ = np.where(carry, 0, B_)
        if int(A_.max()) > (1 << 40):
            return Fraction(10 ** 9), Fraction(1), t     # beta >= 1 by far
        rowA_np += A_.sum(1)
        rowBh += (B_ >> 26).sum(1)
        rowBl += (B_ & M26).sum(1)
        del e2, low, A_, B_, pos, carry, d0, d1, c2
    rows = [int(a) * (1 << 52) + int(bh) * (1 << 26) + int(bl) for a, bh, bl in zip(rowA_np, rowBh, rowBl)]
    beta = Fraction(max(rows), 1 << (52 + t))
    cn = max(int(h_) * (1 << 26) + int(l_) for h_, l_ in zip(cs_hi, cs_lo))
    return beta, Fraction(cn, 1 << t), t


def check(X, dec, d, n, tol, a):
    rec = {}
    I_, J_, Gf = proposal(X, n, tol)
    m = len(I_)
    pts = np.unique(np.concatenate([I_, J_]))
    nE = len(pts)
    rec["E_sha256"] = hashlib.sha256(json.dumps([[int(i), int(j)] for i, j in zip(I_, J_)]).encode()).hexdigest()
    rec.update(contacts=int(m), points_in_contacts=int(nE),
               proposal_matches_census=(m == a["contacts"] and nE == a["points_in_contacts"]))
    UL, UC = nE + m, nE * d
    rec.update(rows=UL, columns=UC)
    loc = -np.ones(n, int)
    loc[pts] = np.arange(nE)
    Xi = np.rint(np.ldexp(X, GRID)).astype(np.int64)          # X0 = Xi / 2^52
    Li, Lf = build_L(Xi, X, pts, loc, I_, J_, d)
    L1 = Li.copy()
    L1.data = Li.data >> 26
    L0 = Li.copy()
    L0.data = Li.data & M26
    LfT = Lf.T.tocsr()
    # ---- right inverses (untrusted), as the census
    beta, cnorm, how = None, None, None
    try:
        Gm = (Lf @ Lf.T).toarray()
        R = np.linalg.inv(Gm)
        del Gm
        beta, cnorm, t = exact_beta(Li, L1, L0, lambda j0, j1: np.asarray(LfT @ R[:, j0:j1]), UL, UL, UC)
        how = "normal equations"
        del R
    except np.linalg.LinAlgError:
        pass
    if (beta is None or not beta < Fraction(1, 2)) and UL * UC <= QR_MAX:
        import scipy.linalg as sl
        try:
            Qf, Rf = sl.qr(Lf.T.toarray(), mode="economic", check_finite=False)
            Cf = sl.solve_triangular(Rf, Qf.T, lower=False, check_finite=False).T
            del Qf, Rf
            b2, c2, t2 = exact_beta(Li, L1, L0, lambda j0, j1: Cf[:, j0:j1], UL, UL, UC)
            if b2 is not None and (beta is None or b2 < beta):
                beta, cnorm, how = b2, c2, "QR"
            del Cf
        except (np.linalg.LinAlgError, ValueError):
            pass
    if beta is None:
        rec["status"] = "FAIL: no right inverse"
        return rec
    rec.update(right_inverse=how, beta=up_float(beta), C_norm_inf=up_float(cnorm))
    if not beta < 1:
        rec["status"] = "FAIL: beta >= 1"
        return rec
    h0 = (1 - beta) / (2 * cnorm)
    # ---- the exactly normalised stored code: distance to X0, and tau'
    D = dec_ints(dec)
    old = flint.ctx.prec
    flint.ctx.prec = APREC
    try:
        N = [sum(v * v for v in row) for row in D]
        nrm = [flint.arb(x).sqrt() for x in N]
        dist1 = Fraction(0)
        for i in pts:
            i = int(i)
            tot = flint.arb(0)
            for k in range(d):
                tot += abs(flint.arb(D[i][k]) / nrm[i] - flint.arb(flint.fmpq(int(Xi[i, k]), 1 << GRID)))
            dist1 = max(dist1, arb_up(tot))
        h = h0 - dist1
        Gd = flint.fmpz_mat(D)
        G = [int(x) for x in (Gd * Gd.transpose()).entries()]
        inE = np.zeros((n, n), bool)
        inE[I_, J_] = True
        inE[J_, I_] = True
        np.fill_diagonal(inE, True)
        # untrusted: the best pair of E and the largest pair outside E (double precision)
        Xh = np.array([[float(Fraction(t_)) for t_ in row] for row in dec])
        Xh /= np.linalg.norm(Xh, axis=1, keepdims=True)
        Gh = Xh @ Xh.T
        ib, jb = int(I_[np.argmax(Gh[I_, J_])]), int(J_[np.argmax(Gh[I_, J_])])
        M_lo = arb_dn(flint.arb(G[ib * n + jb]) / (nrm[ib] * nrm[jb]))
        out = np.where(inE, -np.inf, Gh)
        u_f = float(out.max()) if np.isfinite(out.max()) else None
    finally:
        flint.ctx.prec = old
    if u_f is None:
        taup = Fraction(2)
        rec["pairs_checked"] = 0
    else:
        # U0: the largest inner product outside E in double precision (untrusted); every pair is then checked exactly
        # against U0, and the pairs above it (if any) get Arb upper bounds; U is the largest of these and U0
        U0 = Fraction(u_f)
        Un, Ud = U0.numerator, U0.denominator
        Un2, Ud2 = Un * Un, Ud * Ud
        above = []
        cnt = 0
        for i in range(n):
            UNi = Un2 * N[i]
            row = G[i * n:(i + 1) * n]
            js = (np.nonzero(~inE[i, i + 1:])[0] + (i + 1)).tolist()
            cnt += len(js)
            if U0 >= 0:
                for j in js:
                    g = row[j]
                    if g > 0 and g * g * Ud2 > UNi * N[j]:
                        above.append((i, j))
            else:
                for j in js:
                    g = row[j]
                    if not (g < 0 and g * g * Ud2 >= UNi * N[j]):
                        above.append((i, j))
        rec["pairs_checked"] = cnt
        rec["pairs_above_double_max"] = len(above)
        U = U0
        if above:
            flint.ctx.prec = APREC
            try:
                for i, j in above:
                    U = max(U, arb_up(flint.arb(G[i * n + j]) / (nrm[i] * nrm[j])))
            finally:
                flint.ctx.prec = old
        taup = M_lo - U
    s = sqrt_up(Fraction(d))
    delta = min(h / s, taup / 4)
    rec.update(h0=dn_float(h0), dist_l1=up_float(dist1), h=dn_float(h), tau_prime=dn_float(taup),
               delta=dn_float(delta), delta_by="tau'" if taup / 4 <= h / s else "h")
    rec["status"] = "FLEX" if (h > 0 and taup > 0) else "FAIL"
    return rec


def main():
    out = sys.argv[1]
    only = None
    if "--only" in sys.argv:
        only = set(sys.argv[sys.argv.index("--only") + 1].split(","))
    os.makedirs(out, exist_ok=True)
    s5 = json.load(open(os.path.join(CERT, "v5", "out", "summary_v5.json")))
    keys = sorted((tuple(map(int, k.split(","))) for k, v in s5["class"].items() if v[1] == "flex"),
                  key=lambda k: -k[0] * k[1])
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    for d, n in keys:
        if only and f"{d}_{n}" not in only:
            continue
        js = os.path.join(out, f"d{d}_n{n}.json")
        lk = js[:-5] + ".lock"
        if os.path.exists(js):
            continue
        try:
            fd = os.open(lk, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.close(fd)
        except FileExistsError:
            continue
        r = json.load(open(os.path.join(CERT, "census", "results_v4", f"d{d}_n{n}.json")))
        a = [x for x in r["attempts"] if x.get("status") == "FLEX"][0]
        t0, c0 = time.time(), time.process_time()
        rec = {"d": d, "n": n, "tol": a["tol"], "program_sha256": SHA}
        try:
            pts = snapshot.read_points(SNAP, idx[(d, n)])
            dec = [[t_.strip() for t_ in pp] for pp in pts]
            X = np.array([[float(x) for x in pp] for pp in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            rec.update(check(X, dec, d, n, a["tol"], a))
        except MemoryError:
            rec["status"] = "error: memory"
        except Exception:
            rec["status"] = "error"
            rec["error"] = traceback.format_exc()[-1500:]
        rec["wall"] = round(time.time() - t0, 2)
        rec["cpu_seconds"] = round(time.process_time() - c0, 2)
        with open(js + ".part", "w") as f:
            json.dump(rec, f)
        os.replace(js + ".part", js)
        os.remove(lk)
        print(d, n, rec.get("status"), rec.get("delta"), rec.get("rows"), rec["cpu_seconds"], flush=True)


if __name__ == "__main__":
    main()
