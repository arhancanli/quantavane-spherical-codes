"""The negative certificate (Theorem 4.1) redone in exact arithmetic on a stratified sample (census version 7).

  nice -n 10 python3 -I flexexact.py OUT.jsonl SAMPLE.tsv

SAMPLE.tsv (from flexsample.py) has the columns d, n, category, reason. For each entry, with the tolerance of its census
FLEX attempt:

UNTRUSTED (proposals only): E, the pairs within tau of the largest cosine of the double-precision code X, rebuilt
exactly as cq.certify_flex builds it, with the SHA-256 of its pair list (the same hash as flexgram7.py); right
inverses of L(X) from the normal equations and, when these give beta >= 1/2, from a QR factorisation (as the census).

TRUSTED, exact (FLINT fmpz matrices, Python integers and rationals; Arb balls only for square roots):
  - L(X) has the doubles of X as entries, scaled to integers; each right inverse is rounded to C' = Cint 2^-t (60
    bits for its largest entry; any matrix will do), the one with the smaller exact beta is kept; L C' is formed exactly in column blocks, and beta = the largest row
    sum of |I - L C'| and ||C'||_inf are exact rationals; if beta < 1, h0 = (1 - beta) / (2 ||C'||_inf);
  - h = h0 - max_i |xhat_i - x_i|_1, with xhat_i the stored decimals normalised exactly (Arb, 256 bits);
  - tau' = (largest inner product of xhat over E, from below) - (largest over the pairs outside E, from above), every
    inner product from the exact integer Gram matrix of the stored decimals and the Arb norms;
  - delta = min(h / s, tau' / 4) with a rational s >= sqrt(d).
Then Theorem 4.1 holds with this delta for the exactly normalised stored code. The census values (h of the double
kernel) and flexgram7.py (tau', delta, E hash) are compared afterwards; nothing here uses them.
Entries with |rows| * |columns| of L above MAX_CELLS are not attempted (memory).
"""
import hashlib, json, math, os, sys, time, traceback
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "verify"))
import exactcheck as XC          # sets the library paths
import numpy as np
import flint
import snapshot

SHA = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
APREC = 256
CBITS = 60
BLOCK = 256
MAX_CELLS = 12_000_000


def arb_dn(a):
    lo = a.lower()
    m, e = lo.mid().man_exp()
    return Fraction(int(m)) * (Fraction(2) ** int(e))


def dec_matrix(dec):
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
    ints = []
    for row in fr:
        out = []
        for q in row:
            v = q * sc
            assert v.denominator == 1
            out.append(int(v))
        ints.append(out)
    return ints


def proposal(X, n, tol):
    G = np.empty((n, n))
    for i0 in range(0, n, 256):
        G[i0:i0 + 256] = X[i0:i0 + 256] @ X.T
    np.fill_diagonal(G, -2.0)
    M = float(G.max())
    I_, J_ = np.nonzero(np.triu(G >= M - tol, 1))
    return I_, J_


def check(X, dec, d, n, tol, a):
    rec = {}
    I_, J_ = proposal(X, n, tol)
    m = len(I_)
    pts = np.unique(np.concatenate([I_, J_]))
    nE = len(pts)
    rec["E_sha256"] = hashlib.sha256(json.dumps([[int(i), int(j)] for i, j in zip(I_, J_)]).encode()).hexdigest()
    rec.update(contacts=int(m), points_in_contacts=int(nE),
               proposal_matches_census=(m == a["contacts"] and nE == a["points_in_contacts"]))
    UL, UC = nE + m, nE * d
    rec.update(rows=UL, columns=UC)
    if UL * UC > MAX_CELLS:
        rec["status"] = "not attempted: above MAX_CELLS"
        return rec
    loc = -np.ones(n, int)
    loc[pts] = np.arange(nE)
    # ---- L(X) exactly, scaled by 2^Lx
    Lx, xi = XC.dyadic_exp(X)
    Lm = flint.fmpz_mat(UL, UC)
    for r_, i in enumerate(pts):
        for k in range(d):
            Lm[r_, int(loc[i]) * d + k] = 2 * xi[int(i) * d + k]
    for e, (i, j) in enumerate(zip(I_, J_)):
        r_ = nE + e
        for k in range(d):
            Lm[r_, int(loc[i]) * d + k] = xi[int(j) * d + k]
            Lm[r_, int(loc[j]) * d + k] = xi[int(i) * d + k]
    # ---- a right inverse (untrusted): normal equations in double precision
    import scipy.sparse as sps
    ar = np.arange(d)
    rN = np.repeat(np.arange(nE), d)
    cN = (loc[pts][:, None] * d + ar).ravel()
    vN = 2.0 * X[pts].ravel()
    rC = np.repeat(nE + np.arange(m), 2 * d)
    cC = np.concatenate([loc[I_][:, None] * d + ar, loc[J_][:, None] * d + ar], 1).ravel()
    vC = np.concatenate([X[J_], X[I_]], 1).ravel()
    Lf = sps.csr_matrix((np.concatenate([vN, vC]), (np.concatenate([rN, rC]), np.concatenate([cN, cC]))),
                        shape=(UL, UC))
    def exact_beta(Cf):
        """beta = max row sum of |I - L C'| and ||C'||_inf, exactly, for C' = Cf rounded to 60-bit integers / 2^t."""
        mx = float(np.max(np.abs(Cf)))
        if not (mx > 0 and math.isfinite(mx)):
            return None, None
        t = CBITS - math.frexp(mx)[1]
        Cint = np.rint(np.ldexp(Cf, t)).astype(np.int64)
        A = np.abs(Cint)                         # |Cint| = hi 2^30 + lo: the int64 row sums cannot overflow
        hi, lo = (A >> 30).sum(1), (A & ((1 << 30) - 1)).sum(1)
        del A
        cnum = max(int(h_) * (1 << 30) + int(l_) for h_, l_ in zip(hi, lo))
        den = 1 << (Lx + t)
        rows = [0] * UL
        for c0 in range(0, UL, BLOCK):
            c1 = min(UL, c0 + BLOCK)
            P = Lm * flint.fmpz_mat(Cint[:, c0:c1].tolist())
            ent = P.entries()
            w = c1 - c0
            for i in range(UL):
                base = i * w
                tot = 0
                for jj in range(w):
                    x = int(ent[base + jj])
                    if c0 + jj == i:
                        x -= den
                    tot += x if x >= 0 else -x
                rows[i] += tot
            del P, ent
        return Fraction(max(rows), den), Fraction(cnum, 1 << t)

    # ---- right inverses (untrusted), as the census: normal equations, and QR when that is not good enough
    Gm = (Lf @ Lf.T).toarray()
    beta, cnorm, how = None, None, None
    try:
        Cf = np.asarray(Lf.T @ np.linalg.inv(Gm))
        beta, cnorm = exact_beta(Cf)
        how = "normal equations"
        del Cf
    except np.linalg.LinAlgError:
        pass
    del Gm
    if (beta is None or not beta < Fraction(1, 2)) and UL * UC <= MAX_CELLS:
        import scipy.linalg as sl
        Qf, Rf = sl.qr(Lf.T.toarray(), mode="economic", check_finite=False)
        try:
            Cf = sl.solve_triangular(Rf, Qf.T, lower=False, check_finite=False).T
            b2, c2 = exact_beta(Cf)
            if b2 is not None and (beta is None or b2 < beta):
                beta, cnorm, how = b2, c2, "QR"
        except (np.linalg.LinAlgError, ValueError):
            pass
        del Qf, Rf
    if beta is None:
        rec["status"] = "FAIL: no right inverse"
        return rec
    rec["right_inverse"] = how
    rec["beta"] = XC.up_float(beta)
    rec["C_norm_inf"] = XC.up_float(cnorm)
    if not beta < 1:
        rec["status"] = "FAIL: beta >= 1"
        return rec
    h0 = (1 - beta) / (2 * cnorm)
    # ---- the exactly normalised stored code
    D = dec_matrix(dec)
    old = flint.ctx.prec
    flint.ctx.prec = APREC
    try:
        nrm = [flint.arb(sum(v * v for v in D[i])).sqrt() for i in range(n)]
        dist1 = Fraction(0)
        for i in range(n):
            tot = flint.arb(0)
            for k in range(d):
                xq = flint.fmpq(*float(X[i, k]).as_integer_ratio())
                tot += abs(flint.arb(D[i][k]) / nrm[i] - flint.arb(xq))
            dist1 = max(dist1, XC._arb_up(tot))
        h = h0 - dist1
        Dm = flint.fmpz_mat(D)
        gent = (Dm * Dm.transpose()).entries()
        inv_n = [1 / x for x in nrm]
        Eset = set(zip(map(int, I_), map(int, J_)))
        m_lo, u_out = None, None
        for i in range(n):
            for j in range(i + 1, n):
                v = flint.arb(int(gent[i * n + j])) * inv_n[i] * inv_n[j]
                if (i, j) in Eset:
                    q = arb_dn(v)
                    if m_lo is None or q > m_lo:
                        m_lo = q
                else:
                    q = XC._arb_up(v)
                    if u_out is None or q > u_out:
                        u_out = q
    finally:
        flint.ctx.prec = old
    taup = (m_lo - u_out) if u_out is not None else Fraction(2)
    s = XC.sqrt_up(Fraction(d))
    delta = min(h / s, taup / 4)
    rec.update(h0=XC.dn_float(h0), dist_l1=XC.up_float(dist1), h=XC.dn_float(h), tau_prime=XC.dn_float(taup),
               delta=XC.dn_float(delta), delta_by="tau'" if taup / 4 <= h / s else "h",
               h_census=a.get("flex_radius_l1_stored"))
    rec["status"] = "FLEX" if (h > 0 and taup > 0) else "FAIL"
    return rec


def main():
    outp, sample = sys.argv[1], sys.argv[2]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    done = set()
    if os.path.exists(outp):
        done = {(r["d"], r["n"]) for r in map(json.loads, open(outp))}
    for line in open(sample):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        d, n, cat = int(p[0]), int(p[1]), p[2]
        if (d, n) in done:
            continue
        r = json.load(open(os.path.join(XC.CERT, "census", "results_v4", f"d{d}_n{n}.json")))
        a = [x for x in r["attempts"] if x.get("status") == "FLEX"][0]
        t0, c0 = time.time(), time.process_time()
        rec = {"d": d, "n": n, "category": cat, "reason": p[3] if len(p) > 3 else "", "tol": a["tol"],
               "program_sha256": SHA}
        try:
            pts = snapshot.read_points(XC.SNAP, idx[(d, n)])
            dec = [[t_.strip() for t_ in pp] for pp in pts]
            X = np.array([[float(x) for x in pp] for pp in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            rec.update(check(X, dec, d, n, a["tol"], a))
        except Exception:
            rec["status"] = "error"
            rec["error"] = traceback.format_exc()[-1500:]
        rec["wall"] = round(time.time() - t0, 1)
        rec["cpu_seconds"] = round(time.process_time() - c0, 1)
        with open(outp, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(d, n, cat, rec.get("status"), rec.get("delta"), rec["wall"], flush=True)


if __name__ == "__main__":
    main()
