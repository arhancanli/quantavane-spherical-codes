"""Census version 9, item B: the stand-alone verifier of the deposited certificate data (files of export9.py).

  python3 -I verify_cert9.py OUT.json FILE.npz [FILE.npz ...]

No census kernel, LAPACK or floating-point step enters: the files hold the untrusted data exactly (integers, and
binary64 arrays, which are exact rationals), and every decision is made in exact integer and rational arithmetic and in
Arb balls by released trusted code.

Positive certificates: exactcheck.check_entry (paper_v8/verify/exactcheck.py, released, unchanged; exactcheck2.py for
the three records decided with the second stress enclosure) decides the file.
check_entry obtains its untrusted data by calling cq.certify; for this program only, the module cq it imports is a
reader of the file (StoredData below), which fills the same dictionary from the file and calls nothing else.
Negative certificates: Theorem 4.1 with the stored right inverse C' = Cint / 2^t: beta and ||C'||_inf by
flexexact8.exact_beta (released, unchanged; the stored matrix is handed to it as binary64, which is exact since
|Cint| < 2^52), and h, tau' and delta by check_negative below, a copy of the trusted part of flexexact8.check that
follows the choice of the right inverse.
Each result is compared with the decision of record (paper_v8 outputs), field by field.
"""
import hashlib, json, math, os, sys, time, types
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
sys.path.insert(0, os.path.join(CERT, "paper_v8", "verify"))
import exactcheck as XC                       # sets the library paths
import exactcheck2 as XC2                     # the second stress enclosure (Lemma 3.4(b)), for its three records
sys.path.insert(0, os.path.join(CERT, "paper_v8", "recount"))
import flexexact8 as FX
import numpy as np
import flint

SHA = {"verify_cert9.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
       "exactcheck.py": XC.SHA, "exactcheck2.py": XC2.SHA, "flexexact8.py": FX.SHA}


class StoredData:
    """Stands in for the module cq inside exactcheck.check_entry: certify() fills keep from the file."""

    def __init__(self, meta, arr):
        self.meta, self.arr = meta, arr

    def certify(self, X, d, n, tol, dec=None, keep=None, **kw):
        m, a = self.meta, self.arr
        keep.update(V=[int(v) for v in m["V"]], E=[tuple(e) for e in m["E"]], Esq=[tuple(e) for e in m["Esq"]],
                    extra=[tuple(e) for e in m["extra"]], Z=[[int(z) for z in row] for row in m["Z"]],
                    C=int(m["C"]), K=int(m["K"]), R=a["R"], xs=a["xs"], beta_f=a["beta_f"], Y=a["Y"])
        return {"status": m["status"]}


def check_positive(meta, arr):
    d, n = meta["d"], meta["n"]
    dec = meta["dec"]
    stub = types.ModuleType("cq")
    stub.certify = StoredData(meta, arr).certify
    saved = sys.modules.get("cq")
    sys.modules["cq"] = stub
    try:
        log = {}
        prog = XC2 if meta.get("record_dir") == "out_exact2" else XC
        prog.check_entry(None, dec, d, n, meta["tol"], meta["status"], log)
    finally:
        if saved is not None:
            sys.modules["cq"] = saved
        else:
            del sys.modules["cq"]
    rec = meta["record_v8"]
    keys = [k for k in ("status", "c_interval_30", "r", "beta_K", "stress_min_lower", "kappa", "rho", "tube",
                        "dist_ls", "gap_lower", "code_inside_rho") if k in rec or k in log]
    same = {k: log.get(k) == rec.get(k) for k in keys}
    return log, same


def check_negative(meta, arr):
    """The trusted part of flexexact8.check after the right inverse, with the stored one."""
    d, n, dec = meta["d"], meta["n"], meta["dec"]
    t = int(meta["t"])
    Cint = arr["Cint"]
    assert int(np.max(np.abs(Cint))) < (1 << 52)
    I_ = np.array([p[0] for p in meta["pairs"]], dtype=int)
    J_ = np.array([p[1] for p in meta["pairs"]], dtype=int)
    pts = np.unique(np.concatenate([I_, J_]))
    nE = len(pts)
    m = len(I_)
    UL, UC = nE + m, nE * d
    loc = -np.ones(n, int)
    loc[pts] = np.arange(nE)
    Xi = arr["X0"]
    Li, _ = FX.build_L(Xi, None, pts, loc, I_, J_, d)
    L1 = Li.copy()
    L1.data = Li.data >> 26
    L0 = Li.copy()
    L0.data = Li.data & FX.M26
    Cf = np.ldexp(Cint.astype(np.float64), -t)                 # exact: |Cint| < 2^52
    assert Cf.shape == (UC, UL)
    beta, cnorm, t2 = FX.exact_beta(Li, L1, L0, lambda j0, j1: Cf[:, j0:j1], UL, UL, UC)
    assert t2 == t
    rec = {"rows": UL, "columns": UC, "beta": FX.up_float(beta), "C_norm_inf": FX.up_float(cnorm)}
    if not beta < 1:
        rec["status"] = "FAIL: beta >= 1"
        return rec
    h0 = (1 - beta) / (2 * cnorm)
    D = FX.dec_ints(dec)
    old = flint.ctx.prec
    flint.ctx.prec = FX.APREC
    try:
        N = [sum(v * v for v in row) for row in D]
        nrm = [flint.arb(x).sqrt() for x in N]
        dist1 = Fraction(0)
        for i in pts:
            i = int(i)
            tot = flint.arb(0)
            for k in range(d):
                tot += abs(flint.arb(D[i][k]) / nrm[i] - flint.arb(flint.fmpq(int(Xi[i, k]), 1 << FX.GRID)))
            dist1 = max(dist1, FX.arb_up(tot))
        h = h0 - dist1
        Gd = flint.fmpz_mat(D)
        G = [int(x) for x in (Gd * Gd.transpose()).entries()]
        inE = np.zeros((n, n), bool)
        inE[I_, J_] = True
        inE[J_, I_] = True
        np.fill_diagonal(inE, True)
        Xh = np.array([[float(Fraction(t_)) for t_ in row] for row in dec])
        Xh /= np.linalg.norm(Xh, axis=1, keepdims=True)
        Gh = Xh @ Xh.T
        ib, jb = int(I_[np.argmax(Gh[I_, J_])]), int(J_[np.argmax(Gh[I_, J_])])
        M_lo = FX.arb_dn(flint.arb(G[ib * n + jb]) / (nrm[ib] * nrm[jb]))
        out = np.where(inE, -np.inf, Gh)
        u_f = float(out.max()) if np.isfinite(out.max()) else None
    finally:
        flint.ctx.prec = old
    if u_f is None:
        taup = Fraction(2)
    else:
        U0 = Fraction(u_f)
        Un, Ud = U0.numerator, U0.denominator
        Un2, Ud2 = Un * Un, Ud * Ud
        above = []
        for i in range(n):
            UNi = Un2 * N[i]
            row = G[i * n:(i + 1) * n]
            js = (np.nonzero(~inE[i, i + 1:])[0] + (i + 1)).tolist()
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
        U = U0
        if above:
            flint.ctx.prec = FX.APREC
            try:
                for i, j in above:
                    U = max(U, FX.arb_up(flint.arb(G[i * n + j]) / (nrm[i] * nrm[j])))
            finally:
                flint.ctx.prec = old
        taup = M_lo - U
    s = FX.sqrt_up(Fraction(d))
    delta = min(h / s, taup / 4)
    rec.update(h0=FX.dn_float(h0), dist_l1=FX.up_float(dist1), h=FX.dn_float(h), tau_prime=FX.dn_float(taup),
               delta=FX.dn_float(delta), delta_by="tau'" if taup / 4 <= h / s else "h")
    rec["status"] = "FLEX" if (h > 0 and taup > 0) else "FAIL"
    r8 = meta["record_v8"]
    same = {k: rec.get(k) == r8.get(k) for k in ("status", "rows", "columns", "beta", "C_norm_inf", "h", "tau_prime",
                                                  "delta", "delta_by")}
    return rec, same


def main():
    out = sys.argv[1]
    res = {"program_sha256": SHA, "files": {}}
    for f in sys.argv[2:]:
        t0, c0 = time.time(), time.process_time()
        z = np.load(f, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        arr = {k: z[k] for k in z.files if k != "meta"}
        h = hashlib.sha256(open(f, "rb").read()).hexdigest()
        if meta["kind"] == "positive":
            log, same = check_positive(meta, arr)
        else:
            log, same = check_negative(meta, arr)
        res["files"][os.path.basename(f)] = {"sha256": h, "kind": meta["kind"], "d": meta["d"], "n": meta["n"],
                                             "status": log.get("status"), "same_as_record": same,
                                             "all_same": all(same.values()), "result": log,
                                             "cpu_seconds": round(time.process_time() - c0, 2)}
        print(os.path.basename(f), log.get("status"), "same as record:", all(same.values()),
              [k for k, v in same.items() if not v], round(time.time() - t0, 1), "s", flush=True)
    json.dump(res, open(out, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
