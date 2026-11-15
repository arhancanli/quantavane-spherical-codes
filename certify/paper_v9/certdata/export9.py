"""Census version 9, item B: the certificate data of a representative subset of the census, written as files that the
stand-alone verifier verify_cert9.py decides without the census kernel (no LAPACK, no floating-point step).

  nice -n 10 python3 -I export9.py OUTDIR --pos d_n:STATUS,... --neg d_n,...

Positive certificate (Theorems 3.5 and 3.11; STATUS one of PASS, LOCALMIN, NONOPT): the census certificate is re-run
exactly as exactcheck.py does (cq.certify with keep; UNTRUSTED) and the data that exactcheck.check_entry reads are
written exactly: the stored decimals of the code, the tolerance, the non-rattlers V, the contacts E, the square subset
E' and the extra contacts, the centre z~ = (Z, C) / 2^100 (integers, as decimal strings), and, as binary64 arrays (exact
rationals), the approximate inverse R, the slice reference xi, the coefficients lambda_f and the rattler vectors Y.
Negative certificate (Theorem 4.1): the proposal E (pairs of the double-precision code within tau of its largest
cosine), the reference X0 (integers over 2^52) and the right inverse rounded to integers over 2^t (exactly the matrix
whose exact test flexexact8.py passed: normal equations, QR, or for (17,317) the SVD of flexrefine8.py).
Files: OUTDIR/{pos,neg}_d{d}_n{n}.npz (integer arrays) with the metadata in the same file (json string).
"""
import hashlib, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.path.join(CERT, "paper_v8", "verify"))
import exactcheck as XC                       # sets the library paths
sys.path.insert(0, os.path.join(CERT, "paper_v8", "recount"))
import flexexact8 as FX
import numpy as np
import cq
import snapshot

SHA = {"export9.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
       "exactcheck.py": XC.SHA, "flexexact8.py": FX.SHA}


def code(d, n):
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    pts = snapshot.read_points(SNAP, idx[(d, n)])
    dec = [[t.strip() for t in p] for p in pts]
    X = np.array([[float(x) for x in p] for p in dec])
    return dec, X / np.linalg.norm(X, axis=1, keepdims=True)


def export_pos(d, n, status, out):
    import glob
    rec = None
    for sub in ("out_exact2", "out_exact"):
        fs = sorted(glob.glob(os.path.join(CERT, "paper_v8", "verify", sub, f"d{d}_n{n}_{status}_*.json")))
        if fs:
            rec = json.load(open(fs[0]))
            break
    assert rec is not None and rec["status"] == status, (d, n)
    dec, X = code(d, n)
    keep = {}
    cl = cq.certify(X, d, n, rec["tol"], dec=dec, keep=keep)
    assert cl.get("status") == status, "census certificate not reproduced on re-run"
    meta = {"kind": "positive", "d": d, "n": n, "status": status, "tol": rec["tol"], "K": int(keep["K"]),
            "C": str(int(keep["C"])), "Z": [[str(int(z)) for z in row] for row in keep["Z"]],
            "V": [int(v) for v in keep["V"]], "E": [[int(a), int(b)] for a, b in keep["E"]],
            "Esq": [[int(a), int(b)] for a, b in keep["Esq"]], "extra": [[int(a), int(b)] for a, b in keep["extra"]],
            "dec": dec, "record_v8": rec, "record_dir": sub, "program_sha256": SHA}
    beta = np.asarray(keep["beta_f"], dtype=np.float64) if keep["extra"] else np.zeros(0)
    np.savez_compressed(os.path.join(out, f"pos_d{d}_n{n}.npz"), meta=np.array(json.dumps(meta)),
                        R=np.asarray(keep["R"], dtype=np.float64), xs=np.asarray(keep["xs"], dtype=np.float64),
                        Y=np.asarray(keep["Y"], dtype=np.float64), beta_f=beta)
    return {"U": int(len(keep["V"]) * d + 1), "bytes": os.path.getsize(os.path.join(out, f"pos_d{d}_n{n}.npz"))}


def export_neg(d, n, out):
    rec = None
    for sub in ("out_flex8r", "out_flex8"):
        p = os.path.join(CERT, "paper_v8", "recount", sub, f"d{d}_n{n}.json")
        if os.path.exists(p) and json.load(open(p)).get("status") == "FLEX":
            rec, which = json.load(open(p)), sub
            break
    assert rec is not None
    dec, X = code(d, n)
    tol = rec["tol"]
    I_, J_, _ = FX.proposal(X, n, tol)
    pts = np.unique(np.concatenate([I_, J_]))
    loc = -np.ones(n, int)
    loc[pts] = np.arange(len(pts))
    Xi = np.rint(np.ldexp(X, FX.GRID)).astype(np.int64)
    Li, Lf = FX.build_L(Xi, X, pts, loc, I_, J_, d)
    UL, UC = Lf.shape
    how = rec["right_inverse"]
    if how == "normal equations":
        R = np.linalg.inv((Lf @ Lf.T).toarray())
        Cf = np.asarray(Lf.T.tocsr() @ R)
    elif how == "QR":
        import scipy.linalg as sl
        Qf, Rf = sl.qr(Lf.T.toarray(), mode="economic", check_finite=False)
        Cf = sl.solve_triangular(Rf, Qf.T, lower=False, check_finite=False).T
    else:
        raise ValueError("the SVD right inverse of flexrefine8.py is not exported here")
    mx = float(np.max(np.abs(Cf)))
    t = 52 - __import__("math").frexp(mx)[1]
    Cint = np.rint(np.ldexp(Cf, t)).astype(np.int64)
    meta = {"kind": "negative", "d": d, "n": n, "tol": tol, "t": t, "right_inverse": how, "dec": dec,
            "pairs": [[int(i), int(j)] for i, j in zip(I_, J_)], "record_v8": rec, "record_dir": which,
            "program_sha256": SHA}
    np.savez_compressed(os.path.join(out, f"neg_d{d}_n{n}.npz"), meta=np.array(json.dumps(meta)), X0=Xi, Cint=Cint)
    return {"rows": int(UL), "cols": int(UC), "bytes": os.path.getsize(os.path.join(out, f"neg_d{d}_n{n}.npz"))}


def main():
    out = sys.argv[1]
    args = sys.argv[2:]
    opt = lambda k, dflt: args[args.index(k) + 1] if k in args else dflt
    os.makedirs(out, exist_ok=True)
    log = {}
    for item in filter(None, opt("--pos", "").split(",")):
        dn, st = item.split(":")
        d, n = map(int, dn.split("_"))
        t0 = time.time()
        log[f"pos_{dn}"] = dict(export_pos(d, n, st, out), seconds=round(time.time() - t0, 1))
        print("pos", dn, log[f"pos_{dn}"], flush=True)
    for dn in filter(None, opt("--neg", "").split(",")):
        d, n = map(int, dn.split("_"))
        t0 = time.time()
        log[f"neg_{dn}"] = dict(export_neg(d, n, out), seconds=round(time.time() - t0, 1))
        print("neg", dn, log[f"neg_{dn}"], flush=True)
    json.dump({"program_sha256": SHA, "files": log}, open(os.path.join(out, "export_log.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
