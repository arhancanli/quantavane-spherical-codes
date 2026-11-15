"""Census version 9: the ranks of the strict entries above the size limit (paper_v8, Section 6.8) recomputed with
FLINT's nmod_mat, an established implementation, at a prime different from the one of run_above.py.

  nice -n 10 python3 -I flintrank9.py OUTDIR [--only d_n,...] [--max-mb 2500]

For every entry that paper_v8/above/run_above.py or above_more.py certified STRICT: the exact core is rebuilt by
certify/v5/core.py (released, unchanged; the field Q(sqrt D) recorded by run_above.py), the rigidity matrix M of
Proposition 3.19(d) is formed mod p by run_above.modp_csr (released, unchanged), with p the largest prime below 2^26 in
which D is a square (kfield.prime_for), and K = M^T M mod p is formed from scipy's sparse int64 product (entries of M
below 2^26, so every product is below 2^52 and every sum of at most 2^11 of them below 2^63: exact). The rank of K is
then computed by FLINT (python-flint nmod_mat.rank), which replaces the hand-written blocked LU of run_above.py.
K is dense in FLINT (8 bytes per entry), so an entry is skipped when that exceeds --max-mb.
"""
import hashlib, json, os, sys, time, traceback
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
from kfield import prime_for

SHA = {"flintrank9.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
       "run_above.py": RA.SHA}


def strict_records():
    out = {}
    for sub in ("out", "out_more"):
        dd = os.path.join(CERT, "paper_v8", "above", sub)
        for f in sorted(os.listdir(dd)):
            if f.endswith(".json"):
                j = json.load(open(os.path.join(dd, f)))
                if str(j.get("status", "")).startswith("STRICT"):
                    out[(j["d"], j["n"])] = j
    return out


def one(d, n, rec, max_mb, log):
    disc = int(rec.get("disc", 1))
    cols = int(rec["cols"])
    mb = cols * cols * 8 / 2 ** 20
    log.update(cols=cols, dense_mb=round(mb), v8_rank=rec.get("rank_mod_p"), v8_primes=[t["prime"] for t in rec.get("rank_tries", [])])
    if mb > max_mb:
        log["status"] = "skipped: the dense FLINT matrix exceeds %d MB" % max_mb
        return
    t0 = time.time()
    C = core.build(SNAP, d, n, disc)
    if C["status"] != "ok":
        log["status"] = "core: " + C["status"]
        return
    G, V, E, r, B = (C[k] for k in ("G", "V", "E", "r", "B"))
    nV = len(V)
    assert nV * d == cols and r == d
    # the support of the stress, if run_above.py used Theorem 2.4 on a subset E2 of the contacts
    if "support_contacts" in rec and str(rec.get("status", "")).startswith("STRICT (bar-rigid with a strictly positive stress on a subset"):
        log["status"] = "skipped: certified on a subset of the contacts"
        return
    log["seconds_core"] = round(time.time() - t0, 1)
    p, rt = prime_for(disc, start=1 << 26)
    while True:
        try:
            M = RA.modp_csr(G, B, E, d, nV, p, rt)
            break
        except ValueError:
            p, rt = prime_for(disc, start=p)
    assert M.data.min() >= 0 and M.data.max() < p < 2 ** 26
    colnnz = int(np.diff(M.tocsc().indptr).max())
    assert colnnz <= 2 ** 11
    t1 = time.time()
    Kc = (M.T.tocsr() @ M.tocsc()).tocoo()               # exact in int64 (see the docstring)
    del M
    A = flint.nmod_mat(cols, cols, p)
    vals = np.mod(Kc.data, p)
    for i, j, v in zip(Kc.row.tolist(), Kc.col.tolist(), vals.tolist()):
        if v:
            A[i, j] = v
    nnz = len(Kc.data)
    del Kc, vals
    log["seconds_build"] = round(time.time() - t1, 1)
    t2 = time.time()
    rk = A.rank()
    del A
    g = d * (d - 1) // 2
    log.update(prime=p, sqrt_D_mod_p=rt, K_nonzeros=nnz, max_column_nonzeros=colnnz, flint_rank=int(rk),
               bound=cols - g, seconds_flint_rank=round(time.time() - t2, 1))
    log["status"] = "AGREE" if rk == cols - g and rk == rec.get("rank_mod_p") else "DISAGREE"


def main():
    out = sys.argv[1]
    args = sys.argv[2:]
    opt = lambda k, dflt: args[args.index(k) + 1] if k in args else dflt
    only = set(opt("--only", "").split(",")) if "--only" in args else None
    max_mb = int(opt("--max-mb", "2500"))
    recs = strict_records()
    os.makedirs(out, exist_ok=True)
    for (d, n), rec in sorted(recs.items(), key=lambda kv: kv[1]["cols"]):
        if only and f"{d}_{n}" not in only:
            continue
        js = os.path.join(out, f"d{d}_n{n}.json")
        if os.path.exists(js):
            continue
        log = {"d": d, "n": n, "program_sha256": SHA}
        t0 = time.time()
        try:
            one(d, n, rec, max_mb, log)
        except Exception:
            log["status"] = "error"
            log["error"] = traceback.format_exc()[-1500:]
        log["wall"] = round(time.time() - t0, 1)
        with open(js + ".part", "w") as fo:
            json.dump(log, fo)
        os.replace(js + ".part", js)
        print(f"d{d}_n{n} cols={log.get('cols')} {log['status']} rank={log.get('flint_rank')} {log['wall']}s", flush=True)


if __name__ == "__main__":
    main()
