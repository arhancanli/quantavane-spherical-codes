"""Census version 9, item A1: recognition at higher degree for the hyperstatic entries whose strictness is open
(UNTRUSTED survey; strict9.py decides).

  nice -n 10 python3 -I recog9.py OUTDIR [--only d_n,...] [--shard I/K] [--bits 2000] [--maxdeg 30] [--maxvals 20000]

For every entry of class 'local' in paper_v8/v8_verdicts.tsv (a local minimum certified by the hyperstatic theorem,
strictness open):
 1. the census certificate is re-run (cq.certify, the census tolerance) for the zero P' of the square system of E';
 2. P' is refined by chord steps with exact residuals to the grid 2^-BITS (rec9lib.refine_hp);
 3. the largest deviation of an extra contact from c' at the refined P' is recorded (numerical evidence of equality);
 4. the distinct inner products of P' are clustered;
 5. LLL looks for the minimal polynomial of c' (degree <= MAXDEG) and then for every distinct inner product as an
    element of Q(c'); if some inner product is not in Q(c'), a primitive element theta (an inner product of largest
    degree, or c' plus it) is tried, with c' and every inner product recognised in Q(theta).
Output: OUTDIR/d{d}_n{n}.json with the outcome. Nothing here is trusted.
"""
import hashlib, json, math, os, sys, time, traceback, random

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rec9lib as L                     # sets the library paths
import numpy as np
import cq
import snapshot

SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
CERT = L.CERT
SHA = {"recog9.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
       "rec9lib.py": hashlib.sha256(open(os.path.join(HERE, "rec9lib.py"), "rb").read()).hexdigest()}
SHA.update(cq.CODE_SHA256)


def open_entries(cls="local"):
    """The (d, n) of class cls in the released v8 per-entry file: a verdict ('local') or a certificate
    ('strict_hyper': strict by recognition in version 8)."""
    path = os.path.join(CERT, "paper_v8", "v8_verdicts.tsv")
    rows = [l.rstrip("\n").split("\t") for l in open(path) if not l.startswith("#")]
    head = rows[0]
    out = []
    for r in rows[1:]:
        r = dict(zip(head, r))
        if r["verdict"] == cls or r["certificate"] == cls:
            out.append((int(r["d"]), int(r["n"])))
    return out


def census_tol(d, n):
    j = json.load(open(os.path.join(CERT, "census", "results_v4", f"d{d}_n{n}.json")))
    a = [x for x in j["attempts"] if x.get("status") == "LOCALMIN"][0]
    return a["tol"], a["unknowns"]


def survey(d, n, X, dec, tol, bits, maxdeg, maxvals, log, out=None, probe=False):
    """The survey of one entry; with a dict `out`, also hands the proposal to the exact stage (strict9.py): the refined
    point, the clusters of inner products and their recognised coordinates, the field and the census data."""
    t0 = time.time()
    keep = {}
    cl = cq.certify(X, d, n, tol, dec=dec, keep=keep)
    log["census_status_rerun"] = cl.get("status")
    log["seconds_certify"] = round(time.time() - t0, 1)
    if cl.get("status") != "LOCALMIN":
        log["status"] = "not LOCALMIN on re-run"
        return log
    V, E, Esq, extra = keep["V"], keep["E"], keep["Esq"], keep["extra"]
    nV = len(V)
    log.update(nV=int(nV), contacts=len(E), extra=len(extra), rattlers=int(n - nV))
    sysm, R = keep["sysm"], keep["R"]
    t1 = time.time()
    Z, C, lg = L.refine_hp(sysm, keep["Z"], keep["C"], keep["K"], np.asarray(R, dtype=np.float64), bits, log=log)
    rn = float(np.abs(np.asarray(R, dtype=np.float64)).sum(axis=1).max())
    if out is not None:
        out.update(V=V, E=E, Esq=Esq, extra=extra, Z=Z, C=C, K=bits)
    del keep, R
    acc = int(-lg - math.log2(rn) - 8)
    acc = min(acc, bits - 8)
    log.update(seconds_refine=round(time.time() - t1, 1), accuracy_bits=acc, R_norm_inf=rn)
    K = bits
    # extra contacts against c' at the refined point
    CK = C << K
    dev = 0
    for (a, b) in E:
        dev = max(dev, abs(sum(x * y for x, y in zip(Z[a], Z[b])) - CK))
    log["extra_dev_log2"] = (dev.bit_length() - 2 * K) if dev else -2 * K
    # distinct inner products of P' (clusters at acc - 40 bits)
    cb = max(acc - 40, 64)
    shift = 2 * K - cb
    half = 1 << (shift - 1)
    clusters = {}
    for a in range(nV):
        za = Z[a]
        for b in range(a + 1, nV):
            g = sum(x * y for x, y in zip(za, Z[b]))
            key = (g + half) >> shift
            if key not in clusters:
                clusters[key] = (g, a, b, 1)
            else:
                g0, a0, b0, cnt = clusters[key]
                clusters[key] = (g0, a0, b0, cnt + 1)
    log["distinct_inner_products"] = len(clusters)
    log["c_digits"] = str(C * 10 ** 60 // (1 << K))
    P = acc - 16
    t2 = time.time()
    mp = L.minpoly_lll(C, K, P, maxdeg, probe=probe)
    log["seconds_minpoly"] = round(time.time() - t2, 2)
    if mp is None:
        log["status"] = "no polynomial of degree <= %d" % maxdeg
        return log
    cpoly, cbits = mp
    log.update(c_poly=cpoly, c_degree=len(cpoly) - 1, c_poly_bits=round(cbits, 1))
    if len(clusters) > maxvals:
        log["status"] = "too many distinct inner products"
        return log
    vals = [(g, a, b) for (g, a, b, cnt) in clusters.values()]
    # try Q(c') first
    t3 = time.time()

    def try_field(theta_num, theta_den_bits, D, items):
        tp = L.powers_int(theta_num, theta_den_bits, D - 1, P)
        recs, worst = [], 0.0
        for (gnum, gden) in items:
            co = L.in_field_lll(gnum, gden, tp, P)
            if co is None:
                return None, (gnum, gden)
            recs.append(co)
            worst = max(worst, max(math.log2(abs(p) + 1) + math.log2(abs(q) + 1) for p, q in co))
        return recs, worst

    D = len(cpoly) - 1
    recs, info = try_field(C, K, D, [(g, 2 * K) for g, a, b in vals])
    if recs is not None:
        log.update(status="recognised in Q(c)", field="Q(c)", field_degree=D, field_poly=cpoly,
                   coeff_bits_max=round(info, 1), seconds_values=round(time.time() - t3, 1))
        if out is not None:
            out.update(field_poly=cpoly, theta=(C, K), c_elem=[(0, 1), (1, 1)] + [(0, 1)] * (D - 2) if D >= 2 else None,
                       clusters=[(a, b, rc) for (g, a, b), rc in zip(vals, recs)], cluster_shift=shift,
                       cluster_half=half)
            if D == 1:
                # Q(c) = Q: c itself is the rational root of the linear polynomial
                out["c_elem"] = [(-cpoly[0], cpoly[1])]
        return log
    # an inner product outside Q(c'): look for a primitive element among the inner products
    log["first_value_outside_Qc"] = str(info[0] * 10 ** 40 // (1 << info[1]))
    cand = []
    degs = {}
    for g, a, b in vals[:200]:
        m2 = L.minpoly_lll(g, 2 * K, P, maxdeg, probe=probe)
        degs[(a, b)] = None if m2 is None else len(m2[0]) - 1
        if m2 is not None:
            cand.append((len(m2[0]) - 1, g, m2[0]))
    log["value_degrees"] = sorted(set(v for v in degs.values() if v is not None))
    log["values_without_poly"] = sum(1 for v in degs.values() if v is None)
    cand.sort(key=lambda t: -t[0])
    tried = 0
    for dg, g, gp in cand[:3]:
        for mode in ("g", "c+g"):
            if mode == "g":
                tnum, tden, tpoly = g, 2 * K, gp
            else:
                tnum, tden = (C << K) + g, 2 * K
                m3 = L.minpoly_lll(tnum, tden, P, maxdeg, probe=probe)
                if m3 is None:
                    continue
                tpoly = m3[0]
            Dt = len(tpoly) - 1
            tried += 1
            items = [((C << K), 2 * K)] + [(g2, 2 * K) for g2, a2, b2 in vals]
            recs, info = try_field(tnum, tden, Dt, items)
            if recs is not None:
                log.update(status="recognised in Q(theta)", field="Q(theta), theta = " +
                           ("an inner product" if mode == "g" else "c + an inner product"),
                           field_degree=Dt, field_poly=tpoly, coeff_bits_max=round(info, 1),
                           seconds_values=round(time.time() - t3, 1))
                if out is not None:
                    out.update(field_poly=tpoly, theta=(tnum, tden), c_elem=recs[0],
                               clusters=[(a, b, rc) for (g, a, b), rc in zip(vals, recs[1:])], cluster_shift=shift,
                               cluster_half=half)
                return log
    log["primitive_elements_tried"] = tried
    log["status"] = "an inner product not recognised"
    log["seconds_values"] = round(time.time() - t3, 1)
    return log


def main():
    out = sys.argv[1]
    args = sys.argv[2:]
    opt = lambda k, dflt: args[args.index(k) + 1] if k in args else dflt
    sh, nsh = (int(x) for x in opt("--shard", "0/1").split("/"))
    bits, maxdeg, maxvals = int(opt("--bits", "2000")), int(opt("--maxdeg", "30")), int(opt("--maxvals", "20000"))
    only = set(opt("--only", "").split(",")) if "--only" in args else None
    probe = "--probe" in args
    pass2 = opt("--pass2", None)          # a first-pass output folder: only its unrecognised entries with few values
    maxfew = int(opt("--fewvals", "200"))
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    todo = []
    for d, n in open_entries():
        if only and f"{d}_{n}" not in only:
            continue
        if pass2:
            j1 = json.load(open(os.path.join(pass2, f"d{d}_n{n}.json")))
            if str(j1.get("status", "")).startswith("recognised") or (j1.get("distinct_inner_products") or 10 ** 9) > maxfew:
                continue
        tol, U = census_tol(d, n)
        todo.append((U, d, n, tol))
    todo.sort()
    os.makedirs(out, exist_ok=True)
    for k, (U, d, n, tol) in enumerate(todo):
        if k % nsh != sh:
            continue
        js = os.path.join(out, f"d{d}_n{n}.json")
        if os.path.exists(js):
            continue
        log = {"d": d, "n": n, "tol": tol, "unknowns": U, "bits": bits, "maxdeg": maxdeg, "probe": probe,
               "program_sha256": SHA}
        t0, c0 = time.time(), time.process_time()
        try:
            pts = snapshot.read_points(SNAP, idx[(d, n)])
            dec = [[t.strip() for t in p] for p in pts]
            X = np.array([[float(x) for x in p] for p in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            survey(d, n, X, dec, tol, bits, maxdeg, maxvals, log, probe=probe)
        except Exception:
            log["status"] = "error"
            log["error"] = traceback.format_exc()[-1500:]
        log["wall"] = round(time.time() - t0, 1)
        log["cpu_seconds"] = round(time.process_time() - c0, 1)
        with open(js + ".part", "w") as fo:
            json.dump(log, fo)
        os.replace(js + ".part", js)
        print(f"d{d}_n{n} U={U} {log['status']} deg={log.get('c_degree')} vals={log.get('distinct_inner_products')} "
              f"dev={log.get('extra_dev_log2')} acc={log.get('accuracy_bits')} {log['wall']}s", flush=True)


if __name__ == "__main__":
    main()
