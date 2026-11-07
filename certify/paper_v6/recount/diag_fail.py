"""UNTRUSTED diagnostics (numerical, double precision) of the attempted entries whose first certificate attempt failed
the existence test (paper, Section on limits; the referee's item M3).

  nice -n 10 python3 -I diag_fail.py OUT.jsonl [--shard I/K]

For each entry without a verdict whose first isostatic or hyperstatic attempt (in the order of the census rule) ended in
'FAIL existence (Krawczyk)': the proposal (V, E) at that attempt's tolerance is rebuilt exactly as the census built it
(sphcert.contact_structure; for a hyperstatic graph the square subset E' of cq.choose_square_subset_qr), and at the
stored code (the slice reference) we record
  - the singular values of the square Jacobian DF (smallest, largest, and how many are below 1e-8 and 1e-12 times the
    largest): a numerically singular Jacobian is what Proposition 3.2(a) predicts for a flexing contact framework;
  - the separation of the proposal: the spread of the proposed contacts (largest cosine minus the smallest cosine in E)
    and the gap to the next pair (largest cosine minus the largest cosine of a pair of V outside E), both against tau;
  - from the census record of the attempt: the Newton residuals (did the chord iteration converge?) and the largest
    row sum of |I - C J_0| + |C| Delta in the Krawczyk test (the contraction factor; it must be below 1).
Nothing here is a proof; it classifies why the existence test failed.
"""
import json, os, sys, time
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, CERT)
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import cq
import sphcert as S
import snapshot


def entries():
    s5 = json.load(open(os.path.join(CERT, "v5", "out", "summary_v5.json")))
    out = []
    for key, v in s5["class"].items():
        d, n = map(int, key.split(","))
        r = json.load(open(os.path.join(CERT, "census", "results_v4", f"d{d}_n{n}.json")))
        att = [a for a in r["attempts"] if a.get("kind") in ("isostatic", "hyperstatic")]
        if not att:
            continue
        if v[0] in ("strict", "local", "nonstrict", "notlocal") and v[1] != "strict_second":
            continue
        if att[0]["status"] != "FAIL existence (Krawczyk)":
            continue
        out.append((att[0]["unknowns"], d, n, att[0]))
    return sorted(out)


def diag(X, d, n, a):
    tol = a["tol"]
    V, E, M, G = S.contact_structure(X, d, tol)
    nV, m = len(V), len(E)
    iso = nV * (d - 1) - d * (d - 1) // 2 + 1
    xs = X[V]
    rec = {"nV": int(nV), "contacts": int(m), "iso": int(iso), "kind": a["kind"]}
    if m > iso:
        full = S.System(xs, E, d)
        sub = cq.choose_square_subset_qr(full, xs, d, nV, iso)
        Esq = [E[e] for e in sub]
    else:
        Esq = E
    sysm = S.System(xs, Esq, d)
    J = sysm.jac(xs).toarray()
    sv = np.linalg.svd(J, compute_uv=False)
    rec.update(U=int(J.shape[0]), sv_min=float(sv.min()), sv_max=float(sv.max()),
               rel_sv_min=float(sv.min() / sv.max()),
               n_below_1e8=int(np.sum(sv < 1e-8 * sv.max())), n_below_1e12=int(np.sum(sv < 1e-12 * sv.max())))
    GV = G[np.ix_(V, V)]
    Eset = set(E)
    cos_E = np.array([GV[a_, b_] for a_, b_ in E])
    iu, ju = np.triu_indices(nV, 1)
    mask = np.array([(int(i), int(j)) not in Eset for i, j in zip(iu, ju)])
    nxt = float(GV[iu[mask], ju[mask]].max()) if mask.any() else None
    rec.update(tol=tol, spread_E=float(M - cos_E.min()), gap_next=(float(M - nxt) if nxt is not None else None))
    nr = a.get("newton_residuals") or []
    rec.update(newton_first=nr[0] if nr else None, newton_last=nr[-1] if nr else None,
               krawczyk_rowsum=(a.get("krawczyk_attempts") or [{}])[0].get("max_rowsum"))
    return rec


def main():
    outp = sys.argv[1]
    args = sys.argv[2:]
    sh, nsh = 0, 1
    if "--shard" in args:
        sh, nsh = (int(x) for x in args[args.index("--shard") + 1].split("/"))
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    done = set()
    if os.path.exists(outp):
        done = {(j["d"], j["n"]) for j in map(json.loads, open(outp))}
    for k, (U, d, n, a) in enumerate(entries()):
        if k % nsh != sh or (d, n) in done:
            continue
        t0 = time.time()
        rec = {"d": d, "n": n}
        try:
            pts = snapshot.read_points(SNAP, idx[(d, n)])
            X = np.array([[float(x) for x in p] for p in pts])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            rec.update(diag(X, d, n, a))
        except Exception:
            import traceback
            rec["error"] = traceback.format_exc()[-300:]
        rec["seconds"] = round(time.time() - t0, 1)
        with open(outp, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(d, n, U, rec.get("rel_sv_min"), rec["seconds"], flush=True)


if __name__ == "__main__":
    main()
