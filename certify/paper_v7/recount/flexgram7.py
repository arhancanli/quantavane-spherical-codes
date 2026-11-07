"""The separation tau' of the negative certificate (Theorem 4.1), with no dense BLAS product in any bound.

  nice -n 10 python3 -I flexgram.py OUT.jsonl [--shard I/K]

For every entry with a FLEX certificate, the contact set E is rebuilt exactly as cq.certify_flex built it (the double
precision Gram matrix by NumPy's matrix product in blocks of 256 rows, E = the pairs with G >= M - tau; untrusted, it
only proposes), and checked against the census record (number of pairs and of points in pairs). Then, TRUSTED: every
inner product of the double-precision code X is recomputed as a sum of individual products (elementwise products and
NumPy sums, so the error is at most gamma_d sum |x_ik x_jk| <= gamma_d (1 + 2^-48) in any order), and with the distance
e between X and the exactly normalised stored code (stored_norm_dist_l1 of the census record, an l1 and hence l2 bound
per point) every inner product of the stored code lies within eps = gamma_d (1 + 2^-48) + 2 e (1 + 2^-48) + d 2^-1074 of
the recomputed one. So M(stored) >= M_lo = max over E of (G - eps), every pair outside E has inner product at most
U_out = max over pairs outside E of (G + eps), and every pair with inner product above M(stored) - tau' lies in E for
tau' = M_lo - U_out (rounded down). Theorem 4.1 then holds with delta = min(h / sqrt d, tau' / 4), h the census's
flex_radius_l1_stored.

Version 7: the bound sum_k |x_ik x_jk| <= |x_i| |x_j| <= max_i |x_i|^2 uses the largest squared norm of the double
precision points, computed exactly (rational sum of squares of doubles), instead of an assumed 1 + 2^-48; and the
SHA-256 of the sorted pair list E is recorded, so that the negative-certificate check of version 7 (flexexact.py) can
show that its h and this tau' refer to the same E.
"""
import json, math, os, sys
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.path.join(CERT, "census"))
import numpy as np
import snapshot

U53 = 2.0 ** -53


def check(X, d, n, a):
    tol = a["tol"]
    # the proposal, exactly as certify_flex (untrusted)
    G = np.empty((n, n))
    for i0 in range(0, n, 256):
        G[i0:i0 + 256] = X[i0:i0 + 256] @ X.T
    np.fill_diagonal(G, -2.0)
    M = float(G.max())
    inE = np.triu(G >= M - tol, 1)
    m = int(inE.sum())
    I_, J_ = np.nonzero(inE)
    pts = len(np.unique(np.concatenate([I_, J_])))
    import hashlib
    from fractions import Fraction
    Esha = hashlib.sha256(json.dumps([[int(i), int(j)] for i, j in zip(I_, J_)]).encode()).hexdigest()
    rec = {"contacts": m, "points_in_contacts": pts, "E_sha256": Esha,
           "proposal_matches_census": m == a["contacts"] and pts == a["points_in_contacts"]}
    # trusted bounds: individual products, NumPy sums; |x_i| |x_j| <= N2 = max_i |x_i|^2 (exact) when N2 >= 1, else <= 1
    N2 = max(sum(Fraction(float(v)) ** 2 for v in row) for row in X)
    N2f = max(1.0, float(N2) * (1 + 2.0 ** -50))      # >= max(1, N2): float(N2) has relative error <= 2^-53
    gd = d * U53 / (1 - d * U53)
    e_st = a["stored_norm_dist_l1"]
    # |<x_i,x_j> - <xh_i,xh_j>| <= |x_i - xh_i| |x_j| + |xh_i| |x_j - xh_j| <= e_st (sqrt(N2f) + 1) <= e_st (N2f + 1)
    eps = gd * N2f + e_st * (N2f + 1) + d * 2.0 ** -1074
    eps = math.nextafter(eps * (1 + 2.0 ** -40), math.inf)
    rec["N2_max_minus_1"] = float(N2 - 1)
    m_lo = -math.inf
    u_out = -math.inf
    for i in range(n):
        g = (X[i][None, :] * X).sum(1)          # elementwise products, then sums
        g[i] = -2.0
        e_row = inE[i] | inE[:, i]
        if e_row.any():
            m_lo = max(m_lo, float(g[e_row].max()))
        out = ~e_row
        out[i] = False
        if out.any():
            u_out = max(u_out, float(g[out].max()))
    m_lo = math.nextafter(m_lo - eps, -math.inf)
    u_out = math.nextafter(u_out + eps, math.inf)
    taup = math.nextafter(m_lo - u_out, -math.inf) if u_out > -math.inf else 2.0
    h = a["flex_radius_l1_stored"]
    delta = math.nextafter(min(h / math.sqrt(d), taup / 4) * (1 - 2.0 ** -40), 0)
    rec.update(eps=eps, M_lo=m_lo, U_out=u_out, tau_prime=taup, h=h, delta=delta, tol=tol)
    return rec


def main():
    outp = sys.argv[1]
    args = sys.argv[2:]
    sh, nsh = 0, 1
    if "--shard" in args:
        sh, nsh = (int(x) for x in args[args.index("--shard") + 1].split("/"))
    s5 = json.load(open(os.path.join(CERT, "v5", "out", "summary_v5.json")))
    keys = sorted(tuple(map(int, k.split(","))) for k, v in s5["class"].items() if v[1] == "flex")
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    done = set()
    if os.path.exists(outp):
        done = {(r["d"], r["n"]) for r in map(json.loads, open(outp))}
    for k, (d, n) in enumerate(keys):
        if k % nsh != sh or (d, n) in done:
            continue
        r = json.load(open(os.path.join(CERT, "census", "results_v4", f"d{d}_n{n}.json")))
        a = [x for x in r["attempts"] if x.get("status") == "FLEX"][0]
        pts = snapshot.read_points(SNAP, idx[(d, n)])
        X = np.array([[float(x) for x in p] for p in pts])
        X = X / np.linalg.norm(X, axis=1, keepdims=True)
        rec = {"d": d, "n": n}
        rec.update(check(X, d, n, a))
        with open(outp, "a") as f:
            f.write(json.dumps(rec) + "\n")
        if k % 500 == 0:
            print(d, n, rec["proposal_matches_census"], rec["tau_prime"], rec["delta"], flush=True)


if __name__ == "__main__":
    main()
