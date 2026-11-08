"""How much a sharper rigidity constant would enlarge the certified neighbourhoods (census version 8).
UNTRUSTED, double precision: nothing here is a certificate.

  nice -n 10 python3 -I zetasample.py OUT.jsonl [--max-unknowns 1600]

For every isostatic certificate (census PASS) with at most --max-unknowns unknowns: the stored code X_V normalised,
the certified contacts E, the tangent space T and the complement W of the infinitesimal rotations in T; the smallest
singular value s of l restricted to W (2-norm); and zeta2 = s / sqrt(|E|), which bounds ||l(t)||_inf / ||t||_inf from
below as zeta of Lemma 3.8 does. Reported with the certified zeta, rho and distance of exactcheck.py, and the radius
rho2 = rho * zeta2 / zeta that the same proof would give with zeta2.
"""
import json, math, os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "verify"))
import exactcheck as XC
import numpy as np
import snapshot


def main():
    outp = sys.argv[1]
    umax = int(sys.argv[sys.argv.index("--max-unknowns") + 1]) if "--max-unknowns" in sys.argv else 1600
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    done = set()
    if os.path.exists(outp):
        done = {(r["d"], r["n"]) for r in map(json.loads, open(outp))}
    ex = {}
    for p in glob.glob(os.path.join(os.path.dirname(HERE), "verify", "out_exact*", "*_PASS_*.json")):
        j = json.load(open(p))
        if j["status"] == "PASS":
            ex[(j["d"], j["n"])] = j
    for (d, n), j in sorted(ex.items(), key=lambda t: t[1]["unknowns"]):
        if (d, n) in done or j["unknowns"] > umax:
            continue
        pts = snapshot.read_points(XC.SNAP, idx[(d, n)])
        X = np.array([[float(x) for x in p] for p in pts])
        X /= np.linalg.norm(X, axis=1, keepdims=True)
        G = X @ X.T
        np.fill_diagonal(G, -2.0)
        M = G.max()
        A = G >= M - j["tol"]
        alive = np.ones(n, bool)
        while True:
            deg = (A & alive[None, :] & alive[:, None]).sum(1)
            drop = alive & (deg < d)
            if not drop.any():
                break
            alive &= ~drop
        V = np.flatnonzero(alive)
        XV = X[V]
        nV = len(V)
        E = [(a, b) for a in range(nV) for b in range(a + 1, nV) if A[V[a], V[b]]]
        # an orthonormal basis of W: tangent vectors orthogonal to the rotation fields
        T = []
        for a in range(nV):
            Q, _ = np.linalg.qr(np.c_[XV[a], np.eye(d)])
            for k in range(1, d):
                v = np.zeros(nV * d)
                v[a * d:(a + 1) * d] = Q[:, k]
                T.append(v)
        T = np.array(T).T                                          # (nV d) x (nV (d-1))
        rots = []
        for s_ in range(d):
            for t_ in range(s_ + 1, d):
                Ast = np.zeros((d, d))
                Ast[s_, t_], Ast[t_, s_] = 1.0, -1.0
                rots.append((XV @ Ast.T).ravel())
        Rm = T.T @ np.array(rots).T                                 # rotation fields in the basis of T
        Ur, _, _ = np.linalg.svd(Rm, full_matrices=True)
        Wb = T @ Ur[:, len(rots):]                                  # orthonormal basis of W
        Lm = np.zeros((len(E), nV * d))
        for e, (a, b) in enumerate(E):
            Lm[e, a * d:(a + 1) * d] = XV[b]
            Lm[e, b * d:(b + 1) * d] = XV[a]
        s = np.linalg.svd(Lm @ Wb, compute_uv=False)
        smin = float(s[-1]) if len(s) else 0.0
        zeta2 = smin / math.sqrt(len(E))
        rec = {"d": d, "n": n, "unknowns": j["unknowns"], "contacts": len(E), "nV": int(nV),
               "same_structure": bool(nV == j.get("nV") and len(E) == j.get("contacts")), "sigma_min_W": smin, "zeta2": zeta2,
               "zeta": j.get("zeta"), "rho": j.get("rho"), "dist_ls": j.get("dist_ls"),
               "rho2": (j["rho"] * zeta2 / j["zeta"]) if j.get("zeta") and j.get("rho") else None}
        with open(outp, "a") as f:
            f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
