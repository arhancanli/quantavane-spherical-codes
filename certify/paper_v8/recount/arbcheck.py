"""An independent re-check of census existence proofs with Arb ball arithmetic (python-flint), for a reproducible sample
(the referee's item m15: a cross-check with an established third-party verified arithmetic).

  nice -n 10 python3 -I arbcheck.py OUT.jsonl [--max-unknowns 400] [--every 20]

For each sampled entry certified PASS (isostatic) or LOCALMIN (hyperstatic; the square subsystem of E'), the census
certificate is re-run only to obtain its untrusted data: the proposal (V, E or E'), the slice reference xi and the grid
point z~ = (Z / 2^100, C / 2^100). The square system F of Section 3.1 is then evaluated in Arb, without any part of the
census kernel (sphcert.py's rounding bounds are not used): F(z~) as balls from exact rationals, the Jacobian as balls
over the whole box z~ + [-r, r]^U, an approximate inverse R from NumPy (untrusted), and Krawczyk's test
-R F(z~) + (I - R J(z~ + X)) X in int(X) (Theorem 3.3), all products in Arb matrices. A pass proves again that F has a
unique zero in the box; the enclosure of the common cosine is compared with the census's 30-digit interval.
"""
import json, os, sys, time
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, CERT)
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import flint
import cq
import snapshot

PREC = 200


def arb_check(keep, d):
    sysm, Z, C, K = keep["sysm"], keep["Z"], keep["C"], keep["K"]
    n, m = sysm.n, sysm.m
    E = list(zip(sysm.Ea.tolist(), sysm.Eb.tolist()))
    ST = sysm.ST
    xs = sysm.xs
    U = n * d + 1
    two = flint.arb(2) ** (-K)
    z = [flint.arb(Z[a][k]) * two for a in range(n) for k in range(d)] + [flint.arb(C) * two]
    xi = [[flint.arb(float(xs[a, k])) for k in range(d)] for a in range(n)]      # exact doubles

    def F(zz):
        out = []
        for a in range(n):
            out.append(sum((zz[a * d + k] ** 2 for k in range(d)), flint.arb(0)) - 1)
        for a, b in E:
            out.append(sum((zz[a * d + k] * zz[b * d + k] for k in range(d)), flint.arb(0)) - zz[-1])
        for s, t in ST:   # psi_st = sum_a <A_st xi_a, p_a>, A_st = e_s e_t^T - e_t e_s^T
            out.append(sum((xi[a][t] * zz[a * d + s] - xi[a][s] * zz[a * d + t] for a in range(n)), flint.arb(0)))
        return out

    def J(zz):
        M = flint.arb_mat(U, U)
        for a in range(n):
            for k in range(d):
                M[a, a * d + k] = 2 * zz[a * d + k]
        for e, (a, b) in enumerate(E):
            r = n + e
            for k in range(d):
                M[r, a * d + k] = zz[b * d + k]
                M[r, b * d + k] = zz[a * d + k]
            M[r, U - 1] = flint.arb(-1)
        for l, (s, t) in enumerate(ST):
            r = n + m + l
            for a in range(n):
                M[r, a * d + s] = xi[a][t]
                M[r, a * d + t] = -xi[a][s]
        return M

    Jm = J([flint.arb(x.mid()) for x in z])
    Jf = np.array([[float(Jm[i, j].mid()) for j in range(U)] for i in range(U)])
    R = np.linalg.inv(Jf)                                                         # untrusted
    Ra = flint.arb_mat([[flint.arb(float(x)) for x in row] for row in R])
    Fz = flint.arb_mat([[x] for x in F(z)])
    RF = Ra * Fz
    ybar = max(float(abs(RF[i, 0]).upper()) for i in range(U))
    I = flint.arb_mat([[flint.arb(int(i == j)) for j in range(U)] for i in range(U)])
    for fac in (4.0, 64.0, 1024.0):
        r = max(ybar * fac, 2.0 ** -(PREC - 20))
        box = [x + flint.arb(0, r) for x in z]
        Km = I - Ra * J(box)
        Kv = -RF + Km * flint.arb_mat([[flint.arb(0, r)] for _ in range(U)])
        mx = max(float(abs(Kv[i, 0]).upper()) for i in range(U))
        if mx < r:
            c = z[-1] + Kv[U - 1, 0]
            return {"arb_pass": True, "r": r, "max_over_r": mx / r, "c_ball": c.str(35, radius=True)}, c
    return {"arb_pass": False, "r": r, "max_over_r": mx / r}, None


def main():
    outp = sys.argv[1]
    args = sys.argv[2:]
    umax = int(args[args.index("--max-unknowns") + 1]) if "--max-unknowns" in args else 400
    every = int(args[args.index("--every") + 1]) if "--every" in args else 20
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    cand = []
    import glob
    for f in sorted(glob.glob(os.path.join(CERT, "census", "results_v4", "*.json"))):
        j = json.load(open(f))
        if j["status"] not in ("PASS", "LOCALMIN"):
            continue
        a = [x for x in j["attempts"] if x.get("status") == j["status"]][0]
        if a["unknowns"] <= umax:
            cand.append((j["status"], a["unknowns"], j["d"], j["n"], a["tol"], a["c_interval_30"]))
    cand.sort()
    sample = cand[::every]
    done = set()
    if os.path.exists(outp):
        done = {(r["d"], r["n"]) for r in map(json.loads, open(outp))}
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    for st, U, d, n, tol, c30 in sample:
        if (d, n) in done:
            continue
        t0 = time.time()
        rec = {"d": d, "n": n, "census_status": st, "unknowns": U, "tol": tol}
        pts = snapshot.read_points(SNAP, idx[(d, n)])
        X = np.array([[float(x) for x in p] for p in pts])
        X = X / np.linalg.norm(X, axis=1, keepdims=True)
        keep = {}
        cq.certify(X, d, n, tol, keep=keep)
        info, cball = arb_check(keep, d)
        rec.update(info)
        if cball is not None:
            lo = flint.arb(flint.fmpq(int(c30[0]), 10 ** 30))
            hi = flint.arb(flint.fmpq(int(c30[1]), 10 ** 30))
            rec["overlaps_census_interval"] = bool(cball.overlaps(lo.union(hi)))
        rec["seconds"] = round(time.time() - t0, 1)
        with open(outp, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(d, n, st, U, rec.get("arb_pass"), rec["seconds"], flush=True)
    flint.ctx.prec = old


if __name__ == "__main__":
    main()
