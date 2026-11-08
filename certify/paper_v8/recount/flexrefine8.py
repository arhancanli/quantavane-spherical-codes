"""The negative certificate (Theorem 4.1) decided exactly with a third right inverse, for the entries where the two of
flexexact8.py fail (census version 8).

  nice -n 10 python3 -I flexrefine8.py FLEX8_OUTDIR OUTDIR

For every entry whose record in FLEX8_OUTDIR (flexexact8.py) is not FLEX: the same certificate with the right inverse
C = V S^-1 U^T from a singular value decomposition L ~ U S V^T in double precision (UNTRUSTED, like the others: any
matrix will do). Everything TRUSTED is the code of flexexact8.py, imported unchanged: the reference configuration X0,
L0 = L(X0), beta and ||C'||_inf by exact sparse products (exact_beta), and then h, tau' and delta exactly as in its
check(), whose second half is repeated here verbatim (flexexact8.py is not changed, so that its records keep its SHA-256).
"""
import hashlib, json, os, sys, time, traceback
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import flexexact8 as FX
import numpy as np
import flint

SHA = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()


def check_svd(X, dec, d, n, tol, a):
    rec = {}
    I_, J_, Gf = FX.proposal(X, n, tol)
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
    Xi = np.rint(np.ldexp(X, FX.GRID)).astype(np.int64)
    Li, Lf = FX.build_L(Xi, X, pts, loc, I_, J_, d)
    L1 = Li.copy()
    L1.data = Li.data >> 26
    L0 = Li.copy()
    L0.data = Li.data & FX.M26
    if UL * UC > FX.QR_MAX:
        rec["status"] = "FAIL: too large for a dense decomposition"
        return rec
    U, S, Vt = np.linalg.svd(Lf.toarray(), full_matrices=False)
    rec["sigma_min_double"] = float(S[-1])
    Cf = (Vt.T / S) @ U.T
    del U, Vt
    beta, cnorm, t = FX.exact_beta(Li, L1, L0, lambda j0, j1: Cf[:, j0:j1], UL, UL, UC)
    del Cf
    if beta is None:
        rec["status"] = "FAIL: no right inverse"
        return rec
    rec.update(right_inverse="SVD", beta=FX.up_float(beta), C_norm_inf=FX.up_float(cnorm))
    if not beta < 1:
        rec["status"] = "FAIL: beta >= 1"
        return rec
    h0 = (1 - beta) / (2 * cnorm)
    # ---- from here verbatim as flexexact8.check (FX. prefixes added)
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
        rec["pairs_checked"] = 0
    else:
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
    return rec


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    idx = {(e["d"], e["n"]): e for e in FX.snapshot.load_index()}
    todo = []
    for f in sorted(os.listdir(src)):
        if f.endswith(".json"):
            j = json.load(open(os.path.join(src, f)))
            if j.get("status") != "FLEX":
                todo.append((j["d"], j["n"]))
    print(len(todo), "entries", flush=True)
    for d, n in todo:
        js = os.path.join(out, f"d{d}_n{n}.json")
        if os.path.exists(js):
            continue
        r = json.load(open(os.path.join(FX.CERT, "census", "results_v4", f"d{d}_n{n}.json")))
        a = [x for x in r["attempts"] if x.get("status") == "FLEX"][0]
        t0, c0 = time.time(), time.process_time()
        rec = {"d": d, "n": n, "tol": a["tol"], "program_sha256": SHA, "flexexact8_sha256": FX.SHA}
        try:
            pts = FX.snapshot.read_points(FX.SNAP, idx[(d, n)])
            dec = [[t_.strip() for t_ in pp] for pp in pts]
            X = np.array([[float(x) for x in pp] for pp in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            rec.update(check_svd(X, dec, d, n, a["tol"], a))
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
        print(d, n, rec.get("status"), rec.get("beta"), rec.get("delta"), rec["cpu_seconds"], flush=True)


if __name__ == "__main__":
    main()
