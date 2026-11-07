"""Exact first- and second-order rigidity of the contact core of a table entry (paper, version 5).

  python3 -I certify/v5/rigidity.py SNAPDIR D N OUT.json [--disc DISC]

core.py builds the exact core configuration Q* near the stored code (exact Gram matrix of the core over K = Q(sqrt DISC),
rank r, contacts E at the exact cosine c, rattlers below c in ball arithmetic, certified distance to the stored code).
This program requires r = d (the core spans R^d) and decides, exactly over K:

 (R) the rank of the bar rigidity matrix R of the core (rows <p_a, w_a> and <p_a, w_b> + <p_b, w_a>, coordinates
     w_a = sum_k u_ak p_{b_k}) over F_q (q prime, sqrt DISC -> a square root mod q), a lower bound for its rank over K;
     the infinitesimal rotations and, when the contact graph is disconnected, the rotations of single components (the
     others fixed) are exact kernel vectors (checked), independent mod q. If the rank bound equals |V|d minus their
     number, they span the kernel exactly: no other infinitesimal flex exists.
 (S) a strictly positive equilibrium stress omega in K^E, checked exactly (sum_b omega_ab W_b = lambda_a W_a for every
     core point a, every omega_e > 0). Candidate: the uniform stress, else the exact solution of a square system for
     a stress whose entries off an isostatic subset E' are the (rounded) entries of a linear-programming stress
     (UNTRUSTED choices; only the exact check counts).
 (H) on the flexes found in (R): the matrix of H_omega (paper, stress-energy lemma); whether it vanishes identically,
     and its exact LDL^T pivots.
If no stress is positive on every contact, the contacts E+ that carry a positive entry in some stress are found by a linear
program (UNTRUSTED), and (R) and (S) are decided exactly for (core, E+); by the remark after Theorem 2.5, a bar-rigid E+
with a strictly positive stress still makes Q* a strict local minimum (the other contacts are at c but carry no stress).
Verdicts: STRICT (bar-rigid with a strictly positive stress: Theorem 2.5, the core of Q* is infinitesimally jammed and Q*
is a strict local minimum modulo rotations and rattlers); LOCAL, NOT STRICT (the contact graph has several components,
each bar-rigid with a strictly positive stress: Q* is a local minimum, its nearby minimisers are exactly the independent
rotations of the components, and it is not jammed).
"""
import argparse, hashlib, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import core
import exact2 as X2
from kfield import K, ZERO, ONE, ksum, fq, is_pos_def, prime_for
import flint
import numpy as np

SHA = {}
for _f in ("rigidity.py", "core.py", "exact2.py", "kfield.py", "dist.py"):
    SHA[_f] = hashlib.sha256(open(os.path.join(HERE, _f), "rb").read()).hexdigest()


def components(nV, E):
    par = list(range(nV))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b in E:
        ra, rb = find(a), find(b)
        if ra != rb:
            par[ra] = rb
    comp = {}
    for v in range(nV):
        comp.setdefault(find(v), set()).add(v)
    return sorted(comp.values(), key=min)


def rotations(nV, W, GBBinv, d, comp=None):
    out = []
    for s in range(d):
        for t in range(s + 1, d):
            u = []
            for a in range(nV):
                if comp is not None and a not in comp:
                    u.extend([ZERO] * d)
                    continue
                Sw = [ZERO] * d
                Sw[s] = W[t][a]
                Sw[t] = -W[s][a]
                u.extend(ksum(GBBinv[k][l] * Sw[l] for l in range(d) if not Sw[l].iszero()) for k in range(d))
            out.append(u)
    return out


def apply_R(G, B, E, u, d, nV):
    out = []
    for a in range(nV):
        out.append(ksum(G[a][B[k]] * u[a * d + k] for k in range(d) if not u[a * d + k].iszero()))
    for a, b in E:
        out.append(ksum(G[a][B[k]] * u[b * d + k] for k in range(d) if not u[b * d + k].iszero()) +
                   ksum(G[b][B[k]] * u[a * d + k] for k in range(d) if not u[a * d + k].iszero()))
    return out


def R_rows_modp(G, B, E, d, nV, p, r):
    rows = []
    for a in range(nV):
        rows.append({a * d + k: G[a][B[k]].modp(p, r) for k in range(d)})
    for a, b in E:
        row = {}
        for k in range(d):
            row[b * d + k] = (row.get(b * d + k, 0) + G[a][B[k]].modp(p, r)) % p
            row[a * d + k] = (row.get(a * d + k, 0) + G[b][B[k]].modp(p, r)) % p
        rows.append(row)
    return rows


def check_stress(W, E, d, nV, omega):
    nb = [[] for _ in range(nV)]
    for e, (a, b) in enumerate(E):
        nb[a].append((b, omega[e]))
        nb[b].append((a, omega[e]))
    for a in range(nV):
        s = [ksum(w * W[k][b] for b, w in nb[a]) for k in range(d)]
        wa = [W[k][a] for k in range(d)]
        k0 = next(k for k in range(d) if not wa[k].iszero())
        lam = s[k0] / wa[k0]
        if any(not (s[k] == lam * wa[k]) for k in range(d)):
            return False
    return True


def lp_stress(XV, E, d):
    """UNTRUSTED: a floating-point stress of largest smallest entry (total 1)."""
    from scipy.optimize import linprog
    import scipy.sparse as sps
    nV, m = len(XV), len(E)
    rows, cols, vals = [], [], []
    for e, (a, b) in enumerate(E):
        for k in range(d):
            rows += [a * d + k, b * d + k]
            cols += [e, e]
            vals += [XV[b, k], XV[a, k]]
    for a in range(nV):
        for k in range(d):
            rows.append(a * d + k)
            cols.append(m + a)
            vals.append(XV[a, k])
    for e in range(m):
        rows.append(nV * d)
        cols.append(e)
        vals.append(1.0)
    Aeq = sps.csr_matrix((vals, (rows, cols)), shape=(nV * d + 1, m + nV + 1))
    beq = np.zeros(nV * d + 1)
    beq[-1] = 1.0
    Aub = sps.hstack([-sps.identity(m, format="csr"), sps.csr_matrix((m, nV)), sps.csr_matrix(np.ones((m, 1)))],
                     format="csr")
    cobj = np.zeros(m + nV + 1)
    cobj[-1] = -1
    res = linprog(cobj, A_ub=Aub, b_ub=np.zeros(m), A_eq=Aeq, b_eq=beq, bounds=[(None, None)] * (m + nV + 1),
                  method="highs-ipm", options={"time_limit": 1800.0})
    return (res.x[:m], -res.fun) if res.status == 0 else (None, None)


def lp_support(XV, E, d):
    """UNTRUSTED: the contacts that carry a positive entry in some nonnegative equilibrium stress (the support of a
    stress in the relative interior of the cone): maximise sum z_e with 0 <= z_e <= 1, z_e <= omega_e, omega >= 0 an
    equilibrium stress (unnormalised, so that every contact in the support can reach z_e = 1)."""
    from scipy.optimize import linprog
    import scipy.sparse as sps
    nV, m = len(XV), len(E)
    rows, cols, vals = [], [], []
    for e, (a, b) in enumerate(E):
        for k in range(d):
            rows += [a * d + k, b * d + k]
            cols += [e, e]
            vals += [XV[b, k], XV[a, k]]
    for a in range(nV):
        for k in range(d):
            rows.append(a * d + k)
            cols.append(m + a)
            vals.append(XV[a, k])
    # variables: omega (m), mu (nV), z (m)
    Aeq = sps.csr_matrix((vals, (rows, cols)), shape=(nV * d, 2 * m + nV))
    Aub = sps.hstack([-sps.identity(m, format="csr"), sps.csr_matrix((m, nV)), sps.identity(m, format="csr")],
                     format="csr")
    cobj = np.r_[np.zeros(m + nV), -np.ones(m)]
    bounds = [(0, 1e6)] * m + [(None, None)] * nV + [(0, 1)] * m
    res = linprog(cobj, A_ub=Aub, b_ub=np.zeros(m), A_eq=Aeq, b_eq=np.zeros(nV * d), bounds=bounds,
                  method="highs-ipm", options={"time_limit": 1800.0})
    if res.status != 0:
        return None
    z = res.x[m + nV:]
    return [e for e in range(m) if z[e] > 0.5]


def exact_stress_from_lp(G, B, W, E, d, nV, XV, wf, disc, log):
    """UNTRUSTED construction of an exact candidate stress: choose an isostatic set of unknowns E' (pivoted QR of the
    numerical equilibrium matrix), fix omega off E' to the rounded LP values (rationals), and solve the exact linear
    system sum_b omega_ab (W_b - c W_a) = 0 on a set of independent rows for omega on E' (over K = Q(sqrt disc): in real
    form, omega_e = alpha_e + beta_e sqrt(disc), a rational system of twice the size). The exact check is check_stress."""
    import scipy.linalg as sl
    m = len(E)
    c = float(G[E[0][0]][E[0][1]])
    A = np.zeros((nV * d, m))
    for e, (a, b) in enumerate(E):
        A[a * d:(a + 1) * d, e] += XV[b] - c * XV[a]
        A[b * d:(b + 1) * d, e] += XV[a] - c * XV[b]
    _, Rr, piv = sl.qr(A, pivoting=True, mode="economic")
    dg = np.abs(np.diag(Rr))
    rank = int((dg > 1e-9 * dg[0]).sum())
    Ep = sorted(int(i) for i in piv[:rank])
    cK = G[E[0][0]][E[0][1]]
    rowsQ = {}
    for e, (a, b) in enumerate(E):
        for k in range(d):
            for (x, y) in ((a, b), (b, a)):
                val = W[k][y] - cK * W[k][x]
                if not val.iszero():
                    rowsQ.setdefault(x * d + k, {})
                    rowsQ[x * d + k][e] = rowsQ[x * d + k].get(e, ZERO) + val
    _, _, prow = sl.qr(A[:, Ep].T, pivoting=True, mode="economic")
    sel = sorted(int(i) for i in prow[:rank])
    wq = [fq(int(round(x * 2 ** 50)), 2 ** 50) for x in wf]
    pos = {e: i for i, e in enumerate(Ep)}
    two = 1 if disc == 1 else 2
    M = flint.fmpq_mat(two * rank, two * rank)
    rhs = flint.fmpq_mat(two * rank, 1)
    for i, rr in enumerate(sel):
        acc_a, acc_b = fq(0), fq(0)
        for e, v in rowsQ.get(rr, {}).items():
            if e in pos:
                j = pos[e]
                if two == 1:
                    M[i, j] = v.a
                else:
                    # (a + b s)(alpha + beta s) = (a alpha + D b beta) + (b alpha + a beta) s
                    M[2 * i, 2 * j] = v.a
                    M[2 * i, 2 * j + 1] = v.b * disc
                    M[2 * i + 1, 2 * j] = v.b
                    M[2 * i + 1, 2 * j + 1] = v.a
            else:
                acc_a -= v.a * wq[e]
                acc_b -= v.b * wq[e]
        if two == 1:
            rhs[i, 0] = acc_a
        else:
            rhs[2 * i, 0] = acc_a
            rhs[2 * i + 1, 0] = acc_b
    try:
        sol = M.solve(rhs)
    except Exception as ex:
        log["exact_stress"] = f"square solve failed: {ex!r}"
        return None
    om = [K(wq[e]) for e in range(m)]
    for e in Ep:
        j = pos[e]
        om[e] = K(sol[j, 0]) if two == 1 else K(sol[2 * j, 0], sol[2 * j + 1, 0])
    log["exact_stress_square_size"] = two * rank
    return om


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("snapdir")
    ap.add_argument("d", type=int)
    ap.add_argument("n", type=int)
    ap.add_argument("out")
    ap.add_argument("--disc", type=int, default=1)
    args = ap.parse_args()
    t0 = time.time()
    d, n = args.d, args.n
    log = {"d": d, "n": n, "program": "rigidity.py", "code_sha256": SHA, "disc": args.disc}
    C = core.build(args.snapdir, d, n, args.disc)
    log["seconds_core"] = round(time.time() - t0, 1)
    if C["status"] != "ok":
        log.update(status=C["status"], **C.get("log", {}))
        json.dump(log, open(args.out, "w"), indent=1, default=str)
        print(json.dumps(log, default=str))
        return
    log.update(C["log"])
    G, V, E, c, r, B, W, GBB, GBBinv = (C[k] for k in ("G", "V", "E", "c", "r", "B", "W", "GBB", "GBBinv"))
    XV = C["Xf"][V]
    nV = len(V)
    if r != d:
        log["status"] = "the core does not span R^d (see subspace.py)"
        json.dump(log, open(args.out, "w"), indent=1, default=str)
        print(json.dumps({k: v for k, v in log.items() if k != "code_sha256"}, default=str))
        return
    g = d * (d - 1) // 2
    cols = nV * d
    p, rt = prime_for(args.disc)
    rk = X2.rank_mod(R_rows_modp(G, B, E, d, nV, p, rt), cols, p, sketch_seed=2026)
    log["seconds_rank"] = round(time.time() - t0, 1)
    comps = components(nV, E)
    rot = rotations(nV, W, GBBinv, d)
    flex = []
    for comp in comps[1:]:
        flex.extend(rotations(nV, W, GBBinv, d, comp=comp))
    for u in rot + flex:
        if any(not x.iszero() for x in apply_R(G, B, E, u, d, nV)):
            raise ValueError("a claimed kernel vector is not in the kernel")
    ind = X2.rank_mod(X2.kvec_rows_modp(rot + flex, p, rt), cols, p)
    nflex = len(flex)
    span_exact = (ind == g + nflex) and (cols - rk == g + nflex)
    log.update(prime=p, rank_mod_p=int(rk), cols=cols, kernel_dim_upper=cols - int(rk), components=[len(s) for s in comps],
               kernel_vectors=g + nflex, kernel_vectors_independent=bool(ind == g + nflex), flex_space_exact=bool(span_exact),
               nontrivial_flexes=(cols - int(rk) - g) if span_exact else None)
    # stress
    om = None
    unif = [ONE] * len(E)
    if check_stress(W, E, d, nV, unif):
        om, log["stress"] = unif, "uniform"
    else:
        wf, tmax = lp_stress(XV, E, d)
        log["lp_maxmin_stress"] = tmax
        log["seconds_lp"] = round(time.time() - t0, 1)
        if wf is not None and tmax > 1e-9:
            cand = exact_stress_from_lp(G, B, W, E, d, nV, XV, wf, args.disc, log)
            if cand is not None and check_stress(W, E, d, nV, cand) and all(x.sign() > 0 for x in cand):
                om = cand
                tot = ksum(om)
                mn = min(om, key=lambda x: float(x))
                log["stress"] = "exact (square solve near the linear program's)"
                log["stress_min_relative"] = float(mn) * len(E) / float(tot)
            else:
                log["stress"] = "no exact strictly positive stress found"
        elif wf is not None:
            # no stress is positive on every contact: the stressed contacts E+ (untrusted choice), then, exactly, the
            # rank of the rigidity matrix of (core, E+) and a strictly positive stress of (core, E+)
            sup = lp_support(XV, E, d)
            log["stress"] = "no strictly positive stress on all contacts"
            if sup:
                Ep = [E[e] for e in sup]
                rkp = X2.rank_mod(R_rows_modp(G, B, Ep, d, nV, p, rt), cols, p, sketch_seed=2027)
                log.update(stressed_contacts=len(Ep), rank_mod_p_stressed=int(rkp))
                wf2, t2 = lp_stress(XV, Ep, d)
                log["lp_maxmin_stress_stressed"] = t2
                if rkp == cols - g and wf2 is not None and t2 > 1e-9:
                    cand = exact_stress_from_lp(G, B, W, Ep, d, nV, XV, wf2, args.disc, log)
                    if cand is not None and check_stress(W, Ep, d, nV, cand) and all(x.sign() > 0 for x in cand):
                        tot = ksum(cand)
                        mn = min(cand, key=lambda x: float(x))
                        log["stress"] = "exact, strictly positive on the stressed contacts E+ (bar-rigid)"
                        log["stress_min_relative"] = float(mn) * len(Ep) / float(tot)
                        log["status"] = ("STRICT (a bar-rigid set of contacts with a strictly positive stress; the other "
                                         "contacts carry no stress)")
    positive = om is not None
    if positive and nflex:
        H = X2.H_matrix(GBB, list(range(nV)), E, c, om, flex, d)
        zero = all(x.iszero() for row in H for x in row)
        pd, piv = is_pos_def(H)
        log.update(H_identically_zero=bool(zero), H_pos_def=bool(pd))
    if positive and span_exact and nflex == 0:
        log["status"] = "STRICT (bar-rigid core with a strictly positive stress)"
    elif positive and span_exact and len(comps) > 1:
        log["status"] = "LOCAL, NOT STRICT (components each bar-rigid with a strictly positive stress; not jammed)"
    elif "status" not in log or not str(log.get("status", "")).startswith("STRICT"):
        log["status"] = "undecided"
    log["seconds"] = round(time.time() - t0, 1)
    json.dump(log, open(args.out, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in log.items() if k != "code_sha256"}, default=str))


if __name__ == "__main__":
    main()
