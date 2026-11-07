"""UNTRUSTED survey (numerical only): symmetry of the open hyperstatic entries and the size of the symmetry-reduced
equal-contact system.

  nice -n 10 python3 -I sym_survey.py OUT.jsonl [--shard I/K] [--only d_n,...]

For each entry of class 'local' (LOCALMIN, strictness open): the census certificate is re-run (cq.certify, keep) to get
V, E, E' and the zero P' of the square subsystem in double precision. The automorphisms of P' (permutations of V that
preserve every inner product to 1e-9, realised by an orthogonal map) are searched by backtracking on a base of d points;
up to MAXGEN random automorphisms are kept as generators. Reported: the order of the group they generate (counted by
closure up to a cap), point orbits t, edge orbits k_E of E, dim Fix (configurations with p_{pi(i)} = h p_i), the
dimension c of the centraliser of the representation in so(d), whether the reduced system (t norms, k_E contacts, c
slice equations; dim Fix + 1 unknowns) is square, and its smallest singular value at P'.
Nothing here is a proof; the survey decides only where an exact symmetric certificate is worth attempting.
"""
import json, os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/flintlib"))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, CERT)
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import cq
import snapshot

TOLG = 1e-9
MAXGEN = 24
CAP = 200000


def classes(X, G, tol):
    """Untrusted invariants of each point: rounded power sums of its inner products."""
    inv = np.stack([np.sort(G, 1)[:, -2]] + [(G ** k).sum(1) for k in (2, 3, 4)], 1)
    lab = np.zeros(len(X), int)
    order = np.lexsort(inv.T[::-1])
    cur = 0
    lab[order[0]] = 0
    for a, b in zip(order[:-1], order[1:]):
        if np.max(np.abs(inv[a] - inv[b]) / (1 + np.abs(inv[a]))) > 1e3 * tol:
            cur += 1
        lab[b] = cur
    return lab


def find_auts(X, tol=TOLG, maxfound=MAXGEN, seed=0, tlimit=60.0):
    n, d = X.shape
    G = X @ X.T
    lab = classes(X, G, tol)
    import scipy.linalg as sl
    _, _, piv = sl.qr(X.T, pivoting=True, mode="economic")
    base = [int(i) for i in piv[:d]]
    XB = X[base]
    try:
        XBinv = np.linalg.inv(XB)
    except np.linalg.LinAlgError:
        return [], base
    rng = random.Random(seed)
    found = []
    seen = set()
    t0 = time.time()
    cand0 = [j for j in range(n) if lab[j] == lab[base[0]]]

    def complete(img):
        H = (XBinv @ X[img]).T          # H x_b = x_c  (rows: X[img] = XB H^T)
        if np.max(np.abs(H.T @ H - np.eye(d))) > 1e-6:
            return None
        Y = X @ H.T
        S = Y @ X.T
        perm = np.argmax(S, 1)
        if np.min(S[np.arange(n), perm]) < 1 - 1e-7 or len(set(perm.tolist())) != n:
            return None
        if np.max(np.abs(G[np.ix_(perm, perm)] - G)) > 1e3 * tol:
            return None
        return perm, H

    def dfs(img):
        if len(found) >= maxfound or time.time() - t0 > tlimit:
            return
        k = len(img)
        if k == d:
            r = complete(img)
            if r is not None:
                key = tuple(r[0].tolist())
                if key not in seen:
                    seen.add(key)
                    found.append(r)
            return
        b = base[k]
        cands = [j for j in range(n) if lab[j] == lab[b] and j not in img
                 and all(abs(G[img[i], j] - G[base[i], b]) < 1e3 * tol for i in range(k))]
        rng.shuffle(cands)
        for j in cands:
            dfs(img + [j])
            if len(found) >= maxfound or time.time() - t0 > tlimit:
                return

    # several randomised restarts, each takes the first automorphism it reaches
    for rep in range(maxfound * 3):
        if len(found) >= maxfound or time.time() - t0 > tlimit:
            break
        before = len(found)
        lim = len(found) + 1
        rng.shuffle(cand0)
        for c0 in cand0:
            img = [c0]
            # greedy random descent with backtracking limited to this start
            stack = [img]
            ok = False

            def dfs1(img):
                nonlocal ok
                if ok or time.time() - t0 > tlimit:
                    return
                k = len(img)
                if k == d:
                    r = complete(img)
                    if r is not None:
                        key = tuple(r[0].tolist())
                        if key not in seen:
                            seen.add(key)
                            found.append(r)
                            ok = True
                    return
                b = base[k]
                cands = [j for j in range(n) if lab[j] == lab[b] and j not in img
                         and all(abs(G[img[i], j] - G[base[i], b]) < 1e3 * tol for i in range(k))]
                rng.shuffle(cands)
                for j in cands:
                    dfs1(img + [j])
                    if ok:
                        return
            dfs1(img)
            if ok:
                break
        if len(found) == before:
            break
    return found, base


def closure_order(gens, n, cap=CAP):
    """Order of the permutation group generated (orbit of the identity under right multiplication), up to cap."""
    ident = tuple(range(n))
    seen = {ident}
    frontier = [ident]
    gl = [tuple(g.tolist()) for g in gens]
    while frontier and len(seen) < cap:
        nxt = []
        for p in frontier:
            for g in gl:
                q = tuple(p[i] for i in g)
                if q not in seen:
                    seen.add(q)
                    nxt.append(q)
        frontier = nxt
    return len(seen) if len(seen) < cap else None


def union_orbits(n, gens_perm):
    par = list(range(n))

    def f(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a
    for g in gens_perm:
        for i in range(n):
            a, b = f(i), f(int(g[i]))
            if a != b:
                par[a] = b
    return [f(i) for i in range(n)]


def survey_entry(X, d, n, tol):
    keep = {}
    log = cq.certify(X, d, n, tol, keep=keep)
    if log.get("status") != "LOCALMIN":
        return {"status": "not LOCALMIN on re-run: " + str(log.get("status"))}
    V, E, Esq, extra = keep["V"], keep["E"], keep["Esq"], keep["extra"]
    Xf = keep["Xf"]
    P = Xf / np.linalg.norm(Xf, axis=1, keepdims=True)
    nV = len(V)
    out = {"nV": nV, "contacts": len(E), "iso": nV * (d - 1) - d * (d - 1) // 2 + 1, "extra": len(extra)}
    t0 = time.time()
    found, base = find_auts(P)
    out["aut_seconds"] = round(time.time() - t0, 2)
    gens = [f[0] for f in found]
    Hs = [f[1] for f in found]
    out["n_auts_found"] = len(gens)
    if not gens:
        out["status"] = "no automorphism found"
        return out
    order = closure_order(gens, nV) if nV <= 400 else None
    out["group_order"] = order
    Eset = set(tuple(sorted(e)) for e in E)
    inv = all(tuple(sorted((int(g[a]), int(g[b])))) in Eset for g in gens for (a, b) in E)
    out["E_invariant"] = inv
    po = union_orbits(nV, gens)
    t = len(set(po))
    # edge orbits
    eidx = {e: k for k, e in enumerate(sorted(Eset))}
    elist = sorted(Eset)
    par = list(range(len(elist)))

    def f(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a
    for g in gens:
        for k, (a, b) in enumerate(elist):
            im = tuple(sorted((int(g[a]), int(g[b]))))
            if im in eidx:
                x, y = f(k), f(eidx[im])
                if x != y:
                    par[x] = y
    kE = len(set(f(k) for k in range(len(elist))))
    # Fix: p_{g(i)} = H p_i  (rows: P[g] = P H^T); kernel of the stacked map, dimension nV*d
    import scipy.sparse as sps
    N = nV * d
    blocks = []
    for g, H in zip(gens, Hs):
        # rows (i, k): P[g(i), k] - sum_l H[k, l] P[i, l]
        r, c, v = [], [], []
        for i in range(nV):
            gi = int(g[i])
            for k in range(d):
                row = i * d + k
                r.append(row); c.append(gi * d + k); v.append(1.0)
                for l in range(d):
                    r.append(row); c.append(i * d + l); v.append(-H[k, l])
        blocks.append(sps.csr_matrix((v, (r, c)), shape=(N, N)))
    A = sps.vstack(blocks).tocsr()
    AtA = (A.T @ A).toarray()
    ev, evec = np.linalg.eigh(AtA)
    fdim = int(np.sum(ev < 1e-8 * max(1.0, ev.max())))
    Bf = evec[:, :fdim]                       # orthonormal basis of Fix (numerical)
    out["dim_fix"] = fdim
    # centraliser in so(d): A H = H A for all generators, A skew
    ST = [(s, tt) for s in range(d) for tt in range(s + 1, d)]
    rows = []
    for H in Hs:
        for s, tt in ST:
            Aq = np.zeros((d, d)); Aq[s, tt] = 1; Aq[tt, s] = -1
            rows.append((Aq @ H - H @ Aq).ravel())
    M = np.array(rows).reshape(len(Hs), len(ST), d * d).transpose(1, 0, 2).reshape(len(ST), -1).T
    sv = np.linalg.svd(M, compute_uv=False) if M.size else np.array([])
    cdim = int(np.sum(sv < 1e-8 * max(1.0, sv.max() if sv.size else 1.0))) + max(0, len(ST) - len(sv))
    out["dim_centraliser"] = cdim
    out["point_orbits"] = t
    out["edge_orbits"] = kE
    neq = t + kE + cdim
    nun = fdim + 1
    out["reduced_equations"], out["reduced_unknowns"] = neq, nun
    out["square"] = neq == nun
    # reduced Jacobian at P' (orbit representatives; slice omitted when cdim = 0)
    reps = {}
    for i in range(nV):
        reps.setdefault(po[i], i)
    erep = {}
    for k, e in enumerate(elist):
        erep.setdefault(f(k), e)
    p = P.ravel()
    J = []
    for r in reps.values():
        row = np.zeros(N + 1); row[r * d:(r + 1) * d] = 2 * P[r]; J.append(row)
    for (a, b) in erep.values():
        row = np.zeros(N + 1); row[a * d:(a + 1) * d] = P[b]; row[b * d:(b + 1) * d] += P[a]; row[N] = -1; J.append(row)
    J = np.array(J)
    Jw = np.hstack([J[:, :N] @ Bf, J[:, N:]])
    s = np.linalg.svd(Jw, compute_uv=False)
    out["reduced_sv_min"] = float(s.min()) if s.size else None
    out["reduced_sv_max"] = float(s.max()) if s.size else None
    out["reduced_rank_1e-8"] = int(np.sum(s > 1e-8 * s.max())) if s.size else 0
    # residual of P' in Fix (how symmetric P' is)
    out["Pprime_fix_residual"] = float(np.abs(A @ p).max())
    out["status"] = "ok"
    return out


def main():
    outp = sys.argv[1]
    args = sys.argv[2:]
    sh, nsh = 0, 1
    only = None
    if "--shard" in args:
        sh, nsh = (int(x) for x in args[args.index("--shard") + 1].split("/"))
    if "--only" in args:
        only = set(args[args.index("--only") + 1].split(","))
    s5 = json.load(open(os.path.join(CERT, "v5", "out", "summary_v5.json")))
    opens = sorted(tuple(map(int, k.split(","))) for k, v in s5["class"].items() if v[0] == "local")
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    done = set()
    if os.path.exists(outp):
        for line in open(outp):
            j = json.loads(line)
            done.add((j["d"], j["n"]))
    todo = []
    for d, n in opens:
        if only and f"{d}_{n}" not in only:
            continue
        r = json.load(open(os.path.join(CERT, "census", "results_v4", f"d{d}_n{n}.json")))
        a = [x for x in r["attempts"] if x.get("status") == "LOCALMIN"][0]
        todo.append((a["unknowns"], d, n, a["tol"]))
    todo.sort()
    for k, (U, d, n, tol) in enumerate(todo):
        if k % nsh != sh or (d, n) in done:
            continue
        t0 = time.time()
        rec = {"d": d, "n": n, "U": U, "tol": tol}
        try:
            pts = snapshot.read_points(SNAP, idx[(d, n)])
            X = np.array([[float(x) for x in p] for p in pts])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            rec.update(survey_entry(X, d, n, tol))
        except Exception as e:
            import traceback
            rec["status"] = "error: " + traceback.format_exc()[-400:]
        rec["seconds"] = round(time.time() - t0, 1)
        with open(outp, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(d, n, U, rec.get("status"), rec.get("group_order"), rec.get("square"), rec["seconds"], flush=True)


if __name__ == "__main__":
    main()
