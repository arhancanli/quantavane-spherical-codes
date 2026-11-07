# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""Entries of the archived table near (22,116) whose stored order lists the points orbit by orbit (UNTRUSTED floating
point: a screening, not a certificate).

  python3 -I certify/note/family_table.py SNAPDIR OUT.json [--dmin 18] [--dmax 26] [--nmin 1] [--nmax 250] [--kmax 60]
                                                           [--aut-seconds 600]

For every entry with dmin <= d <= dmax and N <= nmax: read the stored decimals (double precision), normalise, and
cluster the inner products (a new class at every gap larger than 1e-9 in the sorted values). Entries with more than kmax
classes are counted as unstructured and not analysed further. For the others: the symmetry group of the class matrix
(certify/note/aut116.py, individualisation and refinement; with distinct class values and a spanning code these are the
orthogonal maps permuting the code; an entry whose group is not found within --aut-seconds is recorded as such), its
orbits, and whether the stored order lists the points orbit by orbit (at most
2 r maximal runs of consecutive points of one orbit, r the number of orbits; for (22,116) the runs are A, B and the four
points p, q, p, q). For the entries grouped by orbits (nontrivial group, at least two orbits): the
contacts (the largest class), rattlers, the rank of the core, the infinitesimal flexes of the core modulo infinitesimal
rotations (dense SVD of the rigidity matrix), the space of equilibrium stresses on the contacts (SVD) and a strictly
positive one (LP: maximise the smallest entry), and, if there are flexes, the second-order energy H on them for that
stress (Theorem 2.5 of the note: a strictly positive stress with H > 0 on every flex not a rotation proves jamming).
"""
import argparse, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [p for p in os.environ.get("QV_LIBS", "").split(os.pathsep) if p]   # optional library folders
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "census"))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.seterr(all="ignore")
from scipy.optimize import linprog
from snapshot import load_index, read_points
from aut116 import automorphisms, orbit


def classes(X, gap=1e-9):
    G = X @ X.T
    n = len(X)
    iu = np.triu_indices(n, 1)
    v = G[iu]
    order = np.argsort(v)
    sv = v[order]
    cuts = np.flatnonzero(np.diff(sv) > gap)
    lab_sorted = np.zeros(len(sv), dtype=np.int64)
    lab_sorted[cuts + 1] = 1
    lab_sorted = np.cumsum(lab_sorted)
    lab = np.empty_like(lab_sorted)
    lab[order] = lab_sorted
    C = -np.ones((n, n), dtype=np.int64)
    C[iu] = lab
    C.T[iu] = lab
    spread = max(float(sv[lab_sorted == k].max() - sv[lab_sorted == k].min()) for k in range(int(lab_sorted.max()) + 1))
    return C, int(lab_sorted.max()) + 1, spread, G


def rigidity(X, C, kmax_cls):
    n, d = X.shape
    G = X @ X.T
    top = int(C.max())
    E = np.array([(a, b) for a in range(n) for b in range(a + 1, n) if C[a, b] == top], dtype=np.int64)
    mu = float(G[E[:, 0], E[:, 1]].max())
    deg = np.bincount(E.ravel(), minlength=n)
    core = np.flatnonzero(deg > 0)
    rep = {"contacts": int(len(E)), "largest_cosine": mu, "rattlers": int(n - len(core)),
           "degrees": sorted(set(int(x) for x in deg))}
    idx = {int(a): i for i, a in enumerate(core)}
    Xc = X[core]
    m = len(core)
    rep["rank_core"] = int(np.linalg.matrix_rank(Xc, tol=1e-8))
    Ec = np.array([(idx[int(a)], idx[int(b)]) for a, b in E])
    # rigidity matrix: one row per contact and per point (norm)
    R = np.zeros((len(Ec) + m, m * d))
    for r, (a, b) in enumerate(Ec):
        R[r, a * d:(a + 1) * d] = Xc[b]
        R[r, b * d:(b + 1) * d] = Xc[a]
    for a in range(m):
        R[len(Ec) + a, a * d:(a + 1) * d] = Xc[a]
    _, sv, Vt = np.linalg.svd(R, full_matrices=R.shape[0] < R.shape[1])
    tol = 1e-7 * sv.max()
    rk = int((sv > tol).sum())
    nullity = m * d - rk
    rots = []
    for i in range(d):
        for j in range(i + 1, d):
            A = np.zeros((d, d)); A[i, j], A[j, i] = 1, -1
            rots.append((Xc @ A.T).ravel())
    Kr = np.array(rots).T
    dimK = int(np.linalg.matrix_rank(Kr, tol=1e-8))
    rep["flexes_mod_rotations"] = int(nullity - dimK)
    svs = np.sort(sv)
    rep["singular_values_around_gap"] = [float(svs[len(svs) - rk - 1]) if len(svs) > rk else None,
                                         float(svs[len(svs) - rk]) if rk > 0 else None]
    # equilibrium stresses: sum_b w_ab (x_b - mu x_a) = 0 at every point a of the core
    Eq = np.zeros((m * d, len(Ec)))
    for k, (a, b) in enumerate(Ec):
        Eq[a * d:(a + 1) * d, k] += Xc[b] - mu * Xc[a]
        Eq[b * d:(b + 1) * d, k] += Xc[a] - mu * Xc[b]
    # the numerical row space of Eq (orthonormal rows); the stresses are its orthogonal complement
    _, s2, V2 = np.linalg.svd(Eq, full_matrices=False)
    r2 = int((s2 > 1e-7 * s2.max()).sum())
    Vr = V2[:r2]
    rep["stress_dimension"] = int(len(Ec) - r2)
    omega = None
    if len(Ec) > r2:
        ne = len(Ec)
        c = np.zeros(ne + 1); c[-1] = -1
        lp = linprog(c, A_ub=np.hstack([-np.eye(ne), np.ones((ne, 1))]), b_ub=np.zeros(ne),
                     A_eq=np.vstack([np.hstack([Vr, np.zeros((r2, 1))]), np.concatenate([np.ones(ne), [0]])[None, :]]),
                     b_eq=np.concatenate([np.zeros(r2), [1.0]]), bounds=[(None, None)] * (ne + 1), method="highs")
        if lp.status == 0:
            omega = lp.x[:ne]
            omega = omega - Vr.T @ (Vr @ omega)                  # project onto the stresses
            rep["maxmin_stress_times_E"] = float(omega.min() * ne / omega.sum())
            rep["equilibrium_residual"] = float(np.abs(Eq @ omega).max())
    rep["strictly_positive_stress"] = bool(omega is not None and rep["maxmin_stress_times_E"] > 1e-6)
    if rep["flexes_mod_rotations"] > 0 and omega is not None:
        Nul = Vt[rk:].T
        Q, _ = np.linalg.qr(Kr)
        Fl = Nul - Q @ (Q.T @ Nul)
        U, s3, _ = np.linalg.svd(Fl, full_matrices=False)
        Fl = U[:, :int((s3 > 1e-6).sum())]
        H = np.zeros((Fl.shape[1], Fl.shape[1]))
        for i in range(Fl.shape[1]):
            for j in range(i, Fl.shape[1]):
                yi, yj = Fl[:, i].reshape(m, d), Fl[:, j].reshape(m, d)
                a, b = Ec[:, 0], Ec[:, 1]
                H[i, j] = H[j, i] = float(np.sum(omega * (np.sum(yi[a] * yj[b], 1) + np.sum(yj[a] * yi[b], 1)
                                                          - mu * np.sum(yi[a] * yj[a], 1) - mu * np.sum(yi[b] * yj[b], 1))))
        ev = np.linalg.eigvalsh(H)
        rep["H_eigs_min_max"] = [float(ev.min()), float(ev.max())]
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("snapdir")
    ap.add_argument("out")
    ap.add_argument("--dmin", type=int, default=18)
    ap.add_argument("--dmax", type=int, default=26)
    ap.add_argument("--nmax", type=int, default=250)
    ap.add_argument("--kmax", type=int, default=60)
    ap.add_argument("--notes", default="", help="only entries whose references and notes contain this text")
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--aut-seconds", type=int, default=600, help="time limit for the symmetry group of one entry")
    args = ap.parse_args()
    t0 = time.time()
    entries = [e for e in load_index() if args.dmin <= e["d"] <= args.dmax and args.nmin <= e["n"] <= args.nmax
               and args.notes in e["notes"]]
    rng = dict(vars(args), snapdir="SNAPDIR (the archived table, DSpace@MIT 1721.1/153543)")  # 7 Nov: no local paths in outputs
    out = {"program": "family_table.py", "range": rng, "entries": len(entries), "structured": [], "grouped": []}
    for e in entries:
        d, n = e["d"], e["n"]
        X = np.array([[float(t) for t in row] for row in read_points(args.snapdir, e)])
        X /= np.linalg.norm(X, axis=1)[:, None]
        C, K, spread, G = classes(X)
        if K > args.kmax:
            continue
        import signal
        def _alarm(signum, frame):
            raise TimeoutError
        signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(args.aut_seconds)
        try:
            order, base, lengths, gens, cells = automorphisms(C)
        except TimeoutError:
            out["structured"].append({"d": d, "n": n, "notes": e["notes"], "classes": K, "class_spread": spread,
                                      "aut_timeout_seconds": args.aut_seconds})
            json.dump(out, open(args.out, "w"), indent=1)
            print(json.dumps(out["structured"][-1]), flush=True)
            continue
        finally:
            signal.alarm(0)
        lab = -np.ones(n, dtype=np.int64)
        orbs = []
        for a in range(n):
            if lab[a] < 0:
                o = sorted(orbit(a, gens)) if gens else [a]
                for x in o:
                    lab[x] = len(orbs)
                orbs.append(o)
        runs = 1 + int(np.count_nonzero(np.diff(lab)))       # maximal blocks of consecutive points in one orbit
        contiguous = runs <= 2 * len(orbs)
        rec = {"d": d, "n": n, "notes": e["notes"], "classes": K, "class_spread": spread, "aut_order": int(order),
               "orbit_sizes": [len(o) for o in orbs], "runs_in_stored_order": runs,
               "stored_order_grouped_by_orbits": bool(contiguous)}
        out["structured"].append(rec)
        if contiguous and order > 1 and len(orbs) >= 2:
            try:
                rec.update(rigidity(X, C, args.kmax))
            except Exception as ex:
                rec["analysis_error"] = repr(ex)
            out["grouped"].append(rec)
            print(json.dumps(rec), flush=True)
        json.dump(out, open(args.out, "w"), indent=1)
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, open(args.out, "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("structured", "grouped")}))


if __name__ == "__main__":
    main()
