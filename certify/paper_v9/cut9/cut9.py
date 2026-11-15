"""Census version 9, item A3: the negative certificate (Theorem 4.1) at the cut of the sorted cosines that gives the
largest radius, instead of at the tolerance grid of the census.

  nice -n 10 python3 -I cut9.py OUTDIR [--only d_n,...] [--group unflagged|unresolved|sample|rest|all] [--cap 3100]

Groups: the unflagged entries, the unresolved ones, a random sample of 400 of the others (seed 20261108), the rest.

For an entry with the negative certificate in version 8 (paper_v8/recount/out_flex8, out_flex8r):
UNTRUSTED search (double precision): E_k = the k pairs of largest cosine of the double-precision code. A cut is the
window w = M - (c_k + c_{k+1}) / 2 (M the largest cosine), in the middle of the gap below the k-th pair, so that the
threshold proposal at tolerance w is exactly E_k and tau' is the whole gap c_1 - c_{k+1}. The search starts at the
version 8 window and moves by factors of 10 (up while tau' limits the radius, down while h does), then by factors
10^(1/2) and 10^(1/4) around the best window; at each window the cut is moved to the largest gap within a factor 10^(1/8)
of it. At each cut the right inverse from the normal equations gives double-precision estimates of beta, h and
delta = min(h / sqrt d, tau' / 4); a cut is admissible when L has at most --cap rows, no more rows than columns, and
beta < 1/2.
TRUSTED: the best cut found is decided by flexexact8.check (paper_v8/recount, released, unchanged) with the tolerance
w: exact beta, h, tau' and delta. If it fails, the next best cut is tried.
"""
import hashlib, json, math, os, random, sys, time, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
SNAP = os.environ.get("QV_SNAP", os.path.expanduser("~/cohn-table-dspace-153543/files"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
sys.path.insert(0, os.path.join(CERT, "paper_v8", "recount"))
sys.path.insert(0, os.path.join(CERT, "census"))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import scipy.sparse as sps
import flexexact8 as FX
import snapshot

SHA = {"cut9.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(), "flexexact8.py": FX.SHA}


def v8_records():
    out = {}
    for sub in ("out_flex8", "out_flex8r"):
        dd = os.path.join(CERT, "paper_v8", "recount", sub)
        for f in os.listdir(dd):
            if f.endswith(".json"):
                j = json.load(open(os.path.join(dd, f)))
                if j.get("status") == "FLEX":
                    out[(j["d"], j["n"])] = j
    return out


def v8_flags():
    path = os.path.join(CERT, "paper_v8", "v8_verdicts.tsv")
    rows = [l.rstrip("\n").split("\t") for l in open(path) if not l.startswith("#")]
    head = rows[0]
    return {(int(r[0]), int(r[1])): dict(zip(head, r)) for r in rows[1:]}


class Pairs:
    """The pair cosines of the double-precision code, sorted downwards (untrusted)."""

    def __init__(self, X, kmax):
        n = len(X)
        G = np.empty((n, n))
        for i0 in range(0, n, 256):                    # as flexexact8.proposal, so that the cut is the same
            G[i0:i0 + 256] = X[i0:i0 + 256] @ X.T
        np.fill_diagonal(G, -2.0)
        iu = np.triu_indices(n, 1)
        v = G[iu]
        kk = min(kmax + 1, len(v))
        part = np.argpartition(-v, kk - 1)[:kk]
        order = part[np.argsort(-v[part], kind="stable")]
        self.I, self.J, self.c = iu[0][order], iu[1][order], v[order]
        self.M = float(self.c[0])

    def k_of_window(self, w):
        """The number of pairs with cosine >= M - w."""
        return int(np.searchsorted(-self.c, -(self.M - w), side="right"))

    def snap(self, k, fac=10 ** 0.125):
        """The cut with the largest gap c_j - c_{j+1} among the cuts whose window lies within a factor fac of cut k."""
        if k >= len(self.c) - 1:
            return k
        w = self.M - self.c[k]
        lo, hi = self.k_of_window(w / fac), self.k_of_window(w * fac)
        lo, hi = max(1, lo), min(len(self.c) - 1, max(hi, lo + 1))
        gaps = self.c[lo - 1:hi - 1] - self.c[lo:hi]
        return lo + int(np.argmax(gaps)) if len(gaps) else k

    def window(self, k):
        return self.M - (float(self.c[k - 1]) + float(self.c[k])) / 2.0


def estimate(X, P, k, d, cap):
    """Untrusted double-precision beta, h, tau', delta for E = the top k pairs (normal-equations right inverse)."""
    I_, J_ = P.I[:k], P.J[:k]
    pts = np.unique(np.concatenate([I_, J_]))
    nE = len(pts)
    UL, UC = nE + k, nE * d
    if UL > cap or UL > UC:
        return None
    loc = -np.ones(len(X), int)
    loc[pts] = np.arange(nE)
    Xi = np.rint(np.ldexp(X, FX.GRID)).astype(np.int64)
    _, Lf = FX.build_L(Xi, X, pts, loc, I_, J_, d)
    Gm = (Lf @ Lf.T).toarray()
    try:
        R = np.linalg.inv(Gm)
    except np.linalg.LinAlgError:
        return None
    beta = float(np.abs(np.eye(UL) - Gm @ R).sum(1).max())
    LfT = Lf.T.tocsr()
    cn = np.zeros(UC)
    for j0 in range(0, UL, 512):
        cn += np.abs(np.asarray(LfT @ R[:, j0:j0 + 512])).sum(1)
    cnorm = float(cn.max())
    if not (beta < 0.9 and math.isfinite(cnorm) and cnorm > 0):
        return None
    h = (1 - beta) / (2 * cnorm)
    taup = P.M - float(P.c[k]) if k < len(P.c) else 2.0
    return {"k": k, "rows": UL, "cols": UC, "beta": beta, "h": h, "tau": taup,
            "delta": min(h / math.sqrt(d), taup / 4)}


def search(X, d, rec8, cap, log):
    kmax = 4 * cap
    P = Pairs(X, kmax)
    tested = {}

    def ev(k):
        k = max(1, min(int(k), len(P.c) - 1))
        if k not in tested:
            tested[k] = estimate(X, P, k, d, cap)
        return tested[k]

    def at_window(w):
        return ev(P.snap(max(1, P.k_of_window(w))))

    w0 = rec8["tau_prime"]
    best = ev(P.k_of_window(rec8["tol"]))        # the version 8 cut (approximately)
    if best is None:
        # no admissible estimate at the version 8 cut (a nearly singular L): fewer pairs first, then more
        for w2 in [w0 / 10 ** j for j in range(1, 5)] + [w0 * 10 ** j for j in range(1, 4)]:
            best = at_window(w2)
            if best is not None:
                break
    if best is None:
        return None, tested, P
    up = best["h"] / math.sqrt(d) > best["tau"] / 4
    w = P.M - float(P.c[best["k"]])
    for step in (10.0, 10 ** 0.5, 10 ** 0.25):
        while True:
            w2 = w * step if up else w / step
            if w2 <= 0 or w2 > 1.0:
                break
            e = at_window(w2)
            if e is None or e["delta"] <= best["delta"]:
                break
            best, w = e, w2
        # after a step that did not improve, try the other direction at the finer step
        for w2 in (w * step ** 0.5, w / step ** 0.5):
            e = at_window(w2)
            if e is not None and e["delta"] > best["delta"]:
                best, w = e, w2
    log["evaluations"] = len(tested)
    return best, tested, P


def main():
    out = sys.argv[1]
    args = sys.argv[2:]
    opt = lambda k, dflt: args[args.index(k) + 1] if k in args else dflt
    only = set(opt("--only", "").split(",")) if "--only" in args else None
    group = opt("--group", "all")
    cap = int(opt("--cap", "3100"))
    nsample = int(opt("--sample", "400"))
    recs = v8_records()
    flags = v8_flags()
    keys = sorted(recs, key=lambda k: recs[k].get("rows", 0))
    unfl = [k for k in keys if flags[k]["flagged"] == "0"]
    unres = [k for k in keys if flags[k]["verdict"] == "unresolved" and k not in set(unfl)]
    rest = [k for k in keys if k not in set(unfl) and k not in set(unres)]
    random.seed(20261108)
    samp = sorted(random.sample(rest, min(nsample, len(rest))), key=lambda k: recs[k].get("rows", 0))
    sset = set(samp)
    rest2 = [k for k in rest if k not in sset]
    todo = {"unflagged": unfl, "unresolved": unres, "sample": samp, "rest": rest2,
            "all": unfl + unres + samp + rest2}[group]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    os.makedirs(out, exist_ok=True)
    for d, n in todo:
        if only and f"{d}_{n}" not in only:
            continue
        js = os.path.join(out, f"d{d}_n{n}.json")
        lk = js[:-5] + ".lock"
        if os.path.exists(js):
            continue
        try:
            os.close(os.open(lk, os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        rec8 = recs[(d, n)]
        log = {"d": d, "n": n, "group": ("unflagged" if (d, n) in set(unfl) else
                                          "unresolved" if (d, n) in set(unres) else
                                          "sample" if (d, n) in sset else "rest"),
               "v8": {k: rec8.get(k) for k in ("tol", "rows", "columns", "contacts", "h", "tau_prime", "delta", "delta_by")},
               "program_sha256": SHA}
        t0, c0 = time.time(), time.process_time()
        try:
            pts = snapshot.read_points(SNAP, idx[(d, n)])
            dec = [[t_.strip() for t_ in pp] for pp in pts]
            X = np.array([[float(x) for x in pp] for pp in dec])
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            best, tested, P = search(X, d, rec8, cap, log)
            log["seconds_search"] = round(time.time() - t0, 1)
            if best is None:
                log["status"] = "no admissible cut"
            else:
                cands = sorted((e for e in tested.values() if e is not None), key=lambda e: -e["delta"])
                log["estimate"] = best
                for e in cands[:3]:
                    if e["delta"] <= rec8["delta"]:
                        break
                    w = P.window(e["k"])
                    a = {"contacts": e["k"], "points_in_contacts": e["rows"] - e["k"]}
                    r = FX.check(X, dec, d, n, w, a)
                    if r.get("status") == "FLEX" and r["delta"] > rec8["delta"]:
                        log.update(exact=r, window=w, k=e["k"], status="FLEX")
                        break
                    log.setdefault("exact_failures", []).append({"k": e["k"], "status": r.get("status")})
                if "status" not in log:
                    log["status"] = "no larger radius than version 8"
        except MemoryError:
            log["status"] = "error: memory"
        except Exception:
            log["status"] = "error"
            log["error"] = traceback.format_exc()[-1500:]
        log["wall"] = round(time.time() - t0, 1)
        log["cpu_seconds"] = round(time.process_time() - c0, 1)
        with open(js + ".part", "w") as f:
            json.dump(log, f)
        os.replace(js + ".part", js)
        os.remove(lk)
        ex = log.get("exact", {})
        print(d, n, log["status"], "v8 delta %.3g" % rec8["delta"], "new %.3g" % ex.get("delta", float("nan")),
              "k", log.get("k"), "rows", ex.get("rows"), log["cpu_seconds"], flush=True)


if __name__ == "__main__":
    main()
