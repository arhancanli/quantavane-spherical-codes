"""Exact certificate that a table entry is not a local minimum because its contact core is flat (paper, version 5,
proposition on flat cores).

  python3 -I certify/v5/subspace.py SNAPDIR D N OUT.json [--disc DISC]

TRUSTED: core.py builds the exact core configuration Q* near the stored code (exact Gram matrix of the core over
K = Q(sqrt DISC), of rank r, contacts E at the exact cosine c, rattlers below c in ball arithmetic, certified distance to
the stored code). Here, in addition:
 (F) c > 0 and r < d (the core is flat; r is proved exactly by core.py);
 (T) an explicit tilt: e = e_{r+1} Q is a rational unit vector orthogonal to the core of Q* (P* has zero coordinates
     r+1..d); a proper colouring k(a) of the core's contact graph (untrusted, checked); rational (cos, sin) pairs
     ((1-u^2)/(1+u^2), 2u/(1+u^2)) with u_a = s * slope(k(a)); the core points x_a(s) = cos_a x_a + sin_a e, rattlers
     fixed. Every inner product of the tilted configuration is shown to be below c: core-core pairs exactly in K
     (cos_a cos_b G_ab + sin_a sin_b), core-rattler pairs in ball arithmetic, rattler pairs as in core.py.
So M(tilted) < c = M(Q*): Q* is not a local minimum of the largest cosine (the tilt can be made as small as wanted).
"""
import argparse, hashlib, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import core
from kfield import K, ZERO, ONE, ksum, fq
import flint

SHA = {}
for _f in ("subspace.py", "core.py", "exact2.py", "kfield.py", "dist.py"):
    SHA[_f] = hashlib.sha256(open(os.path.join(HERE, _f), "rb").read()).hexdigest()


def pyth(u):
    u = fq(u)
    den = 1 + u * u
    return (1 - u * u) / den, 2 * u / den


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
    log = {"d": d, "n": n, "program": "subspace.py", "code_sha256": SHA, "disc": args.disc}
    C = core.build(args.snapdir, d, n, args.disc)
    if C["status"] != "ok":
        log.update(status=C["status"], **{k: v for k, v in C.get("log", {}).items()})
        json.dump(log, open(args.out, "w"), indent=1, default=str)
        print(json.dumps(log, default=str))
        return
    log.update(C["log"])
    G, V, E, c, r, Q, PQ, dec = C["G"], C["V"], C["E"], C["c"], C["r"], C["Q"], C["PQ"], C["dec"]
    nV = len(V)
    log.update(cos_positive=bool(c.sign() > 0), flat=bool(r < d))
    if not (c.sign() > 0 and r < d):
        log["status"] = "hypotheses not met (c <= 0 or the core spans R^d)"
        json.dump(log, open(args.out, "w"), indent=1, default=str)
        print(json.dumps({k: v for k, v in log.items() if k != "code_sha256"}, default=str))
        return
    # tilt direction: e = e_{r+1} Q (rational, unit, orthogonal to the core of Q*)
    e = [Q[r, k] for k in range(d)]
    if sum(x * x for x in e) != 1:
        raise ValueError("e not a unit vector")
    nb = {a: set() for a in range(nV)}
    for a, b in E:
        nb[a].add(b)
        nb[b].add(a)
    col = {}
    for v in sorted(range(nV), key=lambda v: -len(nb[v])):
        used = {col[u] for u in nb[v] if u in col}
        col[v] = next(k for k in range(nV + 1) if k not in used)
    if any(col[a] == col[b] for a, b in E):
        raise ValueError("colouring not proper")
    ncol = max(col.values()) + 1
    Kf = int(4 / float(c)) + 2
    slopes = [0]
    k = 0
    while len(slopes) < ncol:
        slopes += [Kf ** k, -Kf ** k]
        k += 1
    top = max(abs(x) for x in slopes)
    s = fq(1, 1000 * top * Kf)
    cs = [pyth(s * slopes[col[a]]) for a in range(nV)]
    log.update(colours=ncol, slope_base=Kf, tilt_s=str(s))
    nonneg, worst = 0, None
    for a in range(nV):
        ca, sa = cs[a]
        for b in range(a + 1, nV):
            cb, sb = cs[b]
            diff = G[a][b] * (ca * cb) + K(sa * sb) - c
            if diff.sign() >= 0:
                nonneg += 1
            if worst is None or diff > worst:
                worst = diff
    log.update(core_pairs_not_below_c=nonneg, core_max_minus_c=float(worst))
    # core-rattler pairs after the tilt (balls)
    old = flint.ctx.prec
    flint.ctx.prec = core.PREC
    Ds = flint.arb(args.disc).sqrt()
    cA = core.karb(c, Ds)
    ea = [flint.arb(x) for x in e]
    Vs = set(V)
    bad, worst_out = 0, None
    for rr in range(n):
        if rr in Vs:
            continue
        u = core.stored_unit(dec[rr])[0]
        eu = sum((x * y for x, y in zip(ea, u)), flint.arb(0))
        for a in range(nV):
            ca, sa = cs[a]
            ip = sum((x * y for x, y in zip(PQ[a], u)), flint.arb(0))
            dlt = flint.arb(ca) * ip + flint.arb(sa) * eu - cA
            up = float(dlt.upper())
            worst_out = up if worst_out is None else max(worst_out, up)
            if not (dlt < 0):
                bad += 1
    flint.ctx.prec = old
    log.update(tilted_rattler_pairs_not_below_c=bad, tilted_rattler_max_minus_c_upper=worst_out)
    ok = nonneg == 0 and bad == 0
    log["status"] = "NOT A LOCAL MINIMUM (flat core; explicit tilt below c)" if ok else "tilt check failed"
    log["seconds"] = round(time.time() - t0, 1)
    json.dump(log, open(args.out, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in log.items() if k != "code_sha256"}, default=str))


if __name__ == "__main__":
    main()
