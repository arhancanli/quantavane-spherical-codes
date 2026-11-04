"""Structure of the codes: contacts, rattlers, isostaticity, antipodality, and distance from the table's code.

  python3 -I structure.py > structure.out          (needs numpy; one thread is enough)
Codes measured: every run of finish.py (the active-set polish) whose summary line says polish=converged (the strict test:
LP decrease < 2e-15, KKT residual < 1e-10, no negative-curvature escape), from the finish lane (after 1 Nov 01:35)
and the fast lane (after 1 Nov 00:30); the run's own output file (data/runs/<lane>/), checked to have the summary's
largest cosine. Contact: a pair whose cosine (double precision Gram matrix) is within TOL = 1e-10 of the largest.
Rattler: a point with fewer than d contacts after the other rattlers are removed (repeatedly): on S^{d-1} such a point
can be moved off all its contacts (generically), so it carries no constraint. Isostatic count of the rest:
n_r(d-1) - d(d-1)/2 + 1 (the tangent degrees of freedom less rotations, plus the largest cosine itself). Previous code:
the table's code downloaded on 25-26 Oct (not in this repository), with the same point order (every stage moves points
in place).
Files not held here (see README): a code that did not beat the table of 25 Oct by 1e-10 is the table's code unchanged or
polished in the last digits, and the comparison needs the table's code itself. Their values are read from
data/measured/structure.tsv, unless table_codes.py has put matching copies in cache/, when they are recomputed and
checked against it. --measure writes that file (every file needed).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter
from decimal import Decimal
import numpy as np
from common import lane_runs, code, med, q, FAST_SPLIT, snapshot, output, table_code, is_listed, measured, \
    write_measured, same, note, MEASURE

FIN_SPLIT = 1793482500  # 2026-11-01 01:35 +04
TOL = 1e-10
snap = snapshot()


def load(path, d):
    X = np.array(code(path, d), dtype=float)
    return X / np.linalg.norm(X, axis=1, keepdims=True)


def analyse(X, d):
    n = len(X)
    G = X @ X.T
    np.fill_diagonal(G, -2.0)
    M = G.max()
    A = G >= M - TOL
    ncont = {t: int(np.triu(G >= M - t, 1).sum()) for t in (1e-12, 1e-10, 1e-8)}
    alive = np.ones(n, bool)
    while True:
        deg = (A & alive[None, :] & alive[:, None]).sum(1)
        drop = alive & (deg < d)
        if not drop.any():
            break
        alive &= ~drop
    deg0 = A.sum(1)
    nr = int(alive.sum())
    Er = int(np.triu(A & alive[None, :] & alive[:, None], 1).sum())
    iso = nr * (d - 1) - d * (d - 1) // 2 + 1 if nr else 0
    Gm = G.copy()
    np.fill_diagonal(Gm, 2.0)
    anti = int((Gm.min(1) < -1 + 1e-9).sum())
    return dict(M=M, ncont=ncont, E=ncont[1e-10], zero=int((deg0 == 0).sum()), rattlers=n - nr, nr=nr, Er=Er, iso=iso,
                anti=anti, degmin=int(deg0[alive].min()) if nr else 0, degmax=int(deg0.max()), A=A, alive=alive)


def displacement(X, Y):
    """Point-by-point distance between two codes with the same point order: as stored, and after the best rotation."""
    raw = np.linalg.norm(X - Y, axis=1)
    U, _, Vt = np.linalg.svd(X.T @ Y)
    Z = X @ (U @ Vt)
    al = np.linalg.norm(Z - Y, axis=1)
    return raw.max(), al.max(), np.sqrt((al ** 2).mean())


def compare(a, X, Y, n):
    """The comparison of our code X (analysed in a) with the table's code Y, as in the paper's scripts."""
    raw, al, rms = displacement(Y, X)
    keep = a["alive"]
    if keep.sum() >= X.shape[1]:
        _, alr, rmsr = displacement(Y[keep], X[keep])
    else:
        alr = rmsr = float("nan")
    G0 = Y @ Y.T
    np.fill_diagonal(G0, -2.0)
    M0 = G0.max()
    delta = np.sqrt(2 - 2 * a["M"])
    iu = np.triu_indices(n, 1)
    ours = set(zip(*np.nonzero(np.triu(a["A"], 1))))
    top = np.argsort(-G0[iu])[:len(ours)]
    prev_top = set(zip(iu[0][top], iu[1][top]))
    return dict(raw=raw, al=al, rms=rms, rel=al / delta, alr=alr, rmsr=rmsr, relr=alr / delta, M0=M0,
                prev_cont=int(np.triu(G0 >= M0 - TOL, 1).sum()), overlap=len(ours & prev_top) / len(ours),
                prev_near=int(np.triu(G0 >= M0 - 1e-6, 1).sum()))


F64 = ("M", "raw", "al", "rms", "rel", "alr", "rmsr", "relr", "M0")
INTS = ("nc12", "nc10", "nc8", "zero", "rattlers", "nr", "Er", "iso", "anti", "degmin", "degmax", "has_prev",
        "prev_cont", "prev_near")
A_COLS = ["M", "nc12", "nc10", "nc8", "zero", "rattlers", "nr", "Er", "iso", "anti", "degmin", "degmax"]
P_COLS = ["has_prev", "raw", "al", "rms", "rel", "alr", "rmsr", "relr", "M0", "prev_cont", "overlap", "prev_near"]
stored = measured("structure", {**{c: np.float64 for c in F64}, **{c: int for c in INTS}, "overlap": float})


def flat(a):
    x = {c: a[c] for c in A_COLS if c in a}
    x.update(nc12=a["ncont"][1e-12], nc10=a["ncont"][1e-10], nc8=a["ncont"][1e-8])
    return x


codes = {}
for lane, split in (("finish", FIN_SPLIT), ("fast", FAST_SPLIT)):
    for r in lane_runs(lane):
        if r["polish"] != "converged" or r["mtime"] < split or (lane == "finish" and r["secs"] >= 1700):
            continue
        k = (r["d"], r["n"])
        if k not in codes or r["checked"] < codes[k]["checked"]:
            r["lane"] = lane
            codes[k] = r
print(f"finish.py-converged codes: {len(codes)} entries (finish lane {sum(1 for r in codes.values() if r['lane'] == 'finish')},"
      f" fast lane {sum(1 for r in codes.values() if r['lane'] == 'fast')}), d {min(k[0] for k in codes)}-{max(k[0] for k in codes)}")
low = [r for k, r in codes.items() if snap[k][0] - r["checked"] < Decimal("1e-10")]
print(f"  of them beating the table of 25 Oct by >= 1e-10: {len(codes) - len(low)}; the other {len(low)} (fast lane "
      f"{sum(1 for r in low if r['lane'] == 'fast')}, finish lane {sum(1 for r in low if r['lane'] == 'finish')}) gained "
      f"{min(float(snap[(r['d'], r['n'])][0] - r['checked']) for r in low):.1e} to "
      f"{max(float(snap[(r['d'], r['n'])][0] - r['checked']) for r in low):.1e}")

res, out_rows, src, diff = [], [], Counter(), []
for k, r in sorted(codes.items()):
    d, n = k
    key = f"{d}_{n}"
    p = output(f"runs/{r['lane']}/{r['file']}")
    prev = table_code("needs_opt", r["file"])
    s = stored.get(key)
    if p is None:
        if MEASURE:
            sys.exit(f"--measure: runs/{r['lane']}/{r['file']} is missing (see README)")
        if s is None:
            continue
        a = {c: s[c] for c in A_COLS}
        a["ncont"] = {1e-12: s["nc12"], 1e-10: s["nc10"], 1e-8: s["nc8"]}
        a["E"] = a["ncont"][1e-10]
        src["code from data/measured"] += 1
    else:
        X = load(p, d)
        if X.shape != (n, d):
            continue
        a = analyse(X, d)
        if abs(a["M"] - float(r["checked"])) > 1e-14:
            continue  # the output file is not the code of this summary line
        src["code recomputed"] += 1
        if s is not None:
            diff += [(key, c) for c, v in flat(a).items() if not same(float(v) if c == "M" else v, s[c])]
    a["k"] = k
    if MEASURE and prev is None and is_listed("table_codes", "needs_opt/" + r["file"]):
        sys.exit(f"--measure: the table's code needs_opt/{r['file']} is not in cache/ (see README)")
    if p is not None and prev is not None:
        Y = load(prev, d)
        if Y.shape == X.shape:
            a.update(compare(a, X, Y, n))
            src["comparison recomputed"] += 1
            if s is not None and s["has_prev"]:
                diff += [(key, c) for c in P_COLS[1:] if not same(float(a[c]) if c in F64 else a[c], s[c])]
    elif s is not None and s["has_prev"]:
        a.update({c: s[c] for c in P_COLS[1:]})
        src["comparison from data/measured"] += 1
    if MEASURE:
        row = {"key": key, **flat(a), "has_prev": int("al" in a)}
        row.update({c: (a[c] if "al" in a else "") for c in P_COLS[1:]})
        out_rows.append(row)
    res.append(a)
print(f"measured {len(res)} (output file matches the summary line)")
if MEASURE:
    write_measured("structure", A_COLS + P_COLS, out_rows)
note(f"{dict(src)}; recomputed values that differ from data/measured/structure.tsv: {len(diff)} {diff[:5]}")

sep = sum(1 for a in res if a["ncont"][1e-12] == a["ncont"][1e-8])
print(f"\ncontact count the same at tolerance 1e-12 and 1e-8: {sep} of {len(res)}; "
      f"same at 1e-10 and 1e-8: {sum(1 for a in res if a['ncont'][1e-10] == a['ncont'][1e-8])}")
ex = [a["Er"] - a["iso"] for a in res]
c = Counter(ex)
print(f"contacts among non-rattlers minus isostatic count n_r(d-1)-d(d-1)/2+1: exactly 0 in {c[0]} of {len(res)}; "
      f"below 0 in {sum(1 for x in ex if x < 0)}; above 0 in {sum(1 for x in ex if x > 0)}")
print("  distribution of the excess:", sorted(c.items())[:12], "...", sorted(c.items())[-6:])
print(f"  among the {sum(1 for x in ex if x < 0)} below: median deficit {med([x for x in ex if x < 0])}")
for lo, hi in [(3, 10), (11, 17), (18, 24), (25, 32)]:
    s = [a for a in res if lo <= a["k"][0] <= hi]
    if s:
        print(f"  d {lo}-{hi}: {len(s)} codes, isostatic {sum(1 for a in s if a['Er'] == a['iso'])}, "
              f"with rattlers {sum(1 for a in s if a['rattlers'] > 0)}, median rattler share "
              f"{med([a['rattlers'] / a['k'][1] for a in s]):.3f}")
rat = [a["rattlers"] for a in res]
print(f"rattlers: codes with at least one {sum(1 for x in rat if x > 0)} of {len(res)}; median {med(rat)}, "
      f"max {max(rat)} (share max {max(a['rattlers'] / a['k'][1] for a in res):.3f}); points with no contact at all: "
      f"codes with one {sum(1 for a in res if a['zero'] > 0)}, total {sum(a['zero'] for a in res)}")
iso_r = [a for a in res if a["Er"] == a["iso"]]
print(f"  isostatic codes with rattlers: {sum(1 for a in iso_r if a['rattlers'] > 0)} of {len(iso_r)}")
mdeg = [2 * a["Er"] / a["nr"] / (2 * (a["k"][0] - 1)) for a in res if a["nr"]]
print(f"mean contacts per non-rattler point / 2(d-1): median {med(mdeg):.3f}, q10 {q(mdeg, .1):.3f}, q90 {q(mdeg, .9):.3f}")
print(f"min contacts of a non-rattler point: median {med([a['degmin'] for a in res])}; max contacts of a point: "
      f"median {med([a['degmax'] for a in res])}, largest {max(a['degmax'] for a in res)}")
an = [a for a in res if a["anti"] > 0]
print(f"antipodal pairs: codes with any {len(an)} of {len(res)}; fully antipodal "
      f"{sum(1 for a in res if a['anti'] == a['k'][1])}: {[a['k'] for a in an][:10]}")

wp = [a for a in res if "al" in a]
print(f"\ncompared with the table's code of 25-26 Oct (same point order): {len(wp)} codes")
for name, key in (("max move as stored", "raw"), ("max move after best rotation", "al"), ("rms move after rotation", "rms"),
                  ("max move / minimal distance", "rel"), ("non-rattlers: max move after rotation", "alr"),
                  ("non-rattlers: rms move", "rmsr"), ("non-rattlers: max move / minimal distance", "relr")):
    v = [a[key] for a in wp]
    print(f"  {name}: median {med(v):.2e}, q90 {q(v, .9):.2e}, max {max(v):.2e}")
ov = [a["overlap"] for a in wp]
print(f"  share of our contacts that are among the same number of closest pairs of the table's code: median {med(ov):.3f}, "
      f"q10 {q(ov, .1):.3f}, min {min(ov):.3f}; equal to 1 in {sum(1 for x in ov if x == 1)}")
pc = [a["prev_cont"] for a in wp]
print(f"  pairs within 1e-10 of the table code's own largest cosine: median {med(pc)}, ours {med([a['E'] for a in wp])}; "
      f"table codes with <= 2 such pairs: {sum(1 for x in pc if x <= 2)} of {len(wp)}")
print(f"  pairs within 1e-6 of the table code's largest cosine / our contacts: median "
      f"{med([a['prev_near'] / a['E'] for a in wp]):.3f}")
far = sorted(wp, key=lambda a: -a["rel"])[:5]
print("  largest relative moves:", [(a["k"], f"{a['rel']:.3f}") for a in far])
for lane in ("fast", "finish"):
    v = [a["relr"] for a in wp if codes[a["k"]]["lane"] == lane]
    print(f"  lane {lane}: {len(v)} codes, non-rattler max move / minimal distance: median {med(v):.2e}, q90 {q(v, .9):.2e},"
          f" max {max(v):.2e}; above 0.5: {sum(1 for x in v if x > 0.5)}")
print("  codes whose non-rattler points moved by more than half the minimal distance:",
      sorted((a["k"], round(a["relr"], 2)) for a in wp if a["relr"] > 0.5))
print("\nover-constrained codes (contacts above the isostatic count by more than 20):")
oc = sorted((a for a in res if a["Er"] - a["iso"] > 20), key=lambda a: a["k"])
print(f"  {len(oc)} codes; d: {dict(Counter(a['k'][0] for a in oc))}")
for a in oc[:40]:
    print(f"   d{a['k'][0]} N{a['k'][1]} M {a['M']:.12f} contacts {a['E']} isostatic {a['iso']} rattlers {a['rattlers']}"
          f" antipodal points {a['anti']}")
small = [a for a in res if 0 < a["Er"] - a["iso"] <= 20]
print(f"slightly over (1-20): {len(small)}: {[(a['k'], a['Er'] - a['iso']) for a in small]}")
under = [a for a in res if a["Er"] < a["iso"]]
print(f"under: {[(a['k'], a['Er'] - a['iso'], a['ncont']) for a in under]}")
sepc = [a for a in res if a["ncont"][1e-12] == a["ncont"][1e-8]]
print(f"\nwell-separated codes (same contact count at 1e-12 and 1e-8): {len(sepc)}; isostatic {sum(1 for a in sepc if a['Er'] == a['iso'])},"
      f" more {sum(1 for a in sepc if a['Er'] > a['iso'])}, fewer {sum(1 for a in sepc if a['Er'] < a['iso'])}")
for lo, hi in [(3, 10), (11, 17), (18, 24), (25, 32)]:
    s = [a for a in res if lo <= a["k"][0] <= hi]
    print(f"  d {lo}-{hi}: non-isostatic {sum(1 for a in s if a['Er'] != a['iso'])} of {len(s)} "
          f"({100 * sum(1 for a in s if a['Er'] != a['iso']) / len(s):.0f}%), more contacts {sum(1 for a in s if a['Er'] > a['iso'])}")
r_over = [a["Er"] / a["iso"] for a in res if a["Er"] > a["iso"]]
print(f"over-constrained: contacts / isostatic count median {med(r_over):.2f}, max {max(r_over):.2f}")
print(f"rattlers by d range (codes with any / codes, median share): " + "; ".join(
    f"d {lo}-{hi}: {sum(1 for a in res if lo <= a['k'][0] <= hi and a['rattlers'] > 0)}/{sum(1 for a in res if lo <= a['k'][0] <= hi)}"
    f" {med([a['rattlers'] / a['k'][1] for a in res if lo <= a['k'][0] <= hi]):.3f}" for lo, hi in [(3, 10), (11, 17), (18, 24), (25, 32)]))
