# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""Numbers and tables for the note, from the outputs of the exact programs (UNTRUSTED arithmetic: double precision and
mpmath are used here only for normalisations and displayed decimals; every identity and sign stated as proved in the
note comes from structure.gp and check.gp, which are exact).

  python3 -I certify/note/summarize.py NOTEDIR

Reads data/*.gp, data/stored_d22_n116.txt, out/aut116.json, out/construct.json, out/structure.json, out/check.json
and out/jam116_rerun.json. Computes the representative F of the flex orthogonal to the infinitesimal rotations (from the
exact first-order Gram change S and the stored code: F = X B + Z with Z = (I - Pi) S X Gamma^-1, B + B^T =
Gamma^-1 X^T S X Gamma^-1, Gamma B symmetric, Gamma = X^T X), its norm, the representative that fixes the orbit A,
how each orbit moves, and the second-order quantities normalised by |F|^2. Writes NOTEDIR/note_numbers.tex and the
tables NOTEDIR/tab_*.tex, and out/summary.json.
"""
import hashlib, json, math, os, re, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [p for p in os.environ.get("QV_LIBS", "").split(os.pathsep) if p]   # optional library folders
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.seterr(all="ignore")
import mpmath
from mpmath import mp, mpf

D = os.path.join(HERE, "data")
O = os.path.join(HERE, "out")


def pnum(x):
    """a real number printed by PARI ('1.25 E-17')"""
    return mpmath.mpf(str(x).replace(" E", "E").strip())


def fr(t):
    return Fraction(t.strip())


def gp_matrix(text, name):
    m = re.search(r"^%s = \[(.*?)\];" % name, text, re.M)
    return [[fr(x) for x in row.split(",")] for row in m.group(1).split(";")]


def gp_vec(text, name):
    m = re.search(r"^%s = \[(.*?)\];" % name, text, re.M)
    return [fr(x) for x in m.group(1).split(",")]


def gp_list(text, name):
    return [[int(x) for x in m.group(1).split(",")] for m in re.finditer(r"listput\(%s, \[([-0-9, ]+)\]\);" % name, text)]


def load():
    vt = open(os.path.join(D, "values.gp")).read()
    D1, D2, TC = gp_vec(vt, "D1"), gp_vec(vt, "D2"), gp_matrix(vt, "TC")
    ft = open(os.path.join(D, "flex.gp")).read()
    SV, SL = gp_matrix(ft, "SV"), gp_list(ft, "SL")
    st = open(os.path.join(D, "stress.gp")).read()
    OM, EL = gp_matrix(st, "OM"), gp_list(st, "EL")
    return D1, D2, TC, SV, SL, OM, EL


# ------------------------------------------------------------------------------------------------ exact formatting
def part(a, b):
    """(a + b sqrt145) as a LaTeX fraction; '' if zero."""
    if a == 0 and b == 0:
        return ""
    den = math.lcm(a.denominator, b.denominator)
    A, B = int(a * den), int(b * den)
    if B == 0:
        num = str(A)
    elif A == 0:
        num = ("-" if B < 0 else "") + ("" if abs(B) == 1 else str(abs(B))) + r"\sqrt{145}"
    else:
        num = f"{A}{'+' if B > 0 else '-'}{'' if abs(B) == 1 else abs(B)}" + r"\sqrt{145}"
    if den == 1:
        return num if (B == 0 or A == 0) else f"({num})"
    if A == 0 and B != 0 and B < 0:
        return r"-\frac{%s}{%d}" % (num[1:], den)
    if B == 0 and A < 0:
        return r"-\frac{%d}{%d}" % (-A, den)
    return r"\frac{%s}{%d}" % (num, den)


def tex_elt(c):
    """c: 8 tower coordinates (basis 1, r, s1, r s1, s2, r s2, s1 s2, r s1 s2)."""
    out = []
    for (i, j), mon in zip(((0, 1), (2, 3), (4, 5), (6, 7)), ("", "s_1", "s_2", "s_1s_2")):
        p = part(c[i], c[j])
        if p:
            out.append(p + (r"\,%s" % mon if mon else ""))
    if not out:
        return "0"
    s = out[0]
    for t in out[1:]:
        s += t if t.startswith("-") else "+" + t
    return s


def main():
    notedir = sys.argv[1]
    D1, D2, TC, SV, SL, OM, EL = load()
    aut = json.load(open(os.path.join(O, "aut116.json")))
    st = json.load(open(os.path.join(O, "structure.json")))["results"]
    ck = json.load(open(os.path.join(O, "check.json")))
    jr = json.load(open(os.path.join(O, "jam116_rerun.json")))
    cons = json.load(open(os.path.join(O, "construct.json")))
    assert ck["results"].get("STATUS") == "ALL CHECKS PASSED" and jr["status"].startswith("JAMMED")
    assert cons["model_equals_Xstar_classes"] and cons["Gamma_preserves_model_classes"]
    assert aut["order"] == aut["order_after_random_relabelling"] == 1440
    n, d = 116, 22
    mp.dps = 50
    r = mpmath.sqrt(145)
    s1 = -mpmath.sqrt(mpf(D1[0].numerator) / D1[0].denominator + mpf(D1[1].numerator) / D1[1].denominator * r)
    s2 = -mpmath.sqrt(mpf(D2[0].numerator) / D2[0].denominator + mpf(D2[1].numerator) / D2[1].denominator * r)
    bas = [mpf(1), r, s1, r * s1, s2, r * s2, s1 * s2, r * s1 * s2]
    ev = lambda c: mpmath.fsum(mpf(x.numerator) / x.denominator * b for x, b in zip(c, bas))
    tv = [ev(c) for c in TC]
    c = tv[12]
    # ---------------------------------------------------------------- the stored code (double) and S
    X = np.array([[float(t) for t in line.split(",")] for line in open(os.path.join(D, "stored_d22_n116.txt")).read().split()])
    X /= np.linalg.norm(X, axis=1)[:, None]
    S = np.zeros((n, n))
    for a, b, j, sg in SL:
        S[a, b] = S[b, a] = sg * float(ev(SV[j - 1]))
    Gam = X.T @ X
    Gi = np.linalg.inv(Gam)
    Pi = X @ Gi @ X.T
    Z = (np.eye(n) - Pi) @ S @ X @ Gi
    Sig = Gi @ X.T @ S @ X @ Gi
    # Gamma K + K Gamma = (Sig Gamma - Gamma Sig)/2, K skew, in the eigenbasis of Gamma
    lam, U = np.linalg.eigh(Gam)
    Rr = U.T @ ((Sig @ Gam - Gam @ Sig) / 2) @ U
    K = U @ (Rr / (lam[:, None] + lam[None, :])) @ U.T
    Bm = Sig / 2 + K
    F = X @ Bm + Z
    resS = np.abs(X @ F.T + F @ X.T - S).max()
    sym = np.abs(X.T @ F - (X.T @ F).T).max()
    nF2 = float((F ** 2).sum())
    names = aut["point_orbit_names"]
    A, B, P, Q = names["A"], names["B"], names["P"], names["Q"]
    # the representative fixing A: subtract the rotation Phi with Phi x_a = f_a on A (least squares, then residual)
    idx = [(p, q) for p in range(d) for q in range(p + 1, d)]
    Mrot = np.zeros((len(A) * d, len(idx)))
    for k, (p, q) in enumerate(idx):
        Ph = np.zeros((d, d)); Ph[p, q], Ph[q, p] = 1, -1
        Mrot[:, k] = (X[A] @ Ph.T).reshape(-1)
    coef = np.linalg.lstsq(Mrot, F[A].reshape(-1), rcond=None)[0]
    Phi = np.zeros((d, d))
    for k, (p, q) in enumerate(idx):
        Phi[p, q], Phi[q, p] = coef[k], -coef[k]
    FA = F - X @ Phi.T
    resA = float(np.abs(FA[A]).max())
    unit = 1 / math.sqrt(nF2)
    Fu, FAu = F * unit, FA * unit
    nrm = lambda Y, O_: [float(np.linalg.norm(Y[O_], axis=1).min()), float(np.linalg.norm(Y[O_], axis=1).max())]
    # U = span(P, Q): axis u0, normals nA, nB of span(A), span(B)
    u0 = X[P[0]] + X[P[1]]; u0 /= np.linalg.norm(u0)
    nA = np.linalg.svd(X[A])[2][-1]; nB = np.linalg.svd(X[B])[2][-1]
    Ub = np.vstack([u0, nA, nB])           # rows: an orthonormal basis of U (checked below)
    orthU = float(np.abs(Ub @ Ub.T - np.eye(3)).max())
    def turn(pair):
        """the infinitesimal rotation w (U-coordinates) with f_p = w x p for the two points (FAu)."""
        rows, rhs = [], []
        for v in pair:
            p = Ub @ X[v]
            Cx = np.array([[0, p[2], -p[1]], [-p[2], 0, p[0]], [p[1], -p[0], 0]])   # w x p = Cx w
            rows.append(Cx); rhs.append(Ub @ FAu[v])
        w = np.linalg.lstsq(np.vstack(rows), np.concatenate(rhs), rcond=None)[0]
        res = float(np.abs(np.vstack(rows) @ w - np.concatenate(rhs)).max())
        offU = max(float(np.linalg.norm(FAu[v] - Ub.T @ (Ub @ FAu[v]))) for v in pair)
        return w, res, offU
    wP, resP, offP = turn(P)
    wQ, resQ, offQ = turn(Q)
    # second-order numbers normalised by |F|^2 (F orthogonal to the rotations)
    normM = float(pnum(st["normF2_numeric"]))
    alpha_u = float(pnum(st["alpha_numeric"])) / nF2
    beta_u = float(pnum(st["beta_numeric"])) / nF2
    Hmm_u = float(pnum(st["H_maxmin_numeric"])) / nF2
    vert = []
    k = 1
    while f"vertex_{k}_H_over_normF2" in st:
        vert.append({"zero_types": st[f"vertex_{k}_zero_types"], "H_over_normF2": float(pnum(st[f"vertex_{k}_H_over_normF2"])) * normM / nF2,
                     "ratio": st.get(f"vertex_{k}_ratio_QQ_over_PP")})
        k += 1
    Hs = [v["H_over_normF2"] for v in vert]
    # geometry: <x, u0> on the orbits, angles
    ip_u0 = {nm: [float(X[O_] @ u0 @ np.ones(len(O_)) / len(O_))] for nm, O_ in names.items()}
    summary = {
        "unit_flex_residual_S": resS, "unit_flex_orthogonality_residual": sym, "norm2_F_exactS": nF2,
        "representative_fixing_A_residual": resA,
        "motion_unit_flex_orthogonal_to_rotations": {nm: nrm(Fu, O_) for nm, O_ in names.items()},
        "motion_representative_fixing_A_normalised": {nm: nrm(FAu, O_) for nm, O_ in names.items()},
        "norm_representative_fixing_A_over_unit": float(np.linalg.norm(FAu)),
        "U_basis_orthonormality_residual": orthU,
        "turn_P_axis_U_coords_(u0,nA,nB)": wP.tolist(), "turn_P_residual": resP, "turn_P_off_U": offP,
        "turn_Q_axis_U_coords_(u0,nA,nB)": wQ.tolist(), "turn_Q_residual": resQ, "turn_Q_off_U": offQ,
        "alpha_unit": alpha_u, "beta_unit": beta_u, "threshold": float(pnum(st["threshold_ratio_QQ_over_PP_numeric"])),
        "H_maxmin_over_normF2_total1": Hmm_u, "vertices": vert, "H_range_over_normF2": [min(Hs), max(Hs)],
        "ip_with_u0": ip_u0,
    }
    json.dump(summary, open(os.path.join(O, "summary.json"), "w"), indent=1)
    # ---------------------------------------------------------------- LaTeX: macros
    f6 = lambda x: f"{x:.6f}"
    def sci(x, k=2, up=True):
        """x in the form m.m...\cdot10^{-e} with k decimals, rounded up (upper bounds) or down"""
        e = math.floor(math.log10(x))
        m = x / 10 ** e
        m = math.ceil(m * 10 ** k - 1e-9) / 10 ** k if up else math.floor(m * 10 ** k + 1e-9) / 10 ** k
        return f"{m:.{k}f}" + r"\cdot10^{%d}" % e
    ckr = ck["results"]
    mac = {
        "Contacts": "4{,}042", "MinDeg": ckr["min_contacts_per_point"], "MaxDeg": ckr["max_contacts_per_point"], "RankUnknowns": ckr["rank_unknowns"],
        "RankEquations": f'{int(ckr["rank_equations"]):,}'.replace(",", "{,}"), "RankModp": ckr["rank_mod_p"],
        "PrimeCheck": f'{int(ckr["prime"]):,}'.replace(",", "{,}"),
        "StoredDev": sci(float(pnum(ckr["stored_code_max_deviation"])), 1),
        "ValueGap": f'{math.floor(float(pnum(ckr["smallest_gap_between_values"])) * 1e4) / 1e4:.4f}',
        "CheckSeconds": ckr["seconds_total"], "StructSeconds": st["seconds_total"], "PariVersion": ".".join(x.strip() for x in ck["pari_version"].strip('"[]').split(",")),
        "JamSeconds": round(jr["seconds"]), "JamPrime": f'{jr["prime"]:,}'.replace(",", "{,}"), "JamRank": f'{jr["rank_mod_p"]:,}'.replace(",", "{,}"),
        "JamDistLS": sci(jr["dist_ls_upper"], 1), "JamDistCoord": sci(jr["dist_coord_upper"], 1),
        "AlphaUnit": f"{alpha_u:.4f}", "BetaUnit": f"{beta_u:.4f}", "Threshold": f"{summary['threshold']:.4f}",
        "RatioMax": f"{max(float(pnum(v['ratio'])) for v in vert if v['ratio'] not in (None, 'PP = 0')):.3f}",
        "RatioMaxmin": f"{float(pnum(st['maxmin_ratio_QQ_over_PP'])):.3f}",
        "HMin": f"{min(Hs):.3e}".replace("e-0", r"\cdot10^{-") + "}", "HMax": f"{max(Hs):.3e}".replace("e-0", r"\cdot10^{-") + "}",
        "HMaxmin": f"{Hmm_u:.3e}".replace("e-0", r"\cdot10^{-") + "}",
        "MaxminMin": f'{float(pnum(st["maxmin_min_entry_times_E"])):.4f}',
        "FAB": f"{summary['motion_representative_fixing_A_normalised']['B'][0]:.3f}",
        "FAP": f"{summary['motion_representative_fixing_A_normalised']['P'][0]:.3f}",
        "FAQ": f"{summary['motion_representative_fixing_A_normalised']['Q'][0]:.3f}",
        "FUA": f"{summary['motion_unit_flex_orthogonal_to_rotations']['A'][0]:.4f}",
        "FUB": f"{summary['motion_unit_flex_orthogonal_to_rotations']['B'][0]:.4f}",
        "FUP": f"{summary['motion_unit_flex_orthogonal_to_rotations']['P'][0]:.4f}",
        "FUQ": f"{summary['motion_unit_flex_orthogonal_to_rotations']['Q'][0]:.4f}",
        "FANorm": f"{summary['norm_representative_fixing_A_over_unit']:.3f}",
    }
    # the minimal polynomial of s1 + s2 (check.gp), with integer coefficients
    from math import lcm
    pol = ckr["minimal_polynomial_of_s1_plus_s2"].replace(" ", "")
    terms = re.findall(r"([+-]?[0-9/]*)\*?(y(?:\^(\d+))?)?", pol)
    coef = {}
    for cf, mono, ex in terms:
        if not cf and not mono:
            continue
        e = (int(ex) if ex else 1) if mono else 0
        coef[e] = Fraction(cf + ("1" if cf in ("", "+", "-") else "")) if (cf in ("", "+", "-") or "/" in cf or cf.lstrip("+-").isdigit()) else None
    den = 1
    for v in coef.values():
        den = lcm(den, v.denominator)
    ic = {e: int(v * den) for e, v in coef.items()}
    def mono(e):
        return "" if e == 0 else "y" if e == 1 else "y^{%d}" % e
    ptex = ""
    for e in sorted(ic, reverse=True):
        v = ic[e]
        ptex += ("-" if v < 0 else ("+" if ptex else "")) + f"{abs(v):,}".replace(",", "{,}") + mono(e)
    mac["MinPoly"] = ptex
    mac["RealRoots"] = ckr["real_roots_of_P"]
    mac["AutBacktrack"] = ckr["automorphisms_by_backtracking"]
    mac["AutNodes"] = f'{int(ckr["backtracking_nodes"]):,}'.replace(",", "{,}")
    # degrees by orbit (exact counts of contacts)
    cm = [json.loads(line) for line in open(os.path.join(D, "class_matrix.gp")) if line.strip()]
    degs = {nm: sorted({sum(1 for v in row if v == 12) for row in (cm[a] for a in O_)}) for nm, O_ in names.items()}
    assert all(len(v) == 1 for v in degs.values())
    mac.update(DegA=degs["A"][0], DegB=degs["B"][0], DegP=degs["P"][0], DegQ=degs["Q"][0])
    mac["Basis"] = ", ".join(str(x) for x in json.loads(open(os.path.join(D, "basis.gp")).read()))
    # timings, rounded (they move between runs)
    about = lambda x: int(10 * round(float(x) / 10)) if float(x) >= 10 else int(round(float(x)))
    mac["JamSeconds"] = about(jr["seconds"])
    mac["CheckSeconds"] = about(ckr["seconds_total"])
    mac["StructSeconds"] = about(st["seconds_total"])
    mac["AutSeconds"] = about(aut["seconds"])
    def tower(key):
        return [Fraction(x.strip()) for x in st[key].strip("[]").split(",")]
    mac["AlphaExact"] = tex_elt(tower("alpha_tower"))
    mac["BetaExact"] = tex_elt(tower("beta_tower"))
    mac["ThresholdExact"] = tex_elt(tower("threshold_ratio_tower"))
    mac["HScaledExact"] = tex_elt(tower("H_scaled_min1_tower"))
    mac["HScaled"] = f'{float(pnum(st["H_scaled_min1_numeric"])):.4f}'
    mac["StressTotalScaled"] = f'{float(pnum(st["stress_scaled_total_weight_numeric"])):.1f}'
    with open(os.path.join(notedir, "note_numbers.tex"), "w") as f:
        f.write("% written by certify/note/summarize.py from the outputs of the exact programs; do not edit\n")
        for k_, v in mac.items():
            f.write("\\newcommand{\\%s}{%s}\n" % ("N" + k_, v))
    # ---------------------------------------------------------------- LaTeX: the values, with their orbits
    orbs = aut["orbitals"]
    fl_names = st["flex_orbit_names"].strip("[]").replace('"', "").split(", ")
    rows = []
    for kk in range(13):
        oo = [o for o in orbs if o["class"] == kk]
        desc = ", ".join(f"{o['orbits'][0]}{o['orbits'][1]}\\,{o['size']}{'$^*$' if o['split'] else ''}" for o in oo) if kk < 12 else "14 orbits (Table~\\ref{tab:stress})"
        rows.append(f"$t_{{{kk}}}$ & ${f6(float(tv[kk]))}$ & ${tex_elt(TC[kk])}$ & {sum(o['size'] for o in oo)} & {desc} \\\\")
    with open(os.path.join(notedir, "tab_values.tex"), "w") as f:
        f.write("% written by certify/note/summarize.py; do not edit\n\\begin{tabular}{llllp{4.6cm}}\n\\toprule\n"
                "& value & exact value & pairs & orbits of $\\mathrm{Aut}(X^*)$ on these pairs \\\\\n\\midrule\n")
        f.write("\n".join(rows).replace(" 4042 ", " 4{,}042 ") + "\n\\bottomrule\n\\end{tabular}\n")
    # ---------------------------------------------------------------- LaTeX: the flex
    frows = []
    for j, nm in enumerate(fl_names):
        cls_, rest = nm.split(":")
        val = SV[j]
        if all(x == 0 for x in val):
            continue
        frows.append(f"$t_{{{cls_}}}$ & {rest[0]}{rest[1]} & {rest[2:]} & ${tex_elt(val)}$ & ${float(ev(val)):+.5f}$ \\\\")
    with open(os.path.join(notedir, "tab_flex.tex"), "w") as f:
        f.write("% written by certify/note/summarize.py; do not edit\n\\begin{tabular}{lllll}\n\\toprule\n"
                "class & points & pairs & $\\pm$ value of $S$ on the two halves & numerically \\\\\n\\midrule\n")
        f.write("\n".join(frows) + "\n\\bottomrule\n\\end{tabular}\n")
    # ---------------------------------------------------------------- LaTeX: contact orbits and stresses
    cnames = st["contact_orbit_names"].strip("[]").replace('"', "").split(", ")
    sc = st["stress_scaled_min1_tower"]
    scv = [[Fraction(x.strip()) for x in grp.split(",")] for grp in re.findall(r"\[([^\[\]]+)\]", sc)]
    sigA = [Fraction(x.strip()) for x in st["sigma_A_on_types"].strip("[]").split(",")]
    srows = []
    for t, nm in enumerate(cnames):
        split = nm.endswith("s")
        base = nm[:-1] if split else nm
        srows.append(f"{t + 1} & {base[0]}{base[1]} & {base[2:]} & {'yes' if split else 'no'} & ${float(ev(scv[t])):.5f}$ & "
                     f"${(sigA[t] * 20).numerator if (sigA[t] * 20).denominator == 1 else sigA[t] * 20}$ \\\\")
    with open(os.path.join(notedir, "tab_stress.tex"), "w") as f:
        f.write("% written by certify/note/summarize.py; do not edit\n\\begin{tabular}{rlrlll}\n\\toprule\n"
                "& points & contacts & splits under $M_{10}$ & $\\omega^*$ (smallest entry 1) & $20\\,\\sigma_A$ \\\\\n\\midrule\n")
        f.write("\n".join(srows) + "\n\\bottomrule\n\\end{tabular}\n")
    erows = [f"{t + 1} & ${tex_elt(scv[t])}$ \\\\" for t in range(len(cnames))]
    with open(os.path.join(notedir, "tab_stress_exact.tex"), "w") as f:
        f.write("% written by certify/note/summarize.py; do not edit\n\\begin{tabular}{rl}\n\\toprule\n"
                "& $\\omega^*$ on the contact orbit (Table~\\ref{tab:stress}), exactly \\\\\n\\midrule\n")
        f.write("\n".join(erows) + "\n\\bottomrule\n\\end{tabular}\n")
    print(json.dumps({k_: v for k_, v in summary.items() if k_ not in ("vertices",)}, indent=None)[:3000])


if __name__ == "__main__":
    main()
