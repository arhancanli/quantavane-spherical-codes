"""Every number, table and figure of the census sections of the paper (version 4), from the census outputs.

  python3 -I certify/census/summarize.py CENSUSDIR PAPERDIR

CENSUSDIR: certify/census (classify.jsonl, index.tsv, results_v4/, strict_v4/, exact_v4/, validate_v4.tsv,
snapshot_files.json). PAPERDIR: where v4_numbers.tex, v4_tab_*.tex and v4_fig_census.pdf are written. Also writes
CENSUSDIR/summary_v4.json (the status of every entry and every macro).
"""
import collections, glob, json, math, os, sys
from decimal import Decimal
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/mpllib"))
sys.path.insert(0, HERE)
import numpy as np
import run_census
import snapshot

RANGES = [(3, 5), (6, 10), (11, 17), (18, 24), (25, 32)]
COLS = ["strict_iso", "strict_hyper", "local", "flex", "failed", "too_large"]   # (no entry ends as nonopt_iso)
TOLS = [1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]
EXP = {1e-13: "-13", 1e-12: "-12", 1e-11: "-11", 1e-10: "-10", 1e-9: "-9", 1e-8: "-8", 1e-7: "-7"}
TW = {1e-13: "Thirteen", 1e-12: "Twelve", 1e-11: "Eleven", 1e-10: "Ten", 1e-9: "Nine", 1e-8: "Eight", 1e-7: "Seven"}


def sci(x, digits=1):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "--"
    if x == 0:
        return "$0$"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(m, digits) >= 10:
        m, e = m / 10, e + 1
    if 0 <= e <= 2:
        return f"${x:.{max(digits - e, 0)}f}$"
    return f"${m:.{digits}f}\\cdot10^{{{e}}}$"


def num(k):
    return f"{k:,}".replace(",", "{,}")


def word(name):
    """A macro name from a key: letters only."""
    out = ""
    up = True
    for ch in name:
        if ch.isalpha():
            out += ch.upper() if up else ch
            up = False
        else:
            up = True
    return out


def quant(vals, ps=(0, 0.1, 0.5, 0.9, 1)):
    v = sorted(vals)
    if not v:
        return [None] * len(ps)
    return [v[min(len(v) - 1, int(round(p * (len(v) - 1))))] for p in ps]


def main():
    cdir, pdir = sys.argv[1], sys.argv[2]
    os.makedirs(pdir, exist_ok=True)
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    recs = {(r["d"], r["n"]): r for r in map(json.loads, open(os.path.join(cdir, "classify.jsonl")))
            if idx[(r["d"], r["n"])]["listed"]}
    keys = sorted(recs)
    res, strict, exact = {}, {}, {}
    for f in glob.glob(os.path.join(cdir, "results_v4", "*.json")):
        j = json.load(open(f))
        res[(j["d"], j["n"])] = j
    for f in glob.glob(os.path.join(cdir, "strict_v4", "*.json")):
        j = json.load(open(f))
        strict[(j["d"], j["n"])] = j
    for f in glob.glob(os.path.join(cdir, "exact_v4", "strict", "*.json")):
        j = json.load(open(f))
        exact[(j["d"], j["n"])] = j
    M = {}
    M["NCensus"] = num(len(keys))
    # one version of the certifier (cq.py, sphcert.py) for every output, and of strict.py (stopped otherwise); the driver
    # run_census.py has two versions (the second adds the memory watchdog), recorded per output
    cert = lambda j: json.dumps({k: v for k, v in (j.get("code_sha256") or {}).items() if k != "run_census.py"},
                                sort_keys=True)
    vers = {cert(j) for j in res.values()}
    for j in res.values():
        for a in j["attempts"]:
            if a.get("code_sha256") is not None:
                vers.add(json.dumps(a["code_sha256"], sort_keys=True))
    svers = {json.dumps(j.get("code_sha256"), sort_keys=True) for j in strict.values()}
    if len(vers) > 1 or len(svers) > 1:
        raise SystemExit(f"outputs from more than one version of the certifier: {len(vers)} census, {len(svers)} strict")
    drivers = collections.Counter((j.get("code_sha256") or {}).get("run_census.py") for j in res.values())
    M["NDriverVersions"] = num(len(drivers))
    M["NWatchdogEntries"] = num(sum(1 for j in res.values() if j.get("mem_limit_mb", 0) > 0))
    M["NStoppedMemory"] = num(sum(1 for j in res.values() for a in j["attempts"]
                                  if str(a.get("status", "")).startswith("stopped")))
    M["NEntriesStoppedMemory"] = num(sum(1 for j in res.values() if any(str(a.get("status", "")).startswith("stopped")
                                                                         for a in j["attempts"])))
    M["NResults"] = num(sum(1 for k in keys if k in res))
    status, first = {}, {}
    for k in keys:
        j = res.get(k)
        if j is None:
            status[k] = "missing"
            continue
        st = j["status"]
        att = {a.get("status"): a for a in j["attempts"]}
        if st == "PASS":
            status[k], first[k] = "strict_iso", att["PASS"]
        elif st == "LOCALMIN":
            first[k] = att["LOCALMIN"]
            status[k] = "strict_hyper" if strict.get(k, {}).get("status") == "STRICT" else "local"
        elif st == "FLEX":
            status[k], first[k] = "flex", att["FLEX"]
        elif st == "NONOPT":
            status[k], first[k] = "nonopt_iso", att["NONOPT"]
        elif st == "not attempted":
            status[k] = "too_large"
        else:
            status[k] = "failed"
    cnt = collections.Counter(status.values())
    if cnt.get("nonopt_iso", 0):
        raise SystemExit("entries classed nonopt_iso: add the column back to the tables")
    for c in COLS + ["missing"]:
        M["N" + word(c)] = num(cnt.get(c, 0))
    nstrict = cnt["strict_iso"] + cnt["strict_hyper"]
    nlocal_all = cnt["strict_hyper"] + cnt["local"]
    M["NStrictAll"] = num(nstrict)
    M["NLocalAll"] = num(nlocal_all)
    M["NCertAll"] = num(nstrict + cnt["local"])
    M["NDiagAll"] = num(cnt["flex"] + cnt["nonopt_iso"])
    M["NNothing"] = num(cnt["failed"] + cnt["too_large"])
    M["NNothingFlagged"] = num(sum(1 for k in keys if status[k] in ("failed", "too_large") and idx[k]["flagged"]))
    M["PctCert"] = f"{100 * (nstrict + cnt['local']) / len(keys):.1f}"
    M["PctStrict"] = f"{100 * nstrict / len(keys):.1f}"
    # flagged / optimal splits
    for lab, sel in (("Flag", lambda k: idx[k]["flagged"] == 1), ("Unflag", lambda k: idx[k]["flagged"] == 0),
                     ("Opt", lambda k: idx[k]["optimal"] == 1)):
        ks = [k for k in keys if sel(k)]
        M[f"N{lab}All"] = num(len(ks))
        c = collections.Counter(status[k] for k in ks)
        for col in COLS:
            M[f"N{lab}{word(col)}"] = num(c.get(col, 0))
        M[f"N{lab}Cert"] = num(c["strict_iso"] + c["strict_hyper"] + c["local"])
        M[f"N{lab}Strict"] = num(c["strict_iso"] + c["strict_hyper"])
    opt = [k for k in keys if idx[k]["optimal"]]
    M["OptNotStrictList"] = ", ".join(f"$({d},{n})$" for d, n in opt if status[(d, n)] not in ("strict_iso", "strict_hyper"))
    M["OptStrictHyperList"] = ", ".join(f"$({d},{n})$" for d, n in opt if status[(d, n)] == "strict_hyper")
    # ------------------------------------------------------------------ Table 1 and the appendix table
    head = ["% written by certify/census/summarize.py (version 4); do not edit",
            "\\begin{tabular}{lrrrrrrr}", "\\toprule",
            " & & \\multicolumn{3}{c}{exact local minimum certified nearby} & diagnostic only"
            " & \\multicolumn{2}{c}{nothing certified} \\\\",
            "\\cmidrule(lr){3-5}\\cmidrule(lr){6-6}\\cmidrule(lr){7-8}",
            " & & \\multicolumn{2}{c}{strict} & strictness & contact set & certificate & above \\\\",
            "$d$ & entries & isostatic & hyperstatic & open & excluded & failed & size limit \\\\",
            "\\midrule"]
    lines = list(head)
    for lo, hi in RANGES + [(3, 32)]:
        ks = [k for k in keys if lo <= k[0] <= hi]
        c = collections.Counter(status[k] for k in ks)
        if (lo, hi) == (3, 32):
            lines.append("\\midrule")
        lab = "all" if (lo, hi) == (3, 32) else f"{lo}--{hi}"
        lines.append(f"{lab} & {num(len(ks))} & " + " & ".join(num(c.get(x, 0)) for x in COLS) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(pdir, "v4_tab_census.tex"), "w").write("\n".join(lines) + "\n")
    lines = ["% written by certify/census/summarize.py (version 4); do not edit", "\\begin{tabular}{rrrrrrrr}",
             "\\toprule", "$d$ & entries & strict iso. & strict hyper. & open & excluded & failed & too large \\\\",
             "\\midrule"]
    for d in range(3, 33):
        ks = [k for k in keys if k[0] == d]
        c = collections.Counter(status[k] for k in ks)
        lines.append(f"{d} & {num(len(ks))} & " + " & ".join(num(c.get(x, 0)) for x in COLS) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(pdir, "v4_tab_census_dims.tex"), "w").write("\n".join(lines) + "\n")
    # ------------------------------------------------------------------ the isostatic certificates
    passes = [k for k in keys if status[k] == "strict_iso"]
    A = {k: first[k] for k in passes}
    M["NInside"] = num(sum(1 for k in passes if A[k].get("code_inside_rho") is True))
    M["NOutside"] = num(sum(1 for k in passes if A[k].get("code_inside_rho") is False))
    M["NNoRho"] = num(sum(1 for k in passes if A[k].get("rho") is None))
    M["NRhoDefined"] = num(sum(1 for k in passes if A[k].get("rho") is not None))
    M["NUniqueCoversCode"] = num(sum(1 for k in passes if A[k]["r_unique"] > A[k]["dist_nonrattlers"]))
    M["NPassRattlers"] = num(sum(1 for k in passes if A[k]["rattlers"] > 0))
    M["NPassRattlersMoved"] = num(sum(1 for k in passes if (A[k].get("rattler_placement") or {}).get("moved_first_order", 0)
                                      + (A[k].get("rattler_placement") or {}).get("moved_second_order", 0)
                                      + (A[k].get("rattler_placement") or {}).get("moved_jointly", 0) > 0))
    tc = collections.Counter(A[k]["tol"] for k in passes)
    for t in TOLS:
        M["NPassTol" + TW[t]] = num(tc.get(t, 0))
    M["NPassTolTen"] = num(tc.get(1e-10, 0))
    M["NPassTolOther"] = num(len(passes) - tc.get(1e-10, 0))
    qrows = [("kappa_point", "sharpness $\\kappa/\\sqrt d$"), ("rho", "radius $\\rho$"),
             ("dist_ls", "distance $D(X_V)$"), ("dist_nonrattlers", "largest coordinate distance"),
             ("r_unique", "uniqueness radius $r_u$"), ("c_err", "radius of the enclosure of $c^*$"),
             ("krawczyk_beta", "largest row sum of (K)"), ("stress_beta", "$\\beta_S$ of (S)"),
             ("cond_inf_upper", "condition bound"), ("stress_min_rel", "smallest relative stress $|E|\\,\\omega_-$"),
             ("unknowns", "unknowns $U$")]
    lines = ["% written by certify/census/summarize.py (version 4); do not edit", "\\begin{tabular}{lrrrrr}",
             "\\toprule", " & min & 10\\% & median & 90\\% & max \\\\", "\\midrule"]
    for key, lab in qrows:
        vals = [float(A[k][key]) for k in passes if A[k].get(key) is not None]
        q = quant(vals)
        lines.append(lab + " & " + " & ".join(num(int(x)) if key == "unknowns" and x is not None else sci(x) for x in q)
                     + " \\\\")
        M["Q" + word(key) + "Med"] = sci(q[2])
        M["Q" + word(key) + "Min"] = sci(q[0])
        M["Q" + word(key) + "Max"] = sci(q[4])
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(pdir, "v4_tab_quant.tex"), "w").write("\n".join(lines) + "\n")
    nf = [k for k in passes if A[k]["tol"] in (1e-8, 1e-7) and A[k]["dist_ls"] > 1e-8]
    M["NearFlatN"] = num(len(nf))
    M["NearFlatDistMax"] = sci(max([A[k]["dist_ls"] for k in nf] or [0]))
    M["NearFlatCosMax"] = sci(max([abs(int(A[k]["c_mid"]) / 1e40 - A[k]["max_cos_double"]) for k in nf] or [0]))
    M["MaxUnknownsPass"] = num(max(A[k]["unknowns"] for k in passes)) if passes else "--"
    big = max(passes, key=lambda k: A[k]["unknowns"]) if passes else None
    if big:
        M["LargestPassD"], M["LargestPassN"] = str(big[0]), num(big[1])
        M["LargestPassCPU"] = f"{A[big]['cpu_seconds']:.0f}"
    M["MedCPUPass"] = f"{np.median([A[k]['cpu_seconds'] for k in passes]):.2f}" if passes else "--"
    # ------------------------------------------------------------------ hyperstatic
    locs = [k for k in keys if status[k] in ("strict_hyper", "local")]
    L = {k: first[k] for k in locs}
    M["LocTauMax"] = sci(max([L[k]["tau_extra"] for k in locs] or [0]))
    M["LocDispMax"] = sci(max([L[k]["localmin_displacement"] for k in locs] or [0]))
    ex = [L[k]["contacts"] - L[k]["isostatic_count"] for k in locs]
    M["LocExtraMed"] = num(int(np.median(ex))) if ex else "--"
    M["LocExtraMax"] = num(max(ex)) if ex else "--"
    M["LocExtraMin"] = num(min(ex)) if ex else "--"
    M["MaxUnknownsLoc"] = num(max([L[k]["unknowns"] for k in locs] or [0]))
    tl = collections.Counter(L[k]["tol"] for k in locs)
    M["NLocTolTen"] = num(tl.get(1e-10, 0))
    M["NLocTolOther"] = num(len(locs) - tl.get(1e-10, 0))
    # strictness outcomes
    sreason = collections.Counter()
    sdeg = collections.Counter()
    sdist = []
    for k in locs:
        s = strict.get(k)
        if s is None:
            sreason["not run"] += 1
            continue
        st = s["status"]
        if st == "STRICT":
            sdeg[s.get("field_degree")] += 1
            dl = s.get("dist_ls_upper")
            if dl is not None and not math.isfinite(dl) and s.get("dist_coord_upper") is not None:
                # the ball for the sum of squares straddles zero: use sqrt(n d) times the coordinate bound
                dl = math.sqrt(s["nV"] * k[0]) * s["dist_coord_upper"] * (1 + 1e-12)
            if dl is not None and math.isfinite(dl):
                sdist.append(dl)
            continue
        if st.startswith("no polynomial"):
            sreason["no polynomial of degree at most 10"] += 1
        elif st.startswith("cosine of degree"):
            sreason["cosine of degree above 12"] += 1
        elif "not recognised" in st:
            sreason["an inner product not in $\\mathbb Q(c)$"] += 1
        elif "too many" in st or "work limit" in st:
            sreason["too many distinct inner products or too large"] += 1
        else:
            sreason["other: " + st] += 1
    M["NStrictHyperFieldOne"] = num(sdeg.get(1, 0))
    M["NStrictHyperFieldTwo"] = num(sdeg.get(2, 0))
    M["NStrictHyperFieldHigh"] = num(sum(v for kk, v in sdeg.items() if kk and kk > 2))
    M["StrictDistMax"] = sci(max(sdist)) if sdist else "--"
    M["StrictDistMed"] = sci(float(np.median(sdist))) if sdist else "--"
    lines = ["% written by certify/census/summarize.py (version 4); do not edit", "\\begin{tabular}{lr}", "\\toprule",
             "hyperstatic entries certified LOCALMIN & " + num(len(locs)) + " \\\\", "\\midrule",
             "strict: inner products in $\\mathbb Q$ & " + num(sdeg.get(1, 0)) + " \\\\",
             "strict: inner products in a quadratic field & " + num(sdeg.get(2, 0)) + " \\\\",
             "strict: field of degree 3 to 12 & " + num(sum(v for kk, v in sdeg.items() if kk and kk > 2)) + " \\\\",
             "\\midrule"]
    for r_, v in sorted(sreason.items(), key=lambda x: -x[1]):
        lines.append("open: " + r_.replace("_", "\\_") + " & " + num(v) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(pdir, "v4_tab_strict.tex"), "w").write("\n".join(lines) + "\n")
    for r_, v in sreason.items():
        M["NOpen" + word(r_.split(":")[0])[:24]] = num(v)
    # hyperstatic failures
    hf = collections.Counter()
    for k in keys:
        if status[k] != "failed":
            continue
        hy = [a for a in res[k]["attempts"] if a.get("kind") == "hyperstatic"]
        if not hy:
            continue
        sts = " ".join(a.get("status", "") for a in hy)
        if "singular" in sts or "no square subsystem" in sts:
            hf["singular"] += 1
        elif "existence" in sts:
            hf["existence"] += 1
        elif "no positive stress" in sts:
            hf["nostress"] += 1
        else:
            hf["other"] += 1
    M["NHyperFail"] = num(sum(hf.values()))
    M["NHyperFailSingular"] = num(hf["singular"] + hf["existence"])
    M["NHyperFailNoStress"] = num(hf["nostress"])
    M["NHyperFailOther"] = num(hf["other"])
    M["NIsoFail"] = num(sum(1 for k in keys if status[k] == "failed" and any(a.get("kind") == "isostatic"
                                                                            for a in res[k]["attempts"])))
    # ------------------------------------------------------------------ flexibility
    flexes = [k for k in keys if status[k] == "flex"]
    F = {k: first[k] for k in flexes}
    fc = collections.Counter(F[k]["tol"] for k in flexes)
    for t in TOLS:
        M["NFlexTol" + TW[t]] = num(fc.get(t, 0))
    M["NFlexForest"] = num(sum(1 for k in flexes if F[k].get("forest")))
    M["NFlexQR"] = num(sum(1 for k in flexes if F[k].get("flex_right_inverse") == "QR"))
    rr = [F[k]["flex_radius_l1_stored"] for k in flexes]
    q = quant(rr)
    M["FlexRadMin"], M["FlexRadMed"], M["FlexRadMax"] = sci(q[0]), sci(q[2]), sci(q[4])
    M["NFlexSmallRad"] = num(sum(1 for x in rr if x < 1e-6))
    M["NFlexFlagged"] = num(sum(1 for k in flexes if idx[k]["flagged"]))
    M["NFlexUnflagged"] = num(sum(1 for k in flexes if not idx[k]["flagged"]))
    M["NFlexUnflaggedTopTol"] = num(sum(1 for k in flexes if not idx[k]["flagged"] and F[k]["tol"] == 1e-7))
    # iso NONOPT list
    nl = [k for k in keys if status[k] == "nonopt_iso"]
    M["NonoptList"] = ", ".join(f"$({d},{n})$ at $\\tau=10^{{{EXP[first[(d, n)]['tol']]}}}$" for d, n in nl)
    # ------------------------------------------------------------------ the tolerance table
    lines = ["% written by certify/census/summarize.py (version 4); do not edit", "\\begin{tabular}{lrrrrrrrr}",
             "\\toprule", " & \\multicolumn{4}{c}{class of the graph proposed at $\\tau$} &"
             " \\multicolumn{4}{c}{certificate obtained at $\\tau$} \\\\",
             "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
             "$\\tau$ & isostatic & hyperstatic & hypostatic & empty & strict (iso.) & local (hyper.) & of which strict"
             " & excluded \\\\", "\\midrule"]
    for t in TOLS:
        cl = collections.Counter(run_census.klass({x["tol"]: x for x in recs[k]["structs"]}[t]) for k in keys)
        npass = sum(1 for k in passes if A[k]["tol"] == t)
        nloc = sum(1 for k in locs if L[k]["tol"] == t)
        nstr = sum(1 for k in locs if L[k]["tol"] == t and status[k] == "strict_hyper")
        lines.append(f"$10^{{{EXP[t]}}}$ & " + " & ".join(num(cl.get(x, 0)) for x in ("iso", "hyper", "hypo", "empty"))
                     + f" & {num(npass)} & {num(nloc)} & {num(nstr)} & {num(fc.get(t, 0))} \\\\")
        for x in ("iso", "hyper", "hypo", "empty"):
            M[f"NClass{word(x)}{TW[t]}"] = num(cl.get(x, 0))
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(pdir, "v4_tab_tol.tex"), "w").write("\n".join(lines) + "\n")
    # entries whose class depends on the tolerance
    varying = 0
    for k in keys:
        s = {x["tol"]: x for x in recs[k]["structs"]}
        if len({run_census.klass(s[t]) for t in TOLS}) > 1:
            varying += 1
    M["NClassVaries"] = num(varying)
    # ------------------------------------------------------------------ the size cap
    sk = [k for k in keys if status[k] == "too_large"]
    M["NTooLargeHyperOnly"] = num(sum(1 for k in sk if all(a["kind"] == "hyper" for a in res[k]["skipped_above_cap"])))
    M["NTooLargeFlexOnly"] = num(sum(1 for k in sk if all(a["kind"] == "flex" for a in res[k]["skipped_above_cap"])))
    M["NTooLargeFlagged"] = num(sum(1 for k in sk if idx[k]["flagged"]))
    M["NTooLargeMixed"] = num(len(sk) - sum(1 for k in sk if all(a["kind"] == "hyper" for a in res[k]["skipped_above_cap"]))
                              - sum(1 for k in sk if all(a["kind"] == "flex" for a in res[k]["skipped_above_cap"])))
    M["SizeCap"] = num(3100)
    flexsizes = [F[k]["unknowns"] for k in flexes]
    M["MaxRowsFlex"] = num(max(flexsizes)) if flexsizes else "--"
    # ------------------------------------------------------------------ time
    cpu = collections.Counter()
    for k in keys:
        if k in res:
            for a in res[k]["attempts"]:
                cpu[a.get("kind", "?")] += a.get("cpu_seconds", 0)
    M["CPUHoursCensus"] = f"{sum(cpu.values()) / 3600:.1f}"
    M["CPUHoursStrict"] = f"{sum(s.get('cpu_seconds', 0) for s in strict.values()) / 3600:.1f}"
    # ------------------------------------------------------------------ class at 1e-10 (Section 5)
    s10 = lambda k: {x["tol"]: x for x in recs[k]["structs"]}[1e-10]
    cl10 = collections.Counter(run_census.klass(s10(k)) for k in keys)
    M["NGraph"] = num(len(keys) - cl10["empty"])
    M["NIsoTen"], M["NHyperTen"], M["NHypoTen"], M["NEmptyTen"] = (num(cl10[x]) for x in ("iso", "hyper", "hypo", "empty"))
    # rattlers among strict isostatic optima, by dimension range
    for lab, (lo, hi) in zip("ABCDE", RANGES):
        sh = [A[k]["rattlers"] / k[1] for k in passes if lo <= k[0] <= hi]
        M[f"RatShare{lab}"] = f"{100 * np.median(sh):.1f}\\%" if sh else "--"
        M[f"RatAny{lab}"] = f"{100 * np.mean([x > 0 for x in sh]):.0f}\\%" if sh else "--"
    # ------------------------------------------------------------------ validation against the exact cosines
    vt = os.path.join(cdir, "validate_v4.tsv")
    if os.path.exists(vt):
        rows = [l.rstrip("\n").split("\t") for l in open(vt)][1:]
        M["NExact"] = num(len(rows))
        tested = [r for r in rows if r[-1]]
        M["NExactTested"] = num(len(tested))
        M["NExactRootInside"] = num(sum(r[-1] == "root inside" for r in tested))
        M["NExactNoSign"] = num(sum(r[-1] == "NO SIGN CHANGE" for r in tested))
        M["ExactNoSignList"] = ", ".join(f"$({r[0]},{r[1]})$" for r in tested if r[-1] == "NO SIGN CHANGE")
    # displayed values: certified value rounded up to 12 decimals against the displayed cosine
    below = equal = above = 0
    above_gaps, above_units, above_keys = [], [], []
    for k in passes + locs:
        a = first[k]
        hi = Fraction(int(a["c_interval_30"][1]), 10 ** 30) + Fraction(a.get("tau_extra") or 0)
        c12 = Fraction(-(-hi.numerator * 10 ** 12 // hi.denominator), 10 ** 12)
        shown = Fraction(Decimal(idx[k]["cos"]))
        if c12 < shown:
            below += 1
        elif c12 == shown:
            equal += 1
        else:
            above += 1
            above_gaps.append(hi - shown)
            above_units.append((c12 - shown) * 10 ** 12)
            above_keys.append(k)
    M["NDispBelow"], M["NDispEqual"], M["NDispAbove"] = num(below), num(equal), num(above)
    M["DispAboveMaxGap"] = sci(float(max(above_gaps))) if above_gaps else "--"
    M["DispAboveMaxUnits"] = num(int(max(above_units))) if above_units else "--"
    M["NDispAboveHyper"] = num(sum(1 for k in above_keys if status[k] in ("local", "strict_hyper")))
    # ------------------------------------------------------------------ exact arithmetic checks
    ex_ok = sorted(k for k, j in exact.items() if j.get("status", "").startswith("STRICT"))
    M["NExactStrict"] = num(len(ex_ok))
    M["ExactStrictList"] = ", ".join(f"$({d},{n})$" for d, n in ex_ok)
    M["NExactTried"] = num(sum(1 for j in exact.values() if j.get("status") != "coordinates not recognised as exact"))
    # ------------------------------------------------------------------ write
    with open(os.path.join(pdir, "v4_numbers.tex"), "w") as f:
        f.write("% written by certify/census/summarize.py (version 4); do not edit\n")
        for k_, v in sorted(M.items()):
            f.write(f"\\newcommand{{\\{k_}}}{{{v}}}\n")
    json.dump({"status": {f"d{k[0]}_n{k[1]}": v for k, v in status.items()}, "macros": M},
              open(os.path.join(cdir, "summary_v4.json"), "w"), indent=0)
    # ------------------------------------------------------------------ figure
    import matplotlib
    matplotlib.use("pdf")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "font.family": "serif"})
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.1))
    exs = [s10(k)["m"] - s10(k)["iso"] for k in keys if s10(k)["nV"] > 0]
    ax[0].hist(np.clip(exs, -50, 200), bins=np.arange(-50.5, 201.5, 5), color="0.3")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("$|E|$ minus isostatic count ($\\tau=10^{-10}$)")
    ax[0].set_ylabel("entries")
    if passes:
        ax[1].hist(np.log10([float(A[k]["stress_min_rel"]) for k in passes]), bins=40, color="0.3")
        ax[1].set_xlabel("$\\log_{10}(|E|\\,\\omega_-)$, strict (isostatic)")
        ax[2].hist(np.log10([float(A[k]["dist_ls"]) for k in passes]), bins=40, color="0.3", alpha=0.8,
                   label="$D(X_V)$")
        ax[2].hist(np.log10([float(A[k]["rho"]) for k in passes if A[k].get("rho")]), bins=40, color="0.7",
                   alpha=0.8, label="$\\rho$")
        ax[2].legend(frameon=False)
        ax[2].set_xlabel("$\\log_{10}$ of distance and radius")
    fig.tight_layout()
    fig.savefig(os.path.join(pdir, "v4_fig_census.pdf"))
    print(json.dumps({k: M[k] for k in list(M)[:60]}, indent=0))


if __name__ == "__main__":
    main()
