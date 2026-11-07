"""The numbers of the paper's sections on certification and contact graphs, from certify/results/ and the codes.

  python3 -I certify/analysis.py [OUTDIR] [--results DIR]   (needs numpy and scipy; writes OUTDIR/cert_numbers.tex, the tables
                                                OUTDIR/tab_*.tex, and prints every number with its definition)
Sets: "converged" = certify/lists/converged.txt, the 377 codes of observations/structure.py (finish.py-converged runs,
433 entries) whose files are in this repository; "uploaded" = codes/, the 1,042 codes we uploaded; "live" = the 447 of
them that the table showed on the page saved on 2 November 2026 at 18:46 (+04) (records/live_table.tsv, credit Canli).
A code counts as certified when certify/results/<set>/ or certify/results/<set>_retry/ holds a certificate with status
PASS (sphcert.py; the retry pass uses other contact tolerances). Statements about contact graphs that are not part of a
certificate (degrees, stresses of codes, colour refinement) are observations about the codes as stored.
"""
import csv, json, math, os, statistics, sys
from collections import Counter, defaultdict
from decimal import Decimal, ROUND_CEILING

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import numpy as np
import sphcert as S


def read_list(name):
    out = []
    for line in open(os.path.join(HERE, "lists", name + ".txt")):
        if line.strip() and not line.startswith("#"):
            p, d, n = line.split()
            out.append((p, int(d), int(n)))
    return out


RESULTS = os.path.join(HERE, "results")
SAME = {}  # uploaded path -> converged path of the same code (lists/uploaded_same.txt, written by dedup.py)
if os.path.exists(os.path.join(HERE, "lists", "uploaded_same.txt")):
    for line in open(os.path.join(HERE, "lists", "uploaded_same.txt")):
        a, b = line.split()
        SAME[a] = b


def result(setname, path):
    if setname == "uploaded" and path in SAME:
        r = result("converged", SAME[path])
        if r is not None:
            r = dict(r, same_as=SAME[path])
        return r
    name = os.path.basename(path).replace(".txt", "")
    main = os.path.join(RESULTS, setname, name + ".json")
    retry = os.path.join(RESULTS, setname + "_retry", name + ".json")
    r = json.load(open(main)) if os.path.exists(main) else None
    if r is not None and r["status"] != "PASS" and os.path.exists(retry):
        r2 = json.load(open(retry))
        r2["first_status"] = r["status"]
        if r2["status"] == "PASS":
            return r2
        r["retry_status"] = r2["status"]
    return r


def category(r):
    if r is None:
        return "not run"
    s = r["status"]
    if s == "PASS":
        return "certified" if r.get("tol", 1e-10) == 1e-10 else "certified (other tolerance)"
    if s.startswith("not square"):
        k = int(s.split("=")[1])
        if r.get("n_nonrattlers") == 0:
            return "no contact graph (every point a rattler)"
        return "more contacts than the isostatic count" if k > 0 else "fewer contacts than the isostatic count"
    if s.startswith("FAIL"):
        return s.lower()
    return s


def med(v):
    return statistics.median(v) if v else float("nan")


def ceil12(x):
    return Decimal(x).quantize(Decimal("1e-12"), rounding=ROUND_CEILING)


def drange(d):
    return "3-10" if d <= 10 else "11-17" if d <= 17 else "18-24" if d <= 24 else "25-32"


def colour_refinement(nv, edges):
    """Number of colour classes of the stable colouring (1-dimensional Weisfeiler-Leman) of the graph."""
    adj = [[] for _ in range(nv)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    col = [len(a) for a in adj]
    k = len(set(col))
    while True:
        S.PACE.tick()
        sig = []
        for v in range(nv):
            if v % 64 == 0:
                S.PACE.tick()
            sig.append((col[v], tuple(sorted(col[u] for u in adj[v]))))
        lab = {s: i for i, s in enumerate(sorted(set(sig)))}
        new = [lab[s] for s in sig]
        k2 = len(lab)
        if k2 == k:
            return k2
        col, k = new, k2


def main():
    global RESULTS
    args = sys.argv[1:]
    if "--results" in args:
        i = args.index("--results")
        RESULTS = args[i + 1]
        del args[i:i + 2]
    out = args[0] if args else os.path.join(HERE, "tables")
    os.makedirs(out, exist_ok=True)
    macros = {}
    P = print

    # ------------------------------------------------------------------ the live set and the uploaded values
    live = {}
    for line in open(os.path.join(ROOT, "records", "live_table.tsv"), encoding="utf-8"):
        if line.startswith("#"):
            continue
        r = line.rstrip("\n").split("\t")
        if len(r) > 3 and "Canli" in r[3]:
            live[(int(r[0]), int(r[1]))] = Decimal(r[2])
    manifest = {(int(r["d"]), int(r["n"])): Decimal(r["uploaded_max_cos"])
                for r in csv.DictReader(open(os.path.join(ROOT, "codes", "MANIFEST.tsv")), delimiter="\t")}

    sets = {"converged": read_list("converged"), "uploaded": read_list("uploaded")}
    macros["NUploadedSameAsConverged"] = len(SAME)
    res = {s: {(d, n): (p, result(s, p)) for p, d, n in lst} for s, lst in sets.items()}
    res["live"] = {k: v for k, v in res["uploaded"].items() if k in live}

    P("== certification outcomes")
    for s in ("converged", "uploaded", "live"):
        cats = Counter(category(r) for _, r in res[s].values())
        P(f"{s}: {len(res[s])} codes: " + "; ".join(f"{k} {v}" for k, v in sorted(cats.items(), key=lambda x: -x[1])))
        cert = [r for _, r in res[s].values() if r and r["status"] == "PASS"]
        macros[f"NCert{s.capitalize()}"] = len(cert)
        macros[f"NAll{s.capitalize()}"] = len(res[s])
        macros[f"NCertOtherTol{s.capitalize()}"] = sum(1 for r in cert if r.get("tol", 1e-10) != 1e-10)
        notsq = [r for _, r in res[s].values() if r and r["status"].startswith("not square")]
        notsq_ne = [r for r in notsq if r["n_nonrattlers"] > 0]
        macros[f"NOver{s.capitalize()}"] = sum(1 for r in notsq_ne if int(r["status"].split("=")[1]) > 0)
        macros[f"NUnder{s.capitalize()}"] = sum(1 for r in notsq_ne if int(r["status"].split("=")[1]) < 0)
        macros[f"NOverNoTol{s.capitalize()}"] = sum(1 for r in notsq_ne if int(r["status"].split("=")[1]) > 0 and
                                                   r.get("retry_status", "").startswith("no tolerance"))
        macros[f"NUnderNoTol{s.capitalize()}"] = sum(1 for r in notsq_ne if int(r["status"].split("=")[1]) < 0 and
                                                    r.get("retry_status", "").startswith("no tolerance"))
        empty = [r for r in notsq if r["n_nonrattlers"] == 0]
        macros[f"NEmpty{s.capitalize()}"] = len(empty)
        P(f"   of the codes not square: no point left after removing rattlers (fewer than d pairs per point within the "
          f"tolerance): {len(empty)}")
        P("   retry statuses of the codes not certified:", Counter(r.get("retry_status", "no retry result")
                                                               for _, r in res[s].values() if r and r["status"] != "PASS"))
        fails = [r for _, r in res[s].values() if r and r["status"].startswith("FAIL")]
        macros[f"NFail{s.capitalize()}"] = len(fails)
        if fails:
            P("   failures:", Counter(r["status"] for r in fails))
        if not cert:
            continue
        dist = [r["dist_nonrattlers"] for r in cert]
        drat = [r["dist_rattlers"] for r in cert if r["rattlers"]]
        cerr = [r["c_err"] for r in cert]
        smin = [r["stress_min"] * r["contacts"] for r in cert]  # smallest stress / mean stress
        secs = [r["seconds"]["total"] for r in cert]
        gaps = [r["gaps"]["min_gap_float"] for r in cert]
        P(f"   certified {len(cert)} (other tolerance {macros[f'NCertOtherTol{s.capitalize()}']}); distance of the "
          f"non-rattlers to the code: median {med(dist):.1e}, max {max(dist):.1e}; rattlers moved: max {max(drat or [0]):.1e};"
          f" radius of c*: median {med(cerr):.1e}, max {max(cerr):.1e}; smallest stress / mean stress: median "
          f"{med(smin):.2e}, min {min(smin):.2e}; smallest certified gap below c*: median {med(gaps):.1e}, min "
          f"{min(gaps):.1e}; seconds (paced): median {med(secs):.1f}, max {max(secs):.1f}")
        macros[f"MaxDist{s.capitalize()}"] = f"{max(dist):.0e}"

    # ------------------------------------------------------------------ the converged runs whose files are not here
    unpub = [json.load(open(os.path.join(RESULTS, "unpublished", f))) for f in
             sorted(os.listdir(os.path.join(RESULTS, "unpublished")))] if os.path.isdir(
        os.path.join(RESULTS, "unpublished")) else []
    macros["NUnpubAll"] = len(unpub)
    macros["NUnpubCert"] = sum(1 for r in unpub if r["status"] == "PASS")
    big = [r for _, r in res["uploaded"].values() if r and r["d"] * r["n"] > 6000]
    macros["NBigUploaded"] = len(big)
    macros["NBigEmpty"] = sum(1 for r in big if r.get("n_nonrattlers") == 0)
    P(f"uploaded codes with more than 6,000 coordinates: {len(big)}, with no contact graph: {macros['NBigEmpty']}")
    P(f"unpublished isostatic converged codes (files listed in observations/data/unpublished_outputs.tsv): "
      f"{len(unpub)}, certified {macros['NUnpubCert']}: " + ", ".join(f"{r['d']}/{r['n']} {r['status']}" for r in unpub))

    # ------------------------------------------------------------------ over all certificates (both sets, distinct files)
    allc = {}
    for s in ("converged", "uploaded"):
        for k, (p, r) in res[s].items():
            if r and r["status"] == "PASS":
                allc[r.get("same_as", p)] = r  # an uploaded code identical to a converged one counts once
    allc = list(allc.values())

    def sci(x):
        m, e = f"{x:.1e}".split("e")
        return f"${m}\\cdot10^{{{int(e)}}}$"

    if allc:
        dist = sorted(r["dist_nonrattlers"] for r in allc)
        cerr = [r["c_err"] for r in allc]
        cpu = [r["cpu_seconds"] for r in allc if "cpu_seconds" in r]
        unknowns = [r["n_nonrattlers"] * r["d"] + 1 for r in allc]
        moved = [r for r in allc if r.get("rattler_placement", {}).get("moved_first_order", 0)
                 + r.get("rattler_placement", {}).get("moved_second_order", 0) > 0]
        second = [r for r in allc if r.get("rattler_placement", {}).get("moved_second_order", 0) > 0]
        smin = [r["stress_min"] * r["contacts"] for r in allc]
        serr = [r["stress_err_max"] * r["contacts"] for r in allc]
        macros["MaxStressErrRel"] = sci(max(serr))
        margin = [r["stress_min"] / r["stress_err_max"] for r in allc]  # smallest stress over the largest error bound
        macros["MedStressMargin"] = sci(med(margin))
        macros["MinStressMargin"] = f"{min(margin):.1f}"
        macros["NStressMarginBelowOne"] = sum(1 for x in margin if x < 1)
        macros["NStressMarginBelowOneWord"] = ["no", "one", "two", "three", "four", "five"][min(5, macros["NStressMarginBelowOne"])]
        P(f"   smallest stress / largest stress error bound, per code: median {med(margin):.2e}, min {min(margin):.2f}")
        exact_used = [r for r in allc if r["gaps"]["decided_exactly"] > 0]
        macros.update(NCertDistinct=len(allc), MedDist=sci(med(dist)), MaxDist=sci(max(dist)),
                      NDistBelowTen=sum(1 for x in dist if x < 1e-10), MaxCerr=sci(max(cerr)), MedCerr=sci(med(cerr)),
                      MaxUnknowns=f"{max(unknowns):,}".replace(",", "{,}"), NRatMoved=len(moved),
                      NRatSecond=len(second), MinStressRel=sci(min(smin)), MedStressRel=sci(med(smin)),
                      NGapExact=len(exact_used), NWithRattlers=sum(1 for r in allc if r["rattlers"] > 0))
        if cpu:
            macros.update(MedCPU=f"{med(cpu):.1f}", MaxCPU=f"{max(cpu):.0f}", NCPU=len(cpu))
        P(f"\n== all {len(allc)} distinct certified files: distance median {med(dist):.1e} max {max(dist):.1e} (below 1e-10: "
          f"{macros['NDistBelowTen']}); c* radius median {med(cerr):.1e} max {max(cerr):.1e}; unknowns up to "
          f"{max(unknowns)}; rattlers moved in {len(moved)} (second order in {len(second)}), codes with rattlers "
          f"{macros['NWithRattlers']}; smallest stress / mean: median {med(smin):.2e} min {min(smin):.2e}; exact gap "
          f"decisions used in {len(exact_used)}; CPU seconds median {med(cpu):.2f} max {max(cpu or [0]):.1f} ({len(cpu)} with CPU)")

    # ------------------------------------------------------------------ c* against our uploaded value and the table
    P("\n== live entries: the certified c* against our uploaded value and the displayed value")
    below, same12, lower12 = [], 0, []
    for k, (p, r) in sorted(res["live"].items()):
        if not r or r["status"] != "PASS":
            continue
        lo, hi = (Decimal(x) for x in r["c_interval_25"])
        up_v = manifest[k]
        if hi < up_v - Decimal("5e-21"):  # the manifest keeps 20 decimals of the checker's 50-digit value, rounded
            below.append(float(up_v - hi))
        if ceil12(hi) == live[k] and ceil12(lo) == live[k]:
            same12 += 1
        elif ceil12(hi) < live[k]:
            lower12.append((k, str(live[k]), r["c_interval_25"][1]))
    P(f"certified live entries whose c* is below our uploaded value: {len(below)}; median gap {med(below):.1e}, max "
      f"{max(below or [0]):.1e}; c* rounds up to the displayed 12 decimals in {same12}; lower at 12 decimals in "
      f"{len(lower12)}: {lower12[:10]}")
    macros["NLiveSameTwelve"] = same12
    macros["NLiveLowerTwelve"] = len(lower12)
    macros["NLiveBelowUploaded"] = len(below)
    macros["MaxLiveBelowUploaded"] = sci(max(below)) if below else "0"
    macros["MedLiveBelowUploaded"] = sci(med(below)) if below else "0"

    # ------------------------------------------------------------------ examples: enclosures of c*
    cand = [(6, 248), (8, 448), (32, 232), (3, 952), (4, 310), (5, 267), (7, 118), (16, 63), (21, 54), (12, 181),
            (19, 248), (10, 154)]
    ex = []
    for k in cand:
        for s_ in ("uploaded", "converged"):
            p_r = res[s_].get(k)
            if p_r and p_r[1] and p_r[1]["status"] == "PASS":
                ex.append((k, s_, p_r[1]))
                break
        if len(ex) == 6:
            break
    with open(os.path.join(out, "tab_examples.tex"), "w") as f:
        f.write("% written by certify/analysis.py; do not edit\n\\begin{tabular}{rrrrrlr}\n\\toprule\n")
        f.write("$d$ & $N$ & $|V|$ & $|E|$ & code & $c^*$ (enclosure, 25 decimals) & distance \\\\\n\\midrule\n")
        for (d, n), s_, r in sorted(ex):
            lo, hi = r["c_interval_25"]
            common = os.path.commonprefix([lo, hi])
            f.write(f"{d} & {n} & {r['n_nonrattlers']} & {r['contacts']} & {'served' if (d, n) in live and s_ == 'uploaded' else s_} & "
                    f"\\texttt{{{lo}}} & {sci(r['dist_nonrattlers'])} \\\\\n")
            f.write(f" & & & & & \\texttt{{{hi}}} & \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    P("\n== examples:", [(k, s_, r["c_interval_25"]) for k, s_, r in ex])

    # ------------------------------------------------------------------ by dimension
    rows = []
    for lo_d, hi_d in ((3, 10), (11, 17), (18, 24), (25, 32)):
        row = [f"{lo_d}--{hi_d}"]
        for s in ("converged", "uploaded", "live"):
            ks = [k for k in res[s] if lo_d <= k[0] <= hi_d]
            c = sum(1 for k in ks if res[s][k][1] and res[s][k][1]["status"] == "PASS")
            row += [len(ks), c]
        rows.append(row)
    tot = ["all"] + [sum(r[i] for r in rows) for i in range(1, 7)]
    with open(os.path.join(out, "tab_cert.tex"), "w") as f:
        f.write("% written by certify/analysis.py; do not edit\n\\begin{tabular}{lrrrrrr}\n\\toprule\n")
        f.write(" & \\multicolumn{2}{c}{converged codes} & \\multicolumn{2}{c}{uploaded codes} & "
                "\\multicolumn{2}{c}{served by the table} \\\\\n\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\n")
        f.write("$d$ & codes & certified & codes & certified & codes & certified \\\\\n\\midrule\n")
        for r in rows + [tot]:
            if r[0] == "all":
                f.write("\\midrule\n")
            f.write(" & ".join(str(x) for x in r) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    P("\n== by dimension (converged codes, certified; uploaded, certified; live, certified)")
    for r in rows + [tot]:
        P("  ", r)

    # ------------------------------------------------------------------ contact graphs of the converged codes
    P("\n== contact graphs of the converged codes (contacts within 1e-10 of the largest cosine, double precision)")
    cls = defaultdict(Counter)
    cls_nd = defaultdict(Counter)
    deg_stats, rat_share, rat_any = [], defaultdict(list), defaultdict(list)
    exact_d, refine_disc, refine_tot = [], 0, 0
    stress_rows = []
    over_refine = []
    stable, stable14, antipodal, cert_disc = defaultdict(list), defaultdict(list), [], []
    for (d, n), (p, r) in sorted(res["converged"].items()):
        S.PACE.tick()
        X = S.read_code(os.path.join(ROOT, p), d, n)
        V, E, M, G = S.contact_structure(X, d, 1e-10)
        S.PACE.tick()
        nV, m = len(V), len(E)
        iso = nV * (d - 1) - d * (d - 1) // 2 + 1
        kind = "isostatic" if m == iso else "more" if m > iso else "fewer"
        iu = np.triu_indices(n, 1)
        gaps_all = M - G[iu]
        stable[kind].append(int((gaps_all <= 1e-12).sum()) == int((gaps_all <= 1e-8).sum()))
        stable14[kind].append(int((gaps_all <= 1e-14).sum()) == int((gaps_all <= 1e-9).sum()))
        Gm = G.copy()
        np.fill_diagonal(Gm, 2.0)
        if Gm.min() < -1 + 1e-9:
            antipodal.append((d, n))
        cls[drange(d)][kind] += 1
        q = n / d
        cls_nd["N/d<=8" if q <= 8 else "8<N/d<=16" if q <= 16 else "16<N/d<=32" if q <= 32 else "N/d>32"][kind] += 1
        rat_share[drange(d)].append((n - nV) / n)
        rat_any[drange(d)].append(n - nV > 0)
        deg = np.zeros(nV, int)
        for a, b in E:
            deg[a] += 1
            deg[b] += 1
        k = colour_refinement(nV, E)
        S.PACE.tick()
        refine_tot += 1
        refine_disc += k == nV
        if r and r["status"] == "PASS":
            cert_disc.append(k == nV)
        if kind == "more":
            over_refine.append(k == nV)
        if kind == "isostatic":
            deg_stats.append((d, nV, int(deg.min()), int(deg.max()), float(deg.mean()), float((deg == d).mean())))
            exact_d.append((deg == d).mean())
            if r and r["status"] == "PASS" and "stress_quantiles_over_mean" in r:
                stress_rows.append((d, r["stress_quantiles_over_mean"]))
    tot_cls = Counter()
    for c in cls.values():
        tot_cls.update(c)
    macros.update(NIsoConv=tot_cls["isostatic"], NMoreConv=tot_cls["more"], NFewerConv=tot_cls["fewer"])
    macros["NIsoAllRuns"] = tot_cls["isostatic"] + macros["NUnpubAll"]
    iso_cert = sum(1 for (d, n), (p, r) in res["converged"].items() if r and r["status"] == "PASS" and r.get("tol", 1e-10) == 1e-10)
    macros["NIsoCertConv"] = iso_cert
    macros["NIsoAllCert"] = iso_cert + macros["NUnpubCert"]
    for key in ("3-10", "11-17", "18-24", "25-32"):
        c = cls[key]
        P(f"  d {key}: {sum(c.values())} codes: isostatic {c['isostatic']}, more contacts {c['more']}, fewer {c['fewer']};"
          f" with rattlers {sum(rat_any[key])}, median rattler share {med(rat_share[key]):.3f}")
    for key in ("N/d<=8", "8<N/d<=16", "16<N/d<=32", "N/d>32"):
        c = cls_nd[key]
        P(f"  {key}: {sum(c.values())} codes: isostatic {c['isostatic']}, more {c['more']}, fewer {c['fewer']}")
    with open(os.path.join(out, "tab_iso.tex"), "w") as f:
        f.write("% written by certify/analysis.py; do not edit\n\\begin{tabular}{lrrrrr}\n\\toprule\n")
        f.write("$d$ & codes & isostatic & more contacts & fewer contacts & with rattlers \\\\\n\\midrule\n")
        for key in ("3-10", "11-17", "18-24", "25-32"):
            c = cls[key]
            f.write(f"{key.replace('-', '--')} & {sum(c.values())} & {c['isostatic']} & {c['more']} & {c['fewer']} & "
                    f"{sum(rat_any[key])} \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    macros.update(NIsoMinDegD=sum(1 for x in deg_stats if x[2] == x[0]),
                  ShareExactlyD=f"{100 * med(exact_d):.1f}\\%", ShareExactlyDMax=f"{100 * max(exact_d):.0f}\\%",
                  MaxDegRatioMed=f"{med([x[3] / (2 * (x[0] - 1)) for x in deg_stats]):.2f}",
                  MaxDegRatioMax=f"{max(x[3] / (2 * (x[0] - 1)) for x in deg_stats):.2f}",
                  StableIso=sum(stable["isostatic"]), StableMore=sum(stable["more"]),
                  StableFourteenIso=sum(stable14["isostatic"]), StableFourteenMore=sum(stable14["more"]),
                  NRefDiscMore=sum(over_refine), NRefDiscRest=refine_disc - sum(over_refine),
                  NRestConv=refine_tot - len(over_refine))
    for key, name in (("N/d<=8", "LowNd"), ("N/d>32", "HighNd")):
        macros[f"NAll{name}"] = sum(cls_nd[key].values())
        macros[f"NIso{name}"] = cls_nd[key]["isostatic"]
        macros[f"NMore{name}"] = cls_nd[key]["more"]
    for key, name in (("3-10", "A"), ("11-17", "B"), ("18-24", "C"), ("25-32", "D")):
        macros[f"RatShare{name}"] = f"{100 * med(rat_share[key]):.1f}\\%"
    P(f"  isostatic codes: minimum degree equals d in {sum(1 for x in deg_stats if x[2] == x[0])} of {len(deg_stats)}; "
      f"share of points with exactly d contacts: median {med(exact_d):.3f}, range {min(exact_d):.3f}-{max(exact_d):.3f};"
      f" largest degree / (2(d-1)): median {med([x[3] / (2 * (x[0] - 1)) for x in deg_stats]):.2f}, max "
      f"{max(x[3] / (2 * (x[0] - 1)) for x in deg_stats):.2f}")
    for kind in ("isostatic", "more", "fewer"):
        P(f"  {kind}: same number of pairs within 1e-12 and within 1e-8 of the largest cosine in {sum(stable[kind])} of "
          f"{len(stable[kind])}; within 1e-14 and within 1e-9 in {sum(stable14[kind])}")
    P(f"  codes with a pair of antipodal points (cosine < -1 + 1e-9): {antipodal}")
    P(f"  colour refinement of the contact graph (non-rattlers) is discrete (every vertex its own class) in "
      f"{refine_disc} of {refine_tot} codes; among the codes with more contacts than the isostatic count: "
      f"{sum(over_refine)} of {len(over_refine)}")
    macros["NRefineDiscrete"] = refine_disc
    macros["NCertDiscrete"] = sum(cert_disc)
    macros["NCertConvGraph"] = len(cert_disc)
    macros["NAntipodal"] = len(antipodal)
    P(f"  certified converged codes whose colour refinement is discrete: {sum(cert_disc)} of {len(cert_disc)}")
    macros["NRefineAll"] = refine_tot
    if stress_rows:  # quantiles 0, 1%, 10%, 50%, 90%, 99%, 100% of the certified stress / mean stress, per code
        qs = list(zip(*[x[1] for x in stress_rows]))
        names = ("min", "1%", "10%", "median", "90%", "99%", "max")
        P(f"  stresses of the {len(stress_rows)} certified isostatic converged codes, relative to the mean stress 1/|E|; "
          "for each quantile, its median over the codes (and range): " + "; ".join(
              f"{nm} {med(q):.3g} ({min(q):.2g}-{max(q):.3g})" for nm, q in zip(names, qs)))
        macros.update(StressMedMedian=f"{med(qs[3]):.2f}", StressMedTen=f"{med(qs[2]):.2f}",
                      StressMedNinety=f"{med(qs[4]):.2f}", StressMedMax=f"{med(qs[6]):.1f}",
                      StressMaxMax=f"{max(qs[6]):.0f}")

    with open(os.path.join(out, "cert_numbers.tex"), "w") as f:
        f.write("% written by certify/analysis.py; do not edit\n")
        for k, v in sorted(macros.items()):
            v = f"{v:,}".replace(",", "{,}") if isinstance(v, int) else v
            f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
    P(f"\nwrote {out}/cert_numbers.tex, tab_cert.tex, tab_iso.tex")


if __name__ == "__main__":
    main()
