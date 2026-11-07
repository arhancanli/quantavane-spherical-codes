"""The sentences of the paper (version 4) whose content is a list of entries, from the census outputs.

  python3 -I certify/census/make_text.py CENSUSDIR PAPERDIR      # writes PAPERDIR/v4_text.tex
"""
import collections, glob, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snapshot


def sci(x, digits=1):
    if x is None:
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


def ent(k):
    return f"$({k[0]},{k[1]})$"


def tau(t):
    return f"$\\tau=10^{{{int(round(math.log10(t)))}}}$"


def load(cdir, sub):
    out = {}
    for f in glob.glob(os.path.join(cdir, sub, "*.json")):
        j = json.load(open(f))
        if "d" in j:
            out[(j["d"], j["n"])] = j
    return out


def main():
    cdir, pdir = sys.argv[1], sys.argv[2]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    summ = json.load(open(os.path.join(cdir, "summary_v4.json")))["status"]
    status = {tuple(int(x) for x in k[1:].split("_n")): v for k, v in summ.items()}
    res, strict, exs = load(cdir, "results_v4"), load(cdir, "strict_v4"), load(cdir, "exact_v4/strict")
    T = {}
    # provably optimal entries that are not certified strict
    words = {"local": "hyperstatic, strictness open", "flex": "contact set excluded", "failed": "certificate failed",
             "too_large": "above the size limit", "nonopt_iso": "$Q^*$ not a local minimum"}
    opt = sorted(k for k in status if idx[k]["optimal"])
    rest = []
    for k in opt:
        st = status[k]
        if st in ("strict_iso", "strict_hyper"):
            continue
        w = words.get(st, st)
        if st == "local" and k in strict:
            w += " (" + strict[k]["status"].replace("_", " ") + ")"
        rest.append(f"{ent(k)} ({w})")
    exs_ok = {k for k, j in exs.items() if j.get("status", "").startswith("STRICT")}
    exact_opt = [k for k in opt if status[k] not in ("strict_iso", "strict_hyper") and k in exs_ok]
    T["ExactOnlyOptText"] = (("Exact arithmetic alone certifies " + str(len(exact_opt)) + " more entries marked provably "
                              "optimal, " + ", ".join(ent(k) for k in exact_opt) + ", whose systems exceed the size "
                              "limit (Section~\\ref{sec:hyperres}).") if exact_opt else "")
    T["OptNotStrictText"] = (("the others are " + "; ".join(rest) + ".") + (
        " Of these, " + ", ".join(ent(k) for k in exact_opt) + " are certified strict by exact arithmetic alone "
        "(Section~\\ref{sec:hyperres})." if exact_opt else "")) if rest else \
        "every provably optimal entry in the census is certified strict."
    sh = [k for k in opt if status[k] == "strict_hyper"]
    T["OptStrictHyperText"] = ("the provably optimal entries " + ", ".join(ent(k) for k in sh) + ".") if sh else \
        "no provably optimal entry."
    # isostatic NONOPT, grouped by entry
    groups = collections.OrderedDict()
    for k in sorted(res):
        for a in res[k]["attempts"]:
            if a.get("status") == "NONOPT":
                groups.setdefault(k, []).append(a)
    cls = {"flex": "contact set excluded", "strict_iso": "strict", "local": "strictness open",
           "strict_hyper": "strict", "failed": "certificate failed", "too_large": "above the size limit"}
    rows = []
    for k, atts in groups.items():
        taus = ", ".join(tau(a["tol"]).strip("$").replace("\\tau=", "") for a in atts)
        a = atts[0]
        q = a.get("stress_quantiles_over_mean") or [None]
        far = "; this concerns a configuration rather far from the code" if a["dist_ls"] > 1e-6 else ""
        rows.append(f"{ent(k)} at $\\tau={taus}$ ($Q^*$ within {sci(a['dist_ls'])} of the stored code{far}; smallest "
                    f"relative stress {sci(q[0])})")
    T["NonoptText"] = ("Proposition~\\ref{prop:nonopt} applies to " + " and ".join(rows) + "; " +
                       ("both entries carry" if len(rows) > 1 else "the entry carries") +
                       " the flexibility statement at another tolerance and are counted there.") if rows else \
        "Proposition~\\ref{prop:nonopt} applies to no entry."
    # exact arithmetic (exact_strict.py)
    tried = {k: j for k, j in exs.items() if j.get("status") != "coordinates not recognised as exact"}
    ok = [k for k, j in tried.items() if j["status"].startswith("STRICT")]
    reasons = collections.defaultdict(list)
    phr = {"the uniform stress is not in equilibrium": "the uniform stress is not an equilibrium stress (the exact search "
                                                      "for another one covers only rational coordinates)",
           "rigidity matrix above the memory limit": "the rigidity matrix does not fit the memory limit",
           "no exact positive stress found": "no exact positive stress was found",
           "bar-rigidity not shown": "the rank of the rigidity matrix modulo $q$ is too small"}
    for k, j in sorted(tried.items()):
        if not j["status"].startswith("STRICT"):
            r = j["status"]
            r = "bar-rigidity not shown" if r.startswith("bar-rigidity") else r
            reasons[phr.get(r, r)].append(k)
    cpu = max([j.get("cpu_seconds", 0) for j in tried.values()] or [0])
    s = (f"The exact check was tried on {len(exs)} entries (every hyperstatic certificate, and the entries marked "
         f"provably optimal that the isostatic certificate does not cover); for {len(tried)} of them the coordinates "
         f"or the inner products were recognised as exact, and {len(ok)} of these are confirmed, each in at most "
         f"{cpu:.0f} seconds of CPU time")
    if reasons:
        s += ". For the other " + str(len(tried) - len(ok)) + " it gives no answer: " + "; ".join(
            f"for {len(v)}" + (" (" + ", ".join(ent(k) for k in v) + ")" if len(v) <= 3 else "") + f", {r}"
            for r, v in sorted(reasons.items(), key=lambda x: -len(x[1])))
    summ_status = status
    above = [k for k in ok if summ_status.get(k) in ("too_large", "failed", "flex")]
    if above:
        s += (". Of the confirmed entries, " + ", ".join(ent(k) for k in above) + " have no certificate from the census "
              "rule (their systems are above its size limit), so for them this is the only certificate; "
              "Table~\\ref{tab:census} counts them under the census rule")
    T["ExactStrictText"] = s + "."
    # validation
    vt = os.path.join(cdir, "validate_v4.tsv")
    if os.path.exists(vt):
        rows = [l.rstrip("\n").split("\t") for l in open(vt)][1:]
        ns = [r for r in rows if r[-1] == "NO SIGN CHANGE"]
        T["ExactValidationRest"] = "There is no case without a sign change." if not ns else \
            ("The interval shows no sign change for " + ", ".join(f"$({r[0]},{r[1]})$" for r in ns) + ".")
    # exact_big
    big = []
    for f in sorted(glob.glob(os.path.join(cdir, "exact_v4", "big", "*.json"))):
        big.append(json.load(open(f)))
    if big:
        parts = []
        lab = {"iso": "isostatic certificates (existence, stress and the gaps of the non-rattlers)",
               "hyper": "hyperstatic certificates (existence for $E'$, the combined stress, $\\chi$, $\\kappa$, $\\rho$, "
                        "the displacement and the gaps)",
               "flex": "flexibility certificates (the row sums of $|I-LC|$ and the radius $h$)",
               "strict": "strictness certificates with rational inner products (the Gram matrix, a strictly positive "
                         "stress and bar-rigidity, with no interval arithmetic)"}
        for mode in ("iso", "hyper", "flex", "strict"):
            js = [j for j in big if j.get("mode") == mode]
            if not js:
                continue
            good = [j for j in js if j["status"].startswith(("PASS", "LOCALMIN confirmed", "FLEX confirmed",
                                                              "STRICT confirmed")) and "FAIL" not in j["status"]]
            mx = max(js, key=lambda j: j.get("unknowns", j.get("rows", 0)) or 0)
            p = f"{len(good)} of {len(js)} {lab[mode]}"
            if len(good) < len(js):
                p += " (not confirmed: " + ", ".join(f"$({j['d']},{j['n']})$, {j['status']}" for j in js if j not in good) + ")"
            parts.append(p)
        secs = [j.get("seconds", 0) for j in big]
        agree, diffs = 0, []
        for j in big:
            if j.get("mode") == "iso" and "c_lo_e30" in j:
                a = [x for x in res[(j["d"], j["n"])]["attempts"] if x.get("status") == "PASS"][0]
                lo, hi = (int(x) for x in a["c_interval_30"])
                dl, dh = abs(j["c_lo_e30"] - lo), abs(j["c_hi_e30"] - hi)
                diffs.append(max(dl, dh))
                agree += (dl == 0 and dh == 0)
        cmp_ = (f" For the {len(diffs)} isostatic re-checks the exact 30-digit enclosures of $c^*$ coincide with the "
                f"certifier's in {agree} cases and differ by at most {max(diffs)} units in the 30th digit in the others.") \
            if diffs and agree < len(diffs) else (f" For the {len(diffs)} isostatic re-checks the exact 30-digit "
                                                     "enclosures of $c^*$ coincide with the certifier's." if diffs else "")
        T["ExactCheckText"] = ("The separately written program \\nolinkurl{certify/census/exact_big.py} (it imports none of "
                               "our certifiers, and every quantity that enters a decision is an exact integer or rational) "
                               "re-checked a sample drawn reproducibly from every certificate type (the largest systems and "
                               "random ones; \\nolinkurl{certify/census/exact_sample.py}). It confirmed " + "; ".join(parts) +
                               f". The runs took from {min(secs):.0f} to {max(secs):.0f} seconds." + cmp_)
    T.setdefault("ExactCheckText", "(The exact re-checks of exact\\_big.py have not been run.)")
    T.setdefault("ExactValidationRest", "")
    with open(os.path.join(pdir, "v4_text.tex"), "w") as f:
        f.write("% written by certify/census/make_text.py (version 4) from the census outputs; do not edit\n")
        for k, v in T.items():
            f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
    print(json.dumps(T, indent=1)[:3000])


if __name__ == "__main__":
    main()
