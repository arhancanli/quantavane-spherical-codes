"""CENSUS.md for the paper (version 4): a plain summary of the census, from its outputs.

  python3 -I certify/census/census_md.py CENSUSDIR > CENSUS.md
"""
import collections, glob, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snapshot
import run_census

LAB = [("strict_iso", "exact strict local minimum certified nearby (isostatic)"),
       ("strict_hyper", "exact strict local minimum certified nearby (hyperstatic, exact Gram matrix)"),
       ("local", "exact local minimum certified nearby, strictness open (hyperstatic)"),
       ("flex", "diagnostic: contact set at tau excluded (no nearby local minimum has all its contacts in it)"),
       ("failed", "nothing certified: every attempted certificate failed"),
       ("too_large", "nothing certified: every system above the size limit (3,100 unknowns)")]


def f(k):
    return f"{k:,}"


def main():
    cdir = sys.argv[1]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    S = json.load(open(os.path.join(cdir, "summary_v4.json")))
    M = S["macros"]
    status = {tuple(int(x) for x in k[1:].split("_n")): v for k, v in S["status"].items()}
    keys = sorted(status)
    snap = json.load(open(os.path.join(cdir, "snapshot_files.json")))
    out = ["# Census of Cohn's archived table of spherical codes (paper version 4)", "",
           "Data set: DSpace@MIT handle 1721.1/153543, data as of 31 July 2026 (downloaded 6 November 2026; DSpace "
           "lists the item as last modified 4 August 2026). Every file's MD5 agrees with DSpace's.", "",
           "| file | bytes | SHA-256 |", "|---|---:|---|"]
    for name in ["table.txt"] + [f"data-{k}.txt" for k in range(1, 13)] + ["data.txt"]:
        out.append(f"| {name} | {f(snap[name]['bytes'])} | `{snap[name]['sha256']}` |")
    nfl = sum(idx[k]["flagged"] for k in keys)
    nopt = sum(idx[k]["optimal"] for k in keys)
    out += ["", f"Census: the {f(len(keys))} entries that `table.txt` lists with 3 <= d <= 32 and N <= 1024 ({f(nfl)} "
            f"marked \"needs more optimization\", {nopt} marked \"provably optimal\"). One rule for every entry "
            "(paper, Section 4.2): isostatic certificates at every tolerance where the proposed graph is isostatic, then "
            "hyperstatic ones, then (if none succeeded) the flexibility certificate at the largest tolerance where the graph "
            "is hypostatic or empty; the strictness stage for every hyperstatic certificate. An attempt is made only if the "
            "system it solves has at most 3,100 unknowns (rows).", "",
            "What each class means: the objects of the theorems are exact configurations near the stored codes, never the "
            "stored codes; every distance between them is in the outputs. No entry is proved to have no local minimum "
            "nearby.", "", "## Counts by class", "", "| class | all | marked | unmarked |", "|---|---:|---:|---:|"]
    for c, lab in LAB:
        ks = [k for k in keys if status[k] == c]
        fl = sum(idx[k]["flagged"] for k in ks)
        out.append(f"| {lab} | {f(len(ks))} | {f(fl)} | {f(len(ks) - fl)} |")
    out.append(f"| total | {f(len(keys))} | {f(nfl)} | {f(len(keys) - nfl)} |")
    out += ["", "## Counts by dimension", "",
            "| d | entries | strict iso. | strict hyper. | open | excluded | failed | too large |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for d in range(3, 33):
        ks = [k for k in keys if k[0] == d]
        c = collections.Counter(status[k] for k in ks)
        out.append(f"| {d} | {f(len(ks))} | " + " | ".join(f(c.get(x, 0)) for x, _ in LAB) + " |")
    # provably optimal
    out += ["", "## The entries marked provably optimal", "", "| entry | class |", "|---|---|"]
    for k in keys:
        if idx[k]["optimal"]:
            out.append(f"| ({k[0]},{k[1]}) | {dict(LAB)[status[k]]} |")
    # strictness stage
    strict = {}
    for fn in glob.glob(os.path.join(cdir, "strict_v4", "*.json")):
        j = json.load(open(fn))
        strict[(j["d"], j["n"])] = j
    sc = collections.Counter(j["status"] for j in strict.values())
    out += ["", "## Strictness stage (hyperstatic certificates)", "", "| outcome | entries |", "|---|---:|"]
    for st, v in sc.most_common():
        out.append(f"| {st} | {f(v)} |")
    # exact arithmetic
    exs = {}
    for fn in glob.glob(os.path.join(cdir, "exact_v4", "strict", "*.json")):
        j = json.load(open(fn))
        exs[(j["d"], j["n"])] = j
    ok = sorted(k for k, j in exs.items() if j.get("status", "").startswith("STRICT"))
    out += ["", "## Exact arithmetic alone (exact_strict.py)", "",
            f"Confirmed with no floating-point arithmetic in any decision: " +
            ", ".join(f"({d},{n})" for d, n in ok) + "."]
    big = [json.load(open(fn)) for fn in sorted(glob.glob(os.path.join(cdir, "exact_v4", "big", "*.json")))]
    if big:
        out += ["", "## Independent exact re-checks (exact_big.py, sample of exact_sample.py)", "",
                "| mode | entry | unknowns | status | seconds |", "|---|---|---:|---|---:|"]
        for j in big:
            out.append(f"| {j['mode']} | ({j['d']},{j['n']}) | {j.get('unknowns', j.get('rows', ''))} | {j['status']} | "
                       f"{j.get('seconds', '')} |")
    # tolerance
    recs = {(r["d"], r["n"]): r for r in map(json.loads, open(os.path.join(cdir, "classify.jsonl")))
            if (r["d"], r["n"]) in status}
    out += ["", "## Classes of the proposed graphs by tolerance", "", "| tau | isostatic | hyperstatic | hypostatic | empty |",
            "|---|---:|---:|---:|---:|"]
    for t in [1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]:
        cl = collections.Counter(run_census.klass({x["tol"]: x for x in recs[k]["structs"]}[t]) for k in keys)
        out.append(f"| {t:g} | " + " | ".join(f(cl.get(x, 0)) for x in ("iso", "hyper", "hypo", "empty")) + " |")
    # versions
    shas = collections.Counter()
    for fn in glob.glob(os.path.join(cdir, "results_v4", "*.json")):
        shas[json.dumps(json.load(open(fn)).get("code_sha256"), sort_keys=True)] += 1
    out += ["", "## Programs", "", "Every output file records the SHA-256 of the programs that produced it. The census "
            "outputs (results_v4/) were produced by one version of the certifier (cq.py, sphcert.py); the driver "
            "run_census.py has two versions (the second runs hyperstatic attempts under a memory watchdog):", ""]
    for js, cnt in shas.most_common():
        out.append(f"- {cnt:,} entries: " + ", ".join(f"`{k}` `{v}`" for k, v in json.loads(js).items()))
    st = collections.Counter()
    for fn in glob.glob(os.path.join(cdir, "strict_v4", "*.json")):
        st[json.dumps(json.load(open(fn)).get("code_sha256"), sort_keys=True)] += 1
    for js, cnt in st.most_common():
        out.append(f"- strictness stage, {cnt:,} entries: " + ", ".join(f"`{k}` `{v}`" for k, v in json.loads(js).items()))
    out += ["", f"Macros of the paper: certify/census/summary_v4.json ({len(M)} values)."]
    print("\n".join(out))


if __name__ == "__main__":
    main()
