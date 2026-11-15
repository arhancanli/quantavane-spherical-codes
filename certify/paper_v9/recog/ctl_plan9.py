"""Census version 9: the plan of the controls of strict9.py (--class strict_hyper, the strict certificates of version 8
by recognition). Each entry gets its number of unknowns (recog9.census_tol); those with at most --max-unknowns (2,300)
are run, in two disjoint shards (alternate entries in the order of the unknowns) so that two jobs can share them; the
larger ones are left out for memory: strict9.py's memory grows with the square of the unknowns, and the two shards
of at most 2,300 unknowns already reached peak footprints of 2.6 and 2.8 GB (macOS /usr/bin/time -l).

  nice -n 10 python3 -I ctl_plan9.py [--max-unknowns 2300]
writes ctl_plan.json, ctl_shardA.txt and ctl_shardB.txt next to this file.
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import recog9 as R9


def main():
    args = sys.argv[1:]
    cap = int(args[args.index("--max-unknowns") + 1]) if "--max-unknowns" in args else 2300
    rows = sorted((R9.census_tol(d, n)[1], d, n) for d, n in R9.open_entries("strict_hyper"))
    run = [(U, d, n) for U, d, n in rows if U <= cap]
    shards = {"A": [f"{d}_{n}" for i, (U, d, n) in enumerate(run) if i % 2 == 0],
              "B": [f"{d}_{n}" for i, (U, d, n) in enumerate(run) if i % 2 == 1]}
    for s, L in shards.items():
        open(os.path.join(HERE, f"ctl_shard{s}.txt"), "w").write(",".join(L))
    json.dump({"max_unknowns": cap, "entries": [{"d": d, "n": n, "unknowns": U, "run": U <= cap} for U, d, n in rows]},
              open(os.path.join(HERE, "ctl_plan.json"), "w"), indent=0)
    print(len(rows), "certificates by recognition;", len(run), "with at most", cap, "unknowns are run")


if __name__ == "__main__":
    main()
