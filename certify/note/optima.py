# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""The entries of the archived table in dimensions 3 and 4 marked "Provably optimal": are they infinitesimally jammed?
(UNTRUSTED floating point: a screening, not a certificate; for the remark on low dimensions in the note.)

  python3 -I certify/note/optima.py SNAPDIR OUT.json

A code that is optimal and unique up to isometry is jammed (along a deformation the largest cosine stays optimal, so
every code passed through is congruent to the original, and the Gram matrix, which can only change among finitely many
relabellings of one matrix, is constant). So a uniquely optimal code that is not infinitesimally jammed would answer
the question of Cohn, Jiao, Kumar and Torquato. For each entry marked provably optimal with d = 3 or 4, this program
reads the stored decimals and reports, with family_table.rigidity: the contacts (the inner products within 1e-9 of the
largest), rattlers, and for the core (the points in some contact) the infinitesimal flexes modulo rotations of the
contact framework (dense SVD) and a strictly positive equilibrium stress (linear programming). An entry is
infinitesimally jammed (numerically) if it has no rattlers, no such flex and a strictly positive stress (Lemma 2.2 of
the note).
"""
import json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [p for p in os.environ.get("QV_LIBS", "").split(os.pathsep) if p]   # optional library folders
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "census"))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.seterr(all="ignore")
from snapshot import load_index, read_points
from family_table import classes, rigidity


def main():
    snapdir, out = sys.argv[1], sys.argv[2]
    t0 = time.time()
    res = []
    for e in load_index():
        if e["d"] not in (3, 4) or not e["optimal"]:
            continue
        X = np.array([[float(t) for t in row] for row in read_points(snapdir, e)])
        X /= np.linalg.norm(X, axis=1)[:, None]
        C, K, spread, G = classes(X)
        rep = rigidity(X, C, 60)
        rep.update(d=e["d"], n=e["n"], notes=e["notes"][:80])
        rep["infinitesimally_jammed_numerically"] = bool(rep["rattlers"] == 0 and rep["flexes_mod_rotations"] == 0
                                                          and rep["strictly_positive_stress"])
        res.append(rep)
        print(json.dumps({k: rep[k] for k in ("d", "n", "contacts", "rattlers", "flexes_mod_rotations",
                                              "strictly_positive_stress", "infinitesimally_jammed_numerically")}),
              flush=True)
    summ = {"program": "optima.py", "entries": len(res),
            "without_rattlers": sum(1 for r in res if r["rattlers"] == 0),
            "infinitesimally_jammed_numerically": sum(1 for r in res if r["infinitesimally_jammed_numerically"]),
            "with_rattlers": [(r["d"], r["n"]) for r in res if r["rattlers"]],
            "seconds": round(time.time() - t0, 1), "entries_detail": res}
    json.dump(summ, open(out, "w"), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k != "entries_detail"}))


if __name__ == "__main__":
    main()
