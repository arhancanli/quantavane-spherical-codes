"""Census version 9: an external spot check of the exact descriptions with PARI/GP (cypari, PARI 2.15), an established
computer algebra system that shares no code with FLINT or with our programs.

  nice -n 10 python3 -I pari_check9.py OUT.json FILE.json.gz [FILE.json.gz ...]

For each exact description written by strict9.py --data (field polynomial p, the elements of the distinct inner
products in Q[x]/(p), the index matrix of the Gram matrix, the basis B, the contacts and the element of the cosine):
  - p is irreducible over Q (polisirreducible);
  - every contact entry of G equals the element of the cosine, and the diagonal is 1;
  - G = G_VB G_BB^-1 G_BV over Q[x]/(p) (PARI's t_POLMOD arithmetic; G_BB^-1 by PARI's matrix inverse), so G has rank d.
These are the exact algebraic parts of the certificate (Proposition 3.15); the signs at the real embedding and the
stress test are not repeated here.
"""
import gzip, hashlib, json, os, sys, time
sys.path.insert(0, os.environ.get("QV_PARILIB", os.path.expanduser("~/Library/Caches/quantavane-census/parilib")))
from cypari import pari
pari.allocatemem(512 * 2 ** 20, 512 * 2 ** 20)          # PARI stack of 512 MB (the default 8 MB is too small)

SHA = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()


def polmod(coeffs, P):
    terms = " + ".join(f"({c})*x^{i}" for i, c in enumerate(coeffs)) or "0"
    return f"Mod({terms}, {P})"


def check(path):
    t0 = time.time()
    j = json.load(gzip.open(path, "rt"))
    P = " + ".join(f"({c})*x^{i}" for i, c in enumerate(j["field_poly"]))
    irred = bool(int(pari(f"polisirreducible({P})")))
    els = [polmod(e, P) for e in j["elements"]]
    pari("E9 = [" + ", ".join(els) + "]")
    cid = j["index_matrix"]
    n = len(cid)
    pari("I9 = [" + "; ".join(", ".join(str(x + 1) for x in row) for row in cid) + "]")
    pari(f"G9 = matrix({n}, {n}, i, j, E9[I9[i, j]])")
    diag_one = bool(int(pari(f"prod(i = 1, {n}, G9[i, i] == 1)")))
    cel = polmod(j["c_elem"], P)
    pari(f"C9 = {cel}")
    ok_c = all(bool(int(pari(f"G9[{a + 1}, {b + 1}] == C9"))) for a, b in j["contacts"])
    Bv = "[" + ", ".join(str(b + 1) for b in j["B"]) + "]"
    pari(f"GB9 = vecextract(G9, {Bv}, {Bv}); GVB9 = vecextract(G9, \"1..{n}\", {Bv})")
    rank_d = bool(int(pari("GVB9 * GB9^(-1) * GVB9~ == G9")))
    return {"d": j["d"], "n": j["n"], "degree": len(j["field_poly"]) - 1, "points": n, "irreducible": irred,
            "unit_diagonal": diag_one, "contacts_equal_cosine": ok_c, "gram_rank_d": rank_d,
            "all_ok": irred and diag_one and ok_c and rank_d, "seconds": round(time.time() - t0, 1)}


def main():
    out = sys.argv[1]
    res = []
    for f in sys.argv[2:]:
        r = check(f)
        r["file"] = os.path.basename(f)
        res.append(r)
        print(r, flush=True)
    ok = [r for r in res if r["all_ok"]]
    degs = sorted(set(r["degree"] for r in ok))
    short = f"{len(ok)} recognised codes (fields of degree {degs[0]} to {degs[-1]})" if ok else "no codes"
    long_ = (f"PARI/GP 2.15 (through cypari) confirms, for {len(ok)} of the {len(res)} exact descriptions checked "
             f"(\\nolinkurl{{ranks/pari_check9.py}}), that the field polynomial is irreducible, that every contact entry "
             f"equals the cosine and that $G=G_{{VB}}G_{{BB}}^{{-1}}G_{{BV}}$ over $\\Q[x]/(p)$, with PARI's own "
             f"polynomial and matrix arithmetic.")
    json.dump({"program_sha256": SHA, "pari": ".".join(map(str, pari.version())), "results": res,
               "summary_short": short, "summary_long": long_}, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
