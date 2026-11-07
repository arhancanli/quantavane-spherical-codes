"""External validation: the table's exact entries against the certified enclosures.

  python3 -I certify/census/validate_poly.py RESULTS_DIR OUT.tsv

For every code of the archived table that lists the minimal polynomial p of its cosine and that the census certified
(PASS, LOCALMIN or NONOPT), the enclosure [lo, hi] of the cosine c* of the exact configuration (the decimal interval
of 30 digits written by cq.py, which contains the certified enclosure) is tested exactly: p(lo) and p(hi) are
evaluated in rational arithmetic. p is irreducible, so its roots are simple, and a sign change (or a zero at an end)
proves that a root of p lies in [lo, hi]; no sign change means that [lo, hi] contains no root or at least two.
"""
import json, os, re, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snapshot


def parse_poly(s):
    """'3*x^3-9*x^2-3*x+1' -> {3: 3, 2: -9, 1: -3, 0: 1}"""
    s = s.replace(" ", "")
    if s[0] not in "+-":
        s = "+" + s
    coef = {}
    for sign, c, x, e in re.findall(r"([+-])(\d*)\*?(x?)(?:\^(\d+))?", s):
        if not c and not x:
            continue
        k = int(e) if e else (1 if x else 0)
        v = int(c) if c else 1
        coef[k] = coef.get(k, 0) + (v if sign == "+" else -v)
    return coef


def peval(coef, q):
    return sum(Fraction(v) * q ** k for k, v in coef.items())


def main():
    res, out = sys.argv[1], sys.argv[2]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    rows = []
    for (d, n), e in sorted(idx.items()):
        if not e["poly"] or not e["listed"] or d > 32 or n > 1024:
            continue
        js = os.path.join(res, f"d{d}_n{n}.json")
        status, verdict, lo, hi, kind = "not run", "", "", "", ""
        if os.path.exists(js):
            j = json.load(open(js))
            status = j.get("status", "?")
            att = [a for a in j.get("attempts", []) if a.get("status") in ("PASS", "LOCALMIN", "NONOPT")]
            if att:
                a = att[0]
                kind = a.get("kind", "")
                lo, hi = a["c_interval_30"]
                ql, qh = Fraction(int(lo), 10 ** 30), Fraction(int(hi), 10 ** 30)
                coef = parse_poly(e["poly"])
                pl, ph = peval(coef, ql), peval(coef, qh)
                verdict = "root inside" if pl == 0 or ph == 0 or (pl > 0) != (ph > 0) else "NO SIGN CHANGE"
        rows.append((d, n, e["poly"], e["cos"], status, kind, lo, hi, verdict))
    with open(out, "w") as f:
        f.write("d\tn\tpoly\ttable_cos\tstatus\tkind\tlo_e30\thi_e30\tverdict\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    tested = [r for r in rows if r[-1]]
    print(len(rows), "exact entries;", len(tested), "certified and tested;",
          sum(r[-1] == "root inside" for r in tested), "root inside;",
          sum(r[-1] == "NO SIGN CHANGE" for r in tested), "without sign change")


if __name__ == "__main__":
    main()
