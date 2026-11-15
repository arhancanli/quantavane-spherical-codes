"""Census version 9, item A5: the index of the exact descriptions (exactdata/INDEX.tsv).

  python3 -I index9.py

One line per exactly described strict code: the route (recognition in version 9 or 8, or the exact core), the field
(the defining polynomial of a primitive element, or Q(sqrt D)), the minimal polynomial of the cosine (integer
coefficients, constant first; computed exactly with FLINT as the irreducible factor of the characteristic polynomial of
multiplication by the cosine that vanishes at it), the cosine to 40 digits, the numbers of distinct inner products,
non-rattlers and contacts, the order of the symmetry group found by the numerical survey of version 6 where it was run
(paper_v8/phase2/survey_*.jsonl; a lower bound found by an untrusted search, empty where not surveyed), and the file.
"""
import glob, gzip, json, os, sys
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
import flint


def survey():
    out = {}
    for f in glob.glob(os.path.join(CERT, "paper_v8", "phase2", "survey_*.jsonl")):
        for l in open(f):
            if l.strip():
                j = json.loads(l)
                out[(j["d"], j["n"])] = j.get("group_order")
    return out


def minpoly_of(fp, cel):
    """The minimal polynomial over Q of the element cel (coefficients as strings, power basis) of Q[x]/(fp)."""
    p = flint.fmpq_poly([int(c) for c in fp])
    D = p.degree()
    c = flint.fmpq_poly([flint.fmpq(Fraction(s).numerator, Fraction(s).denominator) for s in cel]) % p
    cols = []
    for k in range(D):
        e = (c * flint.fmpq_poly([0] * k + [1])) % p
        co = e.coeffs() + [flint.fmpq(0)] * (D - len(e.coeffs()))
        cols.append(co)
    M = flint.fmpq_mat(D, D, [cols[j][i] for i in range(D) for j in range(D)])
    ch = M.charpoly()
    # the irreducible factor of the characteristic polynomial that vanishes at c in Q[x]/(p)
    for fac, _ in ch.factor()[1]:
        acc = flint.fmpq_poly([0])
        for co in reversed(fac.coeffs()):
            acc = (acc * c + flint.fmpq_poly([co])) % p
        if acc == 0:
            den = 1
            for q in fac.coeffs():
                den = den * int(q.q) // __import__("math").gcd(den, int(q.q))
            ints = [int(q * den) for q in fac.coeffs()]
            g = 0
            for x in ints:
                g = __import__("math").gcd(g, x)
            ints = [x // g for x in ints]
            return ints if ints[-1] > 0 else [-x for x in ints]
    raise ValueError("no factor vanishes")


def main():
    sv = survey()
    lines = ["d\tn\troute\tfield\tfield_degree\tcosine_minimal_polynomial\tcosine\tdistinct_inner_products\tnonrattlers\t"
             "contacts\tsymmetry_group_order_survey\tfile"]
    for sub, route in (("open", "recognition, version 9"), ("v8recog", "recognition, version 8"), ("core", "exact core")):
        for f in sorted(glob.glob(os.path.join(HERE, sub, "*.json.gz"))):
            j = json.load(gzip.open(f, "rt"))
            d, n = j["d"], j["n"]
            if sub == "core":
                D = j["field_D"]
                field = "Q" if D == 1 else f"Q(sqrt({D}))"
                deg = 1 if D == 1 else 2
                a, b = (Fraction(x) for x in j["cosine"])
                if b == 0:
                    mp = [-a.numerator, a.denominator]
                else:                                     # (x - a)^2 - D b^2
                    q = [a * a - D * b * b, -2 * a, Fraction(1)]
                    den = 1
                    for x in q:
                        den = den * x.denominator // __import__("math").gcd(den, x.denominator)
                    mp = [int(x * den) for x in q]
                cval = str(float(a) + float(b) * D ** 0.5)
            else:
                fp = j["field_poly"]
                field = "Q[x]/(" + " ".join(str(c) for c in fp) + ")"
                deg = len(fp) - 1
                mp = minpoly_of(fp, j["c_elem"])
                cval = j["c_star"]
            lines.append("\t".join(map(str, [d, n, route, field, deg, " ".join(map(str, mp)), cval, len(j["elements"]),
                                             len(j["nonrattlers"]), len(j["contacts"]), sv.get((d, n), ""),
                                             sub + "/" + os.path.basename(f)])))
    open(os.path.join(HERE, "INDEX.tsv"), "w").write(
        "# The exact descriptions of the strict codes (census version 9; exactdata/index9.py). Polynomials: integer\n"
        "# coefficients, constant first. Field: Q[x]/(p) with p's coefficients, constant first, or Q(sqrt(D)).\n" +
        "\n".join(lines) + "\n")
    print(len(lines) - 1, "codes indexed")


if __name__ == "__main__":
    main()
