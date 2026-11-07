\\ SPDX-License-Identifier: MIT
\\ Copyright (c) 2026 Arhan Canli (Quantavane)
\\
\\ field.gp -- the field L = Q(sqrt145, s1, s2) of the note in reduced form (Section 3, "The values"); descriptive, exact.
\\ Run through run_gp.py (DATA = data/), or gp -q field.gp after DATA = "path/to/data/".
\\ From data/values.gp: theta = s1 + s2 has minimal polynomial f (a resultant over Q(sqrt145)); polredabs(f) gives a
\\ reduced defining polynomial g of L and the image of theta; r = sqrt145, s1, s2 are written as polynomials in a root y
\\ of g, chosen at the real embedding of the note (r > 0, s1 < 0, s2 < 0, where theta is the smallest real root of f);
\\ every relation is checked in Q[y]/(g). Also nfdisc(g), and the identity y^2 + y^-2 = (15 + r)/2.

default(realprecision, 60);
say(k, v) = print(k, " = ", v);
need(b, msg) = if (!b, error(msg));
read(Str(DATA, "values.gp"));
Dr1 = D1[1] + D1[2]*r; Dr2 = D2[1] + D2[2]*r;
f = polresultant((x^2 - Dr1 - Dr2)^2 - 4*Dr1*Dr2, r^2 - 145, r);
f = f / pollead(f);
need(polisirreducible(f) && poldegree(f) == 8, "f is not irreducible of degree 8");
R = polredabs(f, 1);
g = subst(R[1], x, y);
th = Mod(subst(lift(R[2]), x, y), g);
need(subst(f, x, th) == 0, "polredabs: wrong image of theta");
K = nfinit(g);
\\ the real embedding: theta is the smallest real root of f
thnum = vecmin(real(select(z -> abs(imag(z)) < 1e-40, polroots(f))));
yr = real(select(z -> abs(imag(z)) < 1e-40, polroots(g)));
yi = select(z -> abs(subst(lift(th), y, z) - thnum) < 1e-40, yr);
need(#yi == 1, "the real root y0 is not unique");
y0 = yi[1];
ev(e) = subst(lift(e), y, y0);
pick(pol, sgn) = {
  my(rt = nfroots(K, pol), s = select(e -> sign(ev(Mod(e, g))) == sgn, rt));
  need(#s == 1, "root not determined by its sign"); Mod(s[1], g);
}
rK = pick(x^2 - 145, 1);
s1K = pick(x^2 - subst(Dr1, r, lift(rK)), -1);
s2K = pick(x^2 - subst(Dr2, r, lift(rK)), -1);
need(s1K + s2K == th, "theta != s1 + s2");
need(Mod(y, g)^2 + Mod(y, g)^-2 == (15 + rK)/2, "y^2 + y^-2 != (15 + sqrt145)/2");
say("field_polredabs", g);
say("field_disc_factor", factor(nfdisc(g)));
say("field_disc", nfdisc(g));
say("field_real_roots", #yr);
say("field_y0", y0);
say("field_sqrt145_in_y", lift(rK));
say("field_s1_in_y", lift(s1K));
say("field_s2_in_y", lift(s2K));
say("field_y2_plus_yinv2", "(15 + sqrt145)/2");
say("STATUS", "FIELD CHECKS PASSED");
