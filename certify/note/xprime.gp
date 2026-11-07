\\ SPDX-License-Identifier: MIT
\\ Copyright (c) 2026 Arhan Canli (Quantavane)
\\
\\ xprime.gp -- the second code X' of the note (Section 7): X* with the two points of Q reflected in the hyperplane
\\ spanned by the orbit A. Exact, in PARI/GP over the field L of degree 8 (as structure.gp). Run through
\\ run_gp.py with DATA = the data directory of X* (data/); writes the data directory of X' (data_xprime/, next to
\\ data/), which structure.gp and check_xprime.gp then read.
\\
\\ For u, v in X*, with n_A the unit normal of span(A) and Pi_A the orthogonal projection onto span(A),
\\   <u, n_A><v, n_A> = <u, v> - <u, Pi_A v> = G[u,v] - G[u,LA] G[LA,LA]^-1 G[LA,v],
\\ LA a set of 21 points of A spanning span(A) (G[LA,LA] invertible, checked). The reflection r_A(x) = x - 2<x,n_A>n_A
\\ changes only the inner products between Q and the other points: G'[u,q] = G[u,q] - 2<u,n_A><q,n_A>. The program
\\ checks that every entry of G' is one of thirteen values t0' < t1' < t2 < ... < t12 = c, where t0', t1' replace the
\\ values t0, t1 of X* on the four pairs between P and Q, determines the class matrix of X', checks that the
\\ generators of Aut(X*) preserve it, and writes values.gp, class_matrix.gp, basis.gp, orbitals.gp (the orbits of
\\ Aut(X*) on pairs, with the classes of X') and aut_generators.gp (the same generators; the Galois relabellings of X' are found by galois_xprime.py).

read(Str(DATA, "values.gp"));
read(Str(DATA, "orbitals.gp"));
read(Str(DATA, "aut_generators.gp"));
CM = readvec(Str(DATA, "class_matrix.gp"));
B = readvec(Str(DATA, "basis.gp"))[1];
n = 116; d = 22;
default(realprecision, 60);
out(k, v) = print(k, " = ", v);
OUTDIR = Str(DATA, "../data_xprime/");

\\ ------------------------------------------------------------------ the field (as in structure.gp)
f1 = (x^2 - D1[1])^2 - 145*D1[2]^2;
rx = (x^2 - D1[1])/D1[2];
Pabs = polresultant(f1, (y - x)^2 - (D2[1] + D2[2]*rx), x); Pabs = Pabs/pollead(Pabs);
if (!polisirreducible(Pabs) || poldegree(Pabs) != 8, error("Pabs"));
TH = Mod(y, Pabs);
g1 = gcd(f1, (TH - x)^2 - (D2[1] + D2[2]*rx));
if (poldegree(g1) != 1, error("gcd"));
s1 = -polcoeff(g1, 0)/polcoeff(g1, 1); s2 = TH - s1; r = (s1^2 - D1[1])/D1[2];
if (r^2 != 145 || s1^2 != D1[1] + D1[2]*r || s2^2 != D2[1] + D2[2]*r, error("relations"));
BAS = [1, r, s1, r*s1, s2, r*s2, s1*s2, r*s1*s2];
BM = matrix(8, 8, i, j, polcoeff(lift(BAS[j]), i - 1)); BMi = BM^-1;
tow(a) = BMi * vector(8, i, polcoeff(lift(a), i - 1))~;
twr(v) = sum(j = 1, 8, v[j]*BAS[j]);
refine() = {my(m = (RI[1] + RI[2])/2, v = subst(Pabs, y, m)); if (v == 0, error("rational root"));
  if (sign(v) == sign(subst(Pabs, y, RI[1])), RI = [m, RI[2]], RI = [RI[1], m]);}
sgn(a) = {
  my(A = lift(a), m, w, M, bnd, val);
  if (type(A) != "t_POL", return(sign(A)));
  if (A == 0, return(0));
  for (it = 1, 4000,
    m = (RI[1] + RI[2])/2; w = (RI[2] - RI[1])/2; M = max(abs(RI[1]), abs(RI[2]));
    bnd = w*sum(i = 1, poldegree(A), i*abs(polcoeff(A, i))*M^(i - 1));
    val = subst(A, y, m);
    if (abs(val) > bnd, return(sign(val)));
    refine());
  error("sign undecided");
}
{
  my(rts = polrootsreal(Pabs), lo, hi);
  lo = floor(rts[1]*10^6)/10^6 - 1/10^6; hi = lo + 3/10^6;
  if (polsturm(Pabs, [lo, hi]) != 1 || subst(Pabs, y, lo)*subst(Pabs, y, hi) >= 0, error("root isolation"));
  RI = [lo, hi];
  for (i = 1, 200, refine());
}
th0 = polrootsreal(Pabs)[1];
num(a) = subst(lift(a), y, th0);
if (sgn(r) <= 0 || sgn(s1) >= 0 || sgn(s2) >= 0, error("embedding"));
T = vector(13, k, twr(TC[k, ]));
c = T[13];
Gm = matrix(n, n, a, b, if (a == b, Mod(1, Pabs), T[CM[a][b] + 1]));

\\ ------------------------------------------------------------------ the reflection in the hyperplane of A
Apts = vector(72, i, i);              \\ the orbit A: the stored points 0..71 (aut116.py, point_orbit_names)
Qpts = [114, 116];                    \\ the orbit Q: the stored points 113, 115
\\ 21 points of A spanning span(A): chosen greedily in floating point, then G[LA,LA] checked invertible exactly
LA = List();
{
  \\ accept a point when its squared distance from the span of the points accepted so far (a Schur complement,
  \\ det ratio) exceeds 1e-6 in floating point; the exact invertibility is checked below
  my(Gf = matrix(n, n, a, b, num(Gm[a, b])), dprev = 1.);
  foreach(Apts, a, if (#LA < 21, my(tst = concat(Vec(LA), [a]), dn = matdet(matrix(#tst, #tst, i, j, Gf[tst[i], tst[j]])));
    if (dn/dprev > 1e-6, listput(LA, a); dprev = dn)));
}
LA = Vec(LA);
if (#LA != 21, error("LA"));
GLA = matrix(21, 21, i, j, Gm[LA[i], LA[j]]);
if (matdet(GLA) == 0, error("G[LA,LA] singular"));
GLAi = GLA^-1;
out("A_spanning_subset_size", #LA);
\\ the A points lie in span(LA): G[a, a'] = G[a, LA] G[LA,LA]^-1 G[LA, a'] for all a, a' in A (rank of G_AA is 21)
VLA = GLAi*matrix(21, n, i, v, Gm[LA[i], v]);
proj(u, v) = sum(i = 1, 21, Gm[u, LA[i]]*VLA[i, v]);
{
  my(ok = 1);
  for (i = 1, 72, for (j = i, 72, if (proj(Apts[i], Apts[j]) != Gm[Apts[i], Apts[j]], ok = 0)));
  if (!ok, error("A does not lie in span(LA)"));
}
out("A_lies_in_span_of_LA", 1);
mn(u, v) = Gm[u, v] - proj(u, v);    \\ <u, n_A><v, n_A>
G2 = Gm;
{
  for (u = 1, n, if (u != Qpts[1] && u != Qpts[2],
    foreach(Qpts, q, my(val = Gm[u, q] - 2*mn(u, q)); G2[u, q] = val; G2[q, u] = val)));
}
\\ the new values on the four pairs between P and Q, and the check that every other entry is an old value
Ppts = [113, 115];
tnew = Set(concat(vector(2, i, vector(2, j, G2[Ppts[i], Qpts[j]]))));
if (#tnew != 2, error("expected two new values"));
t0n = if (sgn(tnew[1] - tnew[2]) < 0, tnew[1], tnew[2]); t1n = if (t0n == tnew[1], tnew[2], tnew[1]);
if (!(sgn(t1n - t0n) > 0 && sgn(T[3] - t1n) > 0), error("order of the new values"));
T2 = concat([t0n, t1n], T[3..13]);
CM2 = vector(n, a, vector(n, b, -1));
{
  for (a = 1, n, for (b = a + 1, n,
    my(k = 0);
    for (j = 1, 13, if (G2[a, b] == T2[j], k = j; break));
    if (!k, error(Str("entry ", a, ",", b, " is not one of the values")));
    CM2[a][b] = k - 1; CM2[b][a] = k - 1));
}
out("t0_new_numeric", num(t0n)); out("t1_new_numeric", num(t1n));
out("t0_new_tower", tow(t0n)~); out("t1_new_tower", tow(t1n)~);
out("t0_new_plus_t1_new_equals_t0_plus_t1", t0n + t1n == T[1] + T[2]);
out("changed_pairs", sum(a = 1, n, sum(b = a + 1, n, CM2[a][b] != CM[a][b])));
out("classes_of_PQ_pairs", vector(2, i, vector(2, j, CM2[Ppts[i]][Qpts[j]])));
out("contacts", sum(a = 1, n, sum(b = a + 1, n, CM2[a][b] == 12)));
\\ the generators of Aut(X*) preserve the class matrix of X'
{
  for (i = 1, #AUTGEN, my(g = AUTGEN[i]);
    for (a = 1, n, for (b = a + 1, n, if (CM2[g[a] + 1][g[b] + 1] != CM2[a][b], error("generator does not preserve X'")))));
}
out("generators_of_AutXstar_preserve_Xprime", #AUTGEN);

\\ ------------------------------------------------------------------ write the data of X'
system(Str("mkdir -p ", OUTDIR));
fv = Str(OUTDIR, "values.gp"); fc = Str(OUTDIR, "class_matrix.gp"); fb = Str(OUTDIR, "basis.gp");
fo = Str(OUTDIR, "orbitals.gp"); fa = Str(OUTDIR, "aut_generators.gp");
system(Str("rm -f ", fv, " ", fc, " ", fb, " ", fo, " ", fa, " ", OUTDIR, "flex.gp ", OUTDIR, "stress.gp"));
write(fv, "\\\\ The values of X' (xprime.gp): as data/values.gp, with t0, t1 replaced by the values t0', t1' on the pairs P x Q.");
write(fv, Str("D1 = ", D1, ";"));
write(fv, Str("D2 = ", D2, ";"));
TC2 = matrix(13, 8, k, j, tow(T2[k])[j]);
write(fv, Str("TC = ", TC2, ";"));
for (a = 1, n, write(fc, Str(CM2[a])));
write(fb, Str(B));
write(fo, "\\\\ The orbits of Aut(X*) on pairs, with the classes of X' (xprime.gp).");
write(fo, "ORBS = List();");
{
  for (j = 1, #ORBS, my(o = ORBS[j], k = CM2[o[5][1][1] + 1][o[5][1][2] + 1]);
    foreach(o[5], q, if (CM2[q[1] + 1][q[2] + 1] != k, error("orbit not constant in X'")));
    write(fo, Str("listput(ORBS, [", k, ", ", o[2], ", ", o[3], ", \"", o[4], "\", ", o[5], "]);")));
}
write(fa, "\\\\ The generators of Aut(X*) (aut116.py), which preserve X' (xprime.gp); the Galois relabellings are in galois.gp.");
write(fa, Str("AUTGEN = ", AUTGEN, ";"));
write(fa, Str("EPSGEN = ", EPSGEN, ";"));
write(fa, Str("CHARGEN = ", CHARGEN, ";"));
out("written", "certify/note/data_xprime/");  \\ 7 Nov: a fixed relative path, never the machine's own
