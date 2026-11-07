\\ SPDX-License-Identifier: MIT
\\ Copyright (c) 2026 Arhan Canli (Quantavane)
\\
\\ check_xprime.gp -- the exact check of check.gp, applied to the second code X' of the note (Section 7): X* with the
\\ two points of Q reflected in the hyperplane spanned by A (data written by xprime.gp into data_xprime/, flex.gp
\\ and stress.gp by structure.gp run on that directory). The same program as check.gp, with two changes: the values
\\ are shown to generate L through c, t7 - t9 = s1 and t3 - t10 = s2 (t1 is not 2 - r/5 here), and the comparison
\\ (N) with the stored code is left out (X' is not in the table). Run: gp -q check_xprime.gp after
\\ DATA = "path/to/data_xprime/" (or QV_NOTE_DATA=data_xprime through run_gp.py). Every claim is decided exactly.
\\
\\ (F) L = Q[y]/(P), P irreducible of degree 8; the thirteen values generate L, increase strictly, and c = t_12 is the
\\     positive root of 30x^2 + 5x - 1.  (G) G is the Gram matrix of 116 unit vectors spanning R^22.  (C) contacts and
\\     no rattlers.  (R) the flexes modulo rotations have dimension at most 1 (rank modulo a prime ideal of degree 1).
\\ (S) the given S is a nonzero flex.  (W) the given omega is a strictly positive equilibrium stress.  (H) the two
\\ formulas for H_omega(F) agree and are positive.  (A) the generators preserve the classes, multiply S by +-1 and fix
\\ omega, and each Galois relabelling (data_xprime/galois.gp) realises its conjugation.  (O) the permutations preserving
\\ the classes, counted by backtracking.

default(realprecision, 60);
n = 116; d = 22;
Tstart = getwalltime();
say(k, v) = print(k, " = ", v);
need(b, msg) = if (!b, error(msg));
\\ used by (O): depth-first search over partial assignments g (0 = unassigned) with candidate bitsets dom
autdfs(dom, g, left) = {
  my(v = 0, best = n + 1, h, bits, u, nd, dead);
  AUTNODES++;
  if (left == 0,
    for (a = 1, n, for (b = a + 1, n, need(CLS[g[a]][g[b]] == CLS[a][b], "complete assignment is not an automorphism")));
    AUTCNT++; return);
  for (w = 1, n, if (!g[w], h = hammingweight(dom[w]); if (h < best, best = h; v = w)));
  if (best == 0, return);
  bits = dom[v];
  while (bits,
    u = valuation(bits, 2) + 1; bits -= 2^(u - 1);
    nd = dom; dead = 0;
    for (w = 1, n, if (!g[w] && w != v, nd[w] = bitand(nd[w], NBRG[CLS[w][v] + 1][u]); if (!nd[w], dead = 1; break)));
    if (!dead, my(g2 = g); g2[v] = u; autdfs(nd, g2, left - 1)));
}

read(Str(DATA, "values.gp"));
CLS = readvec(Str(DATA, "class_matrix.gp"));
BAS0 = readvec(Str(DATA, "basis.gp"))[1];
read(Str(DATA, "flex.gp"));
read(Str(DATA, "stress.gp"));

\\ ------------------------------------------------------------------ (F) the field
\\ theta = s1 + s2; P = Res_x(f1(x), (y - x)^2 - D2(r(x))), f1 the minimal polynomial of s1
f1 = (x^2 - D1[1])^2 - 145*D1[2]^2;
P = polresultant(f1, (y - x)^2 - (D2[1] + D2[2]*(x^2 - D1[1])/D1[2]), x);
P = P/pollead(P);
need(poldegree(P) == 8 && polisirreducible(P), "P is not irreducible of degree 8");
th = Mod(y, P);
{
  my(h = gcd(f1, (th - x)^2 - (D2[1] + D2[2]*(x^2 - D1[1])/D1[2])));
  need(poldegree(h) == 1, "s1 not determined");
  sA = -polcoeff(h, 0)/polcoeff(h, 1);
}
sB = th - sA; rr = (sA^2 - D1[1])/D1[2];
need(rr^2 == 145 && sA^2 == D1[1] + D1[2]*rr && sB^2 == D2[1] + D2[2]*rr, "relations");
say("field_degree", poldegree(P));
say("minimal_polynomial_of_s1_plus_s2", P);
say("real_roots_of_P", polsturm(P));
\\ isolating interval of the real root theta0 = s1 + s2 (s1, s2 < 0): the smallest real root
IV = [0, 0];
bis() = {my(m = (IV[1] + IV[2])/2, v = subst(P, y, m)); need(v != 0, "rational root");
  if (sign(v) == sign(subst(P, y, IV[1])), IV[1] = m, IV[2] = m);}
{
  my(z = polrootsreal(P)[1], lo = floor(z*10^8)/10^8 - 10^-8, hi = lo + 3*10^-8);
  need(polsturm(P, [lo, hi]) == 1 && sign(subst(P, y, lo)) != sign(subst(P, y, hi)), "isolation");
  IV = [lo, hi];
  for (i = 1, 256, bis());
}
\\ sign at theta0: |A(theta0) - A(m)| <= (half width) * sum_i i |a_i| R^(i-1), R >= |theta0|
sgnL(a) = {
  my(A = lift(a), m, h, R, bound, v);
  if (type(A) != "t_POL", return(sign(A)));
  for (it = 1, 2000,
    m = (IV[1] + IV[2])/2; h = (IV[2] - IV[1])/2; R = max(abs(IV[1]), abs(IV[2]));
    bound = h*sum(i = 1, poldegree(A), i*abs(polcoeff(A, i))*R^(i - 1));
    v = subst(A, y, m);
    if (abs(v) > bound, return(sign(v)));
    bis());
  error("sign undecided");
}
need(sgnL(rr) > 0 && sgnL(sA) < 0 && sgnL(sB) < 0, "embedding");
th0 = (IV[1] + IV[2])/2.;
toR(a) = subst(lift(a), y, th0);
BASIS = [1, rr, sA, rr*sA, sB, rr*sB, sA*sB, rr*sA*sB];
el(v) = sum(j = 1, 8, v[j]*BASIS[j]);
TV = vector(13, k, el(TC[k, ]));
c = TV[13];
need(30*c^2 + 5*c - 1 == 0, "c");
need(sgnL(c) > 0, "c > 0");
for (k = 1, 12, need(sgnL(TV[k + 1] - TV[k]) > 0, "values not increasing"));
need(rr == 60*c + 5 && TV[8] - TV[10] == sA && TV[4] - TV[11] == sB, "generation");
say("values_increasing_and_generate_L", 1);
say("cosine", toR(c));

\\ ------------------------------------------------------------------ (G) the Gram matrix
need(#CLS == n && vector(n, a, #CLS[a]) == vector(n, a, n), "class matrix shape");
for (a = 1, n, need(CLS[a][a] == -1, "diagonal"); for (b = 1, n, if (b != a, need(CLS[a][b] == CLS[b][a] && CLS[a][b] >= 0 && CLS[a][b] <= 12, "classes"))));
G = matrix(n, n, a, b, if (a == b, Mod(1, P), TV[CLS[a][b] + 1]));
Bv = vector(d, k, BAS0[k] + 1);
need(#Set(Bv) == d, "basis");
GB = matrix(d, d, i, j, G[Bv[i], Bv[j]]);
for (k = 1, d, need(sgnL(matdet(matrix(k, k, i, j, GB[i, j]))) > 0, Str("leading minor ", k)));
say("G_BB_positive_definite", 1);
GBi = GB^-1;
W = GBi*matrix(d, n, i, j, G[Bv[i], j]);
need(matrix(n, d, a, k, G[a, Bv[k]])*W == G, "G is not G[,B] G_BB^-1 G[B,]");
say("gram_rank_22_identity", 1);
say("seconds_gram", (getwalltime() - Tstart)\1000);

\\ ------------------------------------------------------------------ (C) contacts
isC(a, b) = a != b && CLS[a][b] == 12;
deg = vector(n, a, sum(b = 1, n, isC(a, b)));
say("contacts", sum(a = 1, n, deg[a])/2);
say("min_contacts_per_point", vecmin(deg));
say("max_contacts_per_point", vecmax(deg));
need(vecmin(deg) >= 1, "rattler");

\\ ------------------------------------------------------------------ (R) rank bound modulo a prime ideal of degree 1
inB = vector(n, a, 0); for (k = 1, d, inB[Bv[k]] = k);
NB = select(a -> !inB[a], vector(n, a, a));
\\ unknowns: S[a,b], a < b, {a,b} meets B, not a contact
UID = matrix(n, n); NU = 0;
for (a = 1, n, for (b = a + 1, n, if ((inB[a] || inB[b]) && !isC(a, b), NU++; UID[a, b] = NU; UID[b, a] = NU)));
\\ equations: pairs v <= w outside B with v = w or a contact
EQS = List(); for (i = 1, #NB, for (j = i, #NB, my(v = NB[i], w = NB[j]); if (v == w || isC(v, w), listput(EQS, [v, w]))));
say("rank_unknowns", NU);
say("rank_equations", #EQS);
{
  my(p = 2^30, rho, ok = 0);
  while (!ok, p = nextprime(p + 1);
    if (denominator(content(P)) % p == 0, next);
    my(rts = polrootsmod(P*Mod(1, p))); if (#rts == 0, next);
    rho = lift(rts[1]);
    ok = 1;
    for (k = 1, d, for (a = 1, n, if (denominator(content(lift(W[k, a]))) % p == 0, ok = 0)));
  );
  PR = p; RHO = Mod(rho, p);
}
red(a) = subst(lift(a), y, RHO);
need(subst(P, y, RHO) == 0, "root mod p");
Wp = matrix(d, n, k, a, red(W[k, a]));
say("prime", PR);
{
  my(Mr = matrix(#EQS, NU, i, j, Mod(0, PR)));
  for (i = 1, #EQS, my(v = EQS[i][1], w = EQS[i][2], u);
    for (k = 1, d,
      u = UID[v, Bv[k]]; if (u, Mr[i, u] += Wp[k, w]);
      u = UID[Bv[k], w]; if (u, Mr[i, u] += Wp[k, v]);
      for (l = 1, d, if (l != k, u = UID[Bv[k], Bv[l]]; if (u, Mr[i, u] -= Wp[k, v]*Wp[l, w])))));
  RK = matrank(Mr);
}
say("rank_mod_p", RK);
need(RK >= NU - 1, "rank bound");
say("flex_space_mod_rotations_dim_at_most_1", 1);
say("seconds_rank", (getwalltime() - Tstart)\1000);

\\ ------------------------------------------------------------------ (S) the flex
SVL = vector(matsize(SV)[1], j, el(SV[j, ]));
S = matrix(n, n, a, b, Mod(0, P));
{
  foreach(SL, q, my(a = q[1] + 1, b = q[2] + 1);
    need(a != b && !isC(a, b) && S[a, b] == 0, "flex entry"); S[a, b] = q[4]*SVL[q[3]]; S[b, a] = S[a, b]);
}
need(S != 0, "S = 0");
SB = matrix(d, n, k, b, S[Bv[k], b]);
SBB = matrix(d, d, k, l, S[Bv[k], Bv[l]]);
need(S == W~*SB + SB~*W - W~*SBB*W, "S is not a first-order Gram change");
say("flex_identity_on_all_pairs", 1);
say("flex_nonzero_pairs", #SL);
say("seconds_flex", (getwalltime() - Tstart)\1000);

\\ ------------------------------------------------------------------ (W) the stress
OMV = vector(matsize(OM)[1], t, el(OM[t, ]));
for (t = 1, #OMV, need(sgnL(OMV[t]) > 0, "stress not positive"));
Om = matrix(n, n, a, b, Mod(0, P)); cnt = 0;
{
  foreach(EL, q, my(a = q[1] + 1, b = q[2] + 1); need(isC(a, b) && Om[a, b] == 0, "stress entry");
    Om[a, b] = OMV[q[3]]; Om[b, a] = Om[a, b]; cnt++);
}
need(cnt == sum(a = 1, n, deg[a])/2, "stress does not cover every contact");
Omt = Om - c*matdiagonal(vector(n, a, sum(b = 1, n, Om[a, b])));
need(Omt*matrix(n, d, a, k, G[a, Bv[k]]) == 0, "equilibrium");
say("stress_positive_on_every_contact", 1);
say("stress_equilibrium", 1);

\\ ------------------------------------------------------------------ (H) the second-order form, two ways
PP = W*W~; MM = PP^-1*GBi*PP^-1;
Z = W*S;
H1 = trace(MM*(Z*Omt*Z~));
M = SB~ - (1/2)*W~*SBB;
Phi = M*GBi*M~;
H2 = sum(a = 1, n, sum(b = 1, n, if (Omt[a, b] != 0, Omt[a, b]*Phi[b, a], 0)));
need(H1 == H2, "the two formulas for H differ");
say("H_two_formulas_agree", 1);
HS = sgnL(H1);
say("H_sign", HS);
say("H_numeric", toR(H1));
say("H_tower_coordinates", lift(matsolve(matrix(8, 8, i, j, polcoeff(lift(BASIS[j]), i - 1)), vector(8, i, polcoeff(lift(H1), i - 1))~))~);
need(HS > 0, "H not positive");

\\ ------------------------------------------------------------------ (A) the symmetry data
read(Str(DATA, "aut_generators.gp"));
read(Str(DATA, "galois.gp"));          \\ the Galois relabellings of X' (galois_xprime.py)
need(#GALPERM == 3, "the three Galois relabellings of X' are missing");
{
  for (i = 1, #AUTGEN, my(g = AUTGEN[i]);
    need(#Set(g) == n && vecmin(g) == 0 && vecmax(g) == n - 1, "not a permutation");
    for (a = 1, n, for (b = a + 1, n,
      need(CLS[g[a] + 1][g[b] + 1] == CLS[a][b], "generator does not preserve the classes");
      need(S[g[a] + 1, g[b] + 1] == EPSGEN[i]*S[a, b], "S not eps-equivariant");
      need(Om[g[a] + 1, g[b] + 1] == Om[a, b], "omega not invariant"))));
}
say("generators_preserve_classes_S_up_to_sign_and_omega", #AUTGEN);
need(EPSGEN == CHARGEN[, 1]~, "EPSGEN differs from the first column of CHARGEN");
\\ the conjugations of L over Q(r): the class maps agree with the values, and the relabellings realise them
{
  my(flip = [[1, 1, -1, -1, 1, 1, -1, -1], [1, 1, 1, 1, -1, -1, -1, -1], [1, 1, -1, -1, -1, -1, 1, 1]]);
  for (j = 1, #GALPERM, my(cm = vector(13, k, k - 1), pi = GALPERM[j][2], mp = GALPERM[j][1]);
    for (i = 1, matsize(mp)[1], cm[mp[i, 1] + 1] = mp[i, 2]);
    for (k = 1, 13, need(TC[cm[k] + 1, ] == vector(8, i, flip[j][i]*TC[k, i]), "class map is not the conjugation"));
    need(#Set(pi) == n, "relabelling");
    for (a = 1, n, for (b = a + 1, n, need(CLS[pi[a] + 1][pi[b] + 1] == cm[CLS[a][b] + 1], "relabelling does not realise the conjugation"))));
}
say("galois_conjugations_realised_by_relabellings", #GALPERM);

\\ ------------------------------------------------------------------ (O) the order of the group, independently
\\ Points are given images one at a time (the unassigned point with fewest candidates first). The candidates of a point
\\ are kept as a bitset: at the start, the points whose row of the class matrix has the same multiset of classes; after
\\ assigning v -> u, those of every unassigned w are cut down to the points u' != u with class(u', u) = class(w, v).
\\ Every permutation preserving the classes survives every cut, so it is reached exactly once as a complete assignment;
\\ every complete assignment is checked on all pairs. So the count is the order of the group.
{
  my(NBR = vector(13, k, vector(n, c, 0)), PROF = vector(n, a, vecsort(CLS[a])), dom0);
  for (c = 1, n, for (u = 1, n, if (u != c, NBR[CLS[u][c] + 1][c] += 2^(u - 1))));
  NBRG = NBR; AUTCNT = 0; AUTNODES = 0;
  dom0 = vector(n, v, sum(u = 1, n, if (PROF[u] == PROF[v], 2^(u - 1), 0)));
  autdfs(dom0, vector(n, i, 0), n);
}
say("automorphisms_by_backtracking", AUTCNT);
say("backtracking_nodes", AUTNODES);
need(AUTCNT == 1440, "the group does not have order 1440");
say("seconds_total", (getwalltime() - Tstart)\1000);
say("STATUS", "ALL CHECKS PASSED");
