\\ SPDX-License-Identifier: MIT
\\ Copyright (c) 2026 Arhan Canli (Quantavane)
\\
\\ structure.gp -- the structure of the code X* of (22,116), in exact arithmetic in PARI/GP (run by run_gp.py, or
\\ in gp after setting DATA to the data directory). Everything is exact over the field L of degree 8; signs at the
\\ real embedding are decided by a rational isolating interval of the generator (rigorous). Inputs (data/): values.gp,
\\ class_matrix.gp, basis.gp (classes.py), orbitals.gp and aut_generators.gp (aut116.py). Outputs: data/flex.gp and
\\ data/stress.gp (the certificate data checked by check.gp) and lines "KEY = value" on standard output.
\\
\\ (1) The flex. S runs over the symmetric matrices that are +s_j / -s_j on the two halves (under M10) of each of the
\\     12 orbits of pairs that split and are not contacts, and 0 elsewhere. The image of the first-order map
\\     Y -> X Y^T + Y X^T is the set of S with S = W^T S[B,] + S[,B] W - W^T S[B,B] W, W = G_BB^-1 G[B,]. These
\\     equations at sample pairs give a 1-dimensional solution space; the identity is then checked on all pairs.
\\ (2) The invariant stresses: constant on the 14 orbits of contacts, equilibrium Omega~ G[,B] = 0 (Omega~ = Omega -
\\     c diag(row sums)); a basis of the 3-dimensional space, checked at every point; the stress sigma_A supported in
\\     the 72-point orbit A; the polygon of nonnegative invariant stresses of total weight 1 and its vertices.
\\ (3) H_omega(F) = sum_t omega_t h_t (t the 14 contact orbits), from the Gram matrix of F = M G_BB^-1 X_B,
\\     M = S[,B] - 1/2 W^T S[B,B], <f_a,f_b> = (M G_BB^-1 M^T)_ab; H at the vertices, on sigma_A, and on the stress
\\     of largest smallest entry (the max-min stress), which is written to data/stress.gp.
\\ (4) Exact characters: tr R_g = sum_a Q[g(a),a] with Q = G G^+ = W^T (W W^T)^-1 W; the multiplicities of the
\\     trivial character and of eps (kernel M10) in the tangent space, the rotations and the contact space.
\\ (5) Geometry: the normals n_A, n_B of the hyperplanes spanned by A and B, the axis u0, ranks of G_AA and G_BB.

read(Str(DATA, "values.gp"));
read(Str(DATA, "orbitals.gp"));
read(Str(DATA, "aut_generators.gp"));
CM = readvec(Str(DATA, "class_matrix.gp"));
B = readvec(Str(DATA, "basis.gp"))[1];
n = 116; d = 22; NC = 12;
default(realprecision, 60);
T0 = getwalltime();
out(k, v) = print(k, " = ", v);

\\ ------------------------------------------------------------------ the field
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
\\ the root: s1 + s2 with s1, s2 < 0 is the smallest real root; isolate it in a rational interval
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
RI = [-1, 0];
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
if (30*c^2 + 5*c - 1 != 0, error("c"));

\\ ------------------------------------------------------------------ Gram matrix and W
Gm = matrix(n, n, a, b, if (a == b, Mod(1, Pabs), T[CM[a][b] + 1]));
Bi = vector(d, k, B[k] + 1);
GBB = matrix(d, d, i, j, Gm[Bi[i], Bi[j]]);
GBx = matrix(d, n, i, j, Gm[Bi[i], j]);
GBBi = GBB^-1;
W = GBBi*GBx;
out("seconds_gram", (getwalltime() - T0)\1000);
isB = vector(n, a, 0); for (k = 1, d, isB[Bi[k]] = k);
nonB = select(a -> !isB[a], vector(n, a, a));

\\ ------------------------------------------------------------------ (1) the flex
JF = List(); JC = List();
for (j = 1, #ORBS, if (ORBS[j][3] && ORBS[j][1] != NC, listput(JF, j)); if (ORBS[j][1] == NC, listput(JC, j)));
JF = Vec(JF); JC = Vec(JC);
out("flex_unknowns", #JF);
out("contact_orbits", #JC);
\\ Sj as sparse lists, and their B-rows
Sentry(j) = {my(Mx = matrix(n, n)); foreach(ORBS[JF[j]][5], q, Mx[q[1] + 1, q[2] + 1] = q[3]; Mx[q[2] + 1, q[1] + 1] = q[3]); Mx};
SJ = vector(#JF, j, Sentry(j));
\\ residual of S_j at (v, w): S_vw - [sum_k W_kv S_{b_k w} + sum_k S_{v b_k} W_kw - sum_{k,l} W_kv S_{b_k b_l} W_lw]
resid(Sx, v, w) = {
  my(acc = Sx[v, w]);
  for (k = 1, d, if (Sx[Bi[k], w], acc -= W[k, v]*Sx[Bi[k], w]); if (Sx[v, Bi[k]], acc -= Sx[v, Bi[k]]*W[k, w]));
  for (k = 1, d, for (l = 1, d, if (Sx[Bi[k], Bi[l]], acc += W[k, v]*Sx[Bi[k], Bi[l]]*W[l, w])));
  acc;
}
samp = List();
{
  \\ about a thousand pairs: v in 12 points spread over the non-basis points (with the two of Q), w over all of them
  my(vs = Set(concat(vector(10, i, nonB[1 + ((i - 1)*9) % #nonB]), [114, 116])));
  foreach(vs, v, if (!isB[v], foreach(nonB, w, listput(samp, [v, w]))));
}
Rf = matrix(#samp, #JF, i, j, resid(SJ[j], samp[i][1], samp[i][2]));
Kf = matker(Rf);
out("flex_sample_equations", #samp);
out("flex_sample_nullity", #Kf);
if (#Kf != 1, error("flex nullity in the sample"));
sv = Kf[, 1];
\\ normalise: the value on the first unknown orbit (class 0, PQ) is 1
sv = sv/sv[1];
S = sum(j = 1, #JF, sv[j]*SJ[j]);
\\ the full identity S = W^T S[B,] + S[,B] W - W^T S[B,B] W
SB = matrix(d, n, k, w, S[Bi[k], w]); SBB = matrix(d, d, k, l, S[Bi[k], Bi[l]]);
Rfull = S - (W~*SB + SB~*W - W~*SBB*W);
out("flex_identity_holds_on_all_pairs", Rfull == 0);
out("flex_nonzero_pairs", sum(a = 1, n, sum(b = a + 1, n, S[a, b] != 0)));
out("flex_zero_on_A_pairs", prod(j = 1, #ORBS, if (ORBS[j][4] == "AA", prod(i = 1, #ORBS[j][5], my(q = ORBS[j][5][i]); S[q[1] + 1, q[2] + 1] == 0), 1)));
for (j = 1, #JF, out(Str("flex_value_orbit_", JF[j] - 1, "_class_", ORBS[JF[j]][1], "_", ORBS[JF[j]][4]), [num(sv[j]), tow(sv[j])~]));
out("seconds_flex", (getwalltime() - T0)\1000);

\\ ------------------------------------------------------------------ (3a) the Gram matrix of F and the h_t
Mm = SB~ - (1/2)*W~*SBB;
PhiF = Mm*GBBi*Mm~;
nT = #JC;
TY = matrix(n, n);
for (t = 1, nT, foreach(ORBS[JC[t]][5], q, TY[q[1] + 1, q[2] + 1] = t; TY[q[2] + 1, q[1] + 1] = t));
sizes = vector(nT, t, ORBS[JC[t]][2]);
ht = vector(nT, t, 0);
for (a = 1, n, for (b = a + 1, n, my(t = TY[a, b]); if (t, ht[t] += 2*PhiF[a, b] - c*(PhiF[a, a] + PhiF[b, b]))));
normF2 = trace(PhiF);
out("seconds_PhiF", (getwalltime() - T0)\1000);

\\ ------------------------------------------------------------------ (2) invariant stresses
\\ equilibrium rows: for point a and column k: sum_t omega_t sum_{b: type(a,b) = t} (G[b, B_k] - c G[a, B_k])
eqrows(a) = {my(R = matrix(d, nT)); for (b = 1, n, my(t = TY[a, b]); if (t, for (k = 1, d, R[k, t] += Gm[b, Bi[k]] - c*Gm[a, Bi[k]]))); R};
EQ = matconcat(vector(n, a, eqrows(a))~);
Ks = matker(EQ);
out("invariant_stress_dimension", #Ks);
\\ sigma_A: the stress supported on the orbits of contacts inside A
idxA = select(t -> ORBS[JC[t]][4] == "AA", vector(nT, t, t));
notA = select(t -> ORBS[JC[t]][4] != "AA", vector(nT, t, t));
KA = matker(matrix(#notA, #Ks, i, j, Ks[notA[i], j]));
out("stresses_supported_in_A", #KA);
sigA = Ks*KA[, 1];
sigA = sigA/sigA[idxA[1]];
out("sigma_A_on_types", vector(nT, t, lift(sigA[t])));
out("sigma_A_total_weight", sum(t = 1, nT, sizes[t]*sigA[t]));
out("H_on_sigma_A_is_zero", sum(t = 1, nT, sigA[t]*ht[t]) == 0);
Hof(w) = sum(t = 1, nT, w[t]*ht[t]);
tot(w) = sum(t = 1, nT, sizes[t]*w[t]);
\\ polygon vertices: two active constraints among the 14 (omega_t = 0) and total weight 1
V = List();
{
  for (i = 1, nT, for (j = i + 1, nT,
    my(Mx = matrix(3, #Ks, a, b, if (a == 1, Ks[i, b], a == 2, Ks[j, b], sum(t = 1, nT, sizes[t]*Ks[t, b]))), z, w, ok);
    if (matrank(Mx) < 3, next);
    z = matsolve(Mx, [0, 0, 1]~); w = Ks*z;
    ok = 1; for (t = 1, nT, if (sgn(w[t]) < 0, ok = 0; break));
    if (ok, my(new = 1); for (u = 1, #V, if (V[u] == w, new = 0)); if (new, listput(V, w)))));
}
out("polygon_vertices", #V);
{
  for (u = 1, #V, out(Str("vertex_", u, "_zero_types"), apply(t -> t - 1, select(t -> V[u][t] == 0, vector(nT, t, t))));
                out(Str("vertex_", u, "_H_over_normF2"), num(Hof(V[u]))/num(normF2));
                out(Str("vertex_", u, "_H_sign"), sgn(Hof(V[u]))));
}
\\ the max-min stress: maximise m subject to omega_t >= m, total 1; numerically by checking all triples of
\\ equalities, then exactly at the optimal triple
best = [-1, 0];
{
  forsubset([nT, 3], s,
    my(Mx = matrix(3, #Ks, a, b, if (a < 3, Ks[s[a], b] - Ks[s[3], b], sum(t = 1, nT, sizes[t]*Ks[t, b]))), z, w, m);
    if (matrank(Mx) < 3, next);
    z = matsolve(Mx, [0, 0, 1]~); w = Ks*z; m = vecmin(vector(nT, t, num(w[t])));
    if (abs(m - num(w[s[1]])) < 10^-30 && m > best[1], best = [m, w]));
}
WMM = best[2];
mm = vecmin(vector(nT, t, num(WMM[t])));
out("maxmin_min_entry_times_E", mm*4042);
out("maxmin_positive", prod(t = 1, nT, sgn(WMM[t]) > 0));
out("maxmin_types_at_min", apply(t -> t - 1, select(t -> abs(num(WMM[t]) - mm) < 10^-30, vector(nT, t, t))));
HMM = Hof(WMM);
out("maxmin_H_sign", sgn(HMM));
out("maxmin_H_over_normF2", num(HMM)/num(normF2));
\\ check the max-min stress is an equilibrium stress at every point (it is in the kernel by construction)
out("maxmin_equilibrium", EQ*WMM == 0);
\\ (2b) An exact optimality certificate for the max-min problem (maximise the smallest entry over the invariant
\\ stresses of total weight 1), by linear-programming duality: weights y_t > 0 on a set S of orbits where WMM is
\\ smallest, sum y_t = 1, and sum_{t in S} y_t omega_t = lam * tot(omega) for every invariant stress omega. Then every
\\ invariant stress of total weight 1 has smallest entry <= sum y_t omega_t = lam, and WMM attains lam; the maximisers
\\ are exactly the invariant stresses of total weight 1 that are equal to lam on S and >= lam elsewhere (complementary
\\ slackness). This set is computed exactly: its affine span, and if that is a line, its two endpoints.
{
  my(m0, Z0, TOT, cert = 0, S, y, lam, M, kk, rows, A0, dirs, dv, lo = [], hi = [], tl, th, e1, e2);
  m0 = WMM[select(t -> abs(num(WMM[t]) - mm) < 10^-30, vector(nT, t, t))[1]];
  Z0 = select(t -> WMM[t] == m0, vector(nT, t, t));
  TOT = vector(#Ks[1, ], j, sum(t = 1, nT, sizes[t]*Ks[t, j]));
  forvec(v = vector(#Z0, i, [0, 1]), if (cert, break);
    S = select(i -> v[i], vector(#Z0, i, i)); if (!#S, next); S = vector(#S, i, Z0[S[i]]);
    M = matrix(#TOT, #S + 1, j, i, if (i <= #S, Ks[S[i], j], -TOT[j]));
    kk = matker(M);
    if (matsize(kk)[2] == 1,
      y = vector(#S, i, kk[i, 1]); lam = kk[#S + 1, 1];
      my(sy = sum(i = 1, #S, y[i]));
      if (sy != 0, y = y/sy; lam = lam/sy;
        if (prod(i = 1, #S, sgn(y[i]) > 0), cert = 1))));
  if (!cert, error("no dual certificate for the max-min stress"));
  out("maxmin_dual_support", apply(t -> t - 1, S));
  out("maxmin_dual_weights_numeric", vector(#S, i, num(y[i])));
  out("maxmin_dual_bound_equals_value", lam == m0 && tot(WMM) == 1);
  if (lam != m0 || tot(WMM) != 1, error("dual bound"));
  \\ the set of maximisers: omega = Ks z with Ks[S,] z = lam, TOT z = 1, and Ks z >= lam
  rows = matrix(#S + 1, #TOT, i, j, if (i <= #S, Ks[S[i], j], TOT[j]));
  dirs = matker(rows);
  out("maxmin_face_dimension", matsize(dirs)[2]);
  if (matsize(dirs)[2] == 1,
    dv = Ks*dirs[, 1];
    out("maxmin_face_direction_is_sigma_A", matrank(matconcat([dv, sigA])) == 1);
    for (t = 1, nT, my(s = sgn(dv[t])); if (s > 0, lo = concat(lo, [(m0 - WMM[t])/dv[t]]));
                                         if (s < 0, hi = concat(hi, [(m0 - WMM[t])/dv[t]])));
    if (!#lo || !#hi, error("the face of maximisers is unbounded"));
    tl = lo[1]; for (i = 2, #lo, if (sgn(lo[i] - tl) > 0, tl = lo[i]));
    th = hi[1]; for (i = 2, #hi, if (sgn(hi[i] - th) < 0, th = hi[i]));
    if (sgn(th - tl) < 0, error("empty face"));
    e1 = WMM + tl*dv; e2 = WMM + th*dv;
    out("maxmin_face_is_segment", sgn(th - tl) > 0);
    out("maxmin_face_endpoint_1_orbits_at_min", apply(t -> t - 1, select(t -> e1[t] == m0, vector(nT, t, t))));
    out("maxmin_face_endpoint_2_orbits_at_min", apply(t -> t - 1, select(t -> e2[t] == m0, vector(nT, t, t))));
    my(oth = if (e1 == WMM, e2, if (e2 == WMM, e1, 0)));
    out("maxmin_face_WMM_is_an_endpoint", oth != 0);
    if (oth != 0,
      out("maxmin_face_other_endpoint_orbits_at_min", apply(t -> t - 1, select(t -> oth[t] == m0, vector(nT, t, t))));
      out("maxmin_face_other_endpoint_scaled_min1_numeric", vector(nT, t, num(oth[t]/m0)));
      out("maxmin_face_other_endpoint_scaled_min1_tower", vector(nT, t, tow(oth[t]/m0)~))));
  if (matsize(dirs)[2] == 0, out("maxmin_face_is_point", 1));
  out("maxmin_value_tower", tow(m0)~);
}

\\ H depends only on the stresses of the two isolated contacts p1p2 and q1q2: on the invariant stresses,
\\ omega -> (omega_PP, omega_QQ) has rank 2 and kernel R sigma_A, and H(sigma_A) = 0
tPP = select(t -> ORBS[JC[t]][4] == "PP", vector(nT, t, t))[1];
tQQ = select(t -> ORBS[JC[t]][4] == "QQ", vector(nT, t, t))[1];
LM = matrix(2, #Ks, i, j, if (i == 1, Ks[tPP, j], Ks[tQQ, j]));
out("map_to_PP_QQ_rank", matrank(LM));
kz = matker(LM); out("map_kernel_is_sigma_A", #kz == 1 && matrank(matconcat([Ks*kz[, 1], sigA])) == 1);
zP = matsolve(matconcat([LM; kz[, 1]~]), [1, 0, 0]~); zQ = matsolve(matconcat([LM; kz[, 1]~]), [0, 1, 0]~);
alph = Hof(Ks*zP); bet = Hof(Ks*zQ);
out("alpha_sign", sgn(alph)); out("beta_sign", sgn(bet));
out("H_equals_alpha_PP_plus_beta_QQ_on_basis", prod(i = 1, #Ks, Hof(Ks[, i]) == alph*Ks[tPP, i] + bet*Ks[tQQ, i]));
kap = -alph/bet;
out("alpha_numeric", num(alph)); out("beta_numeric", num(bet));
out("alpha_tower", tow(alph)~); out("beta_tower", tow(bet)~);
out("threshold_ratio_QQ_over_PP_numeric", num(kap)); out("threshold_ratio_tower", tow(kap)~);
{
  out("PP_QQ_proportional_types", [apply(t -> t - 1, select(t -> matrank(matrix(2, #Ks, i, j, if (i == 1, Ks[t, j], Ks[tPP, j]))) == 1, vector(nT, t, t))),
                                   apply(t -> t - 1, select(t -> matrank(matrix(2, #Ks, i, j, if (i == 1, Ks[t, j], Ks[tQQ, j]))) == 1, vector(nT, t, t)))]);
}
{
  for (u = 1, #V, out(Str("vertex_", u, "_ratio_QQ_over_PP"), if (V[u][tPP] == 0, "PP = 0", num(V[u][tQQ]/V[u][tPP]))));
}
out("maxmin_ratio_QQ_over_PP", num(WMM[tQQ]/WMM[tPP]));
out("maxmin_ratio_minus_threshold_sign", sgn(WMM[tQQ]/WMM[tPP] - kap));
out("seconds_stress", (getwalltime() - T0)\1000);

\\ ------------------------------------------------------------------ (4) exact characters
P2 = W*W~; Q = W~*P2^-1*W;
gens = vector(#AUTGEN, i, AUTGEN[i]);
pmul(g, h) = vector(n, i, g[h[i] + 1]);
\\ the group, with the values of the characters with kernels M10, S6, PGL(2,9) (homomorphisms: products along words)
GMAP = Map(); GL = List();
{
  my(id = vector(n, i, i - 1), todo = List(), xx, yy, cx);
  mapput(GMAP, id, [1, 1, 1]); listput(GL, id); listput(todo, id);
  while (#todo, xx = todo[#todo]; listpop(todo); cx = mapget(GMAP, xx);
    for (i = 1, #gens, yy = pmul(gens[i], xx);
      my(cy = vector(3, k, CHARGEN[i, k]*cx[k]));
      if (!mapisdefined(GMAP, yy), mapput(GMAP, yy, cy); listput(GL, yy); listput(todo, yy),
        \\ every edge x -> gen*x of the Cayley graph, not only a spanning tree: the values define homomorphisms
        if (mapget(GMAP, yy) != cy, error("the values on the generators do not define homomorphisms")))));
}
out("group_order", #GL);
out("character_values_consistent_on_all_cayley_edges", #GL*#gens);
{
  \\ the kernels: three distinct subgroups of index 2; the first (M10) is the stabiliser of the point 112 (0-based)
  my(ker = vector(3, k, Set(select(g -> mapget(GMAP, g)[k] == 1, Vec(GL)))), st = Set(select(g -> g[113] == 112, Vec(GL))));
  if (vector(3, k, #ker[k]) != [720, 720, 720] || #Set(ker) != 3, error("kernels"));
  if (ker[1] != st, error("the kernel of eps is not the stabiliser of p1"));
  out("character_kernels_index2_distinct_eps_kernel_is_stabiliser_of_p1", 1);
}
trR(g) = sum(a = 1, n, Q[g[a] + 1, a]);
Ecnt = List(); for (a = 1, n, for (b = a + 1, n, if (CM[a][b] == NC, listput(Ecnt, [a - 1, b - 1])))); Ecnt = Vec(Ecnt);
\\ rows: characters 1, eps(M10), eps(S6), eps(PGL); columns: R^22, T, rotations, R^E
mults = matrix(4, 4); chichi = 0;
{
  for (u = 1, #GL, my(g = GL[u], lam = concat([1], mapget(GMAP, g)), g2 = pmul(g, g), t1 = trR(g), t2 = trR(g2), fx, fe, v);
    fx = sum(a = 1, n, g[a] == a - 1);
    fe = sum(i = 1, #Ecnt, my(p = Ecnt[i], q = [g[p[1] + 1], g[p[2] + 1]]); q == p || q == [p[2], p[1]]);
    v = [t1, fx*(t1 - 1), (t1^2 - t2)/2, fe];
    for (k = 1, 4, for (m = 1, 4, mults[k, m] += lam[k]*v[m]));
    chichi += t1^2);
  mults = mults/#GL; chichi = chichi/#GL;
}
out("chi22_norm_<chi,chi>", lift(chichi));
out("mult_in_R22_[1,eps_M10,eps_S6,eps_PGL]", vector(4, k, lift(mults[k, 1])));
out("mult_trivial_[T,rot,E]", vector(3, m, lift(mults[1, m + 1])));
out("mult_eps_[T,rot,E]", vector(3, m, lift(mults[2, m + 1])));
\\ dimension of the space fixed by the stabiliser of one point in each orbit: (1/|G_a|) sum_{g in G_a} tr R_g
{
  my(reps = [ORBS[1][5][1][1], ORBS[1][5][1][2]], nmA = select(j -> ORBS[j][4] == "AA", vector(#ORBS, j, j))[1],
     nmB = select(j -> ORBS[j][4] == "BB", vector(#ORBS, j, j))[1], pts, st, f);
  pts = [ORBS[nmA][5][1][1], ORBS[nmB][5][1][1], reps[1], reps[2]];
  for (i = 1, #pts, st = select(g -> g[pts[i] + 1] == pts[i], Vec(GL)); f = sum(u = 1, #st, trR(st[u]))/#st;
    out(Str("fixed_space_dim_of_stabiliser_of_point_", pts[i]), lift(f)));
}
out("mult_epsS6_[T,rot,E]", vector(3, m, lift(mults[3, m + 1])));
out("mult_epsPGL_[T,rot,E]", vector(3, m, lift(mults[4, m + 1])));
out("seconds_characters", (getwalltime() - T0)\1000);

\\ ------------------------------------------------------------------ (5) geometry
ptsA = Set(concat(vector(#ORBS, j, if (ORBS[j][4] == "AA", concat(vector(#ORBS[j][5], i, [ORBS[j][5][i][1], ORBS[j][5][i][2]])), []))));
ptsB = Set(concat(vector(#ORBS, j, if (ORBS[j][4] == "BB", concat(vector(#ORBS[j][5], i, [ORBS[j][5][i][1], ORBS[j][5][i][2]])), []))));
out("size_A", #ptsA); out("size_B", #ptsB);
GAA = matrix(#ptsA, #ptsA, i, j, Gm[ptsA[i] + 1, ptsA[j] + 1]);
GBBo = matrix(#ptsB, #ptsB, i, j, Gm[ptsB[i] + 1, ptsB[j] + 1]);
out("rank_G_AA", matrank(GAA)); out("rank_G_BB", matrank(GBBo));
\\ the four small points: P = {112, 114}, Q = {113, 115} (0-based); p1 + p2 + q1 + q2 = 0 means its Gram norm is 0
sm = [113, 115, 114, 116];
out("small_points_sum_is_zero", sum(i = 1, 4, sum(j = 1, 4, Gm[sm[i], sm[j]])) == 0);
\\ normal of span(A) inside span(p1 - p2, q1 - q2): coefficients (alpha, beta) with <alpha(p1-p2) + beta(q1-q2), x_a> = 0
nrm(pts) = {my(a = pts[1] + 1, u = Gm[113, a] - Gm[115, a], v = Gm[114, a] - Gm[116, a]);
  foreach(pts, z, if ((Gm[113, z + 1] - Gm[115, z + 1])*v - (Gm[114, z + 1] - Gm[116, z + 1])*u != 0, error("normal")));
  [v, -u]};
cA = nrm(ptsA); cB = nrm(ptsB);
ipv(u, v) = {my(e = [Mod(1, Pabs), Mod(-1, Pabs)], pp = [113, 115], qq = [114, 116]);
  sum(i = 1, 2, sum(j = 1, 2, e[i]*e[j]*(u[1]*v[1]*Gm[pp[i], pp[j]] + u[1]*v[2]*Gm[pp[i], qq[j]] + u[2]*v[1]*Gm[qq[i], pp[j]] + u[2]*v[2]*Gm[qq[i], qq[j]]))) };
out("n_A_orthogonal_to_n_B", ipv(cA, cB) == 0);

\\ ------------------------------------------------------------------ write the certificate data
fl = Str(DATA, "flex.gp"); st = Str(DATA, "stress.gp");
write(fl, "\\\\ The infinitesimal flex of X* as its first-order Gram change S (symmetric, 0 on the diagonal and contacts).");
write(fl, "\\\\ SV[j] (tower coordinates, basis 1,r,s1,r*s1,s2,r*s2,s1*s2,r*s1*s2): the values; SL: [a, b, j, sign], 0-based");
write(fl, "\\\\ points, meaning S[a,b] = S[b,a] = sign*SV[j]. Written by structure.gp.");
write(fl, Str("SV = [", strjoin(vector(#JF, j, strjoin(apply(z -> Str(z), Vec(tow(sv[j]))), ", ")), "; "), "];"));
write(fl, "SL = List();");
for (j = 1, #JF, if (sv[j] != 0, foreach(ORBS[JF[j]][5], q, write(fl, Str("listput(SL, [", q[1], ",", q[2], ",", j, ",", q[3], "]);")))));
write(st, "\\\\ An equilibrium stress of X*, constant on the 14 orbits of contacts: OM[t] (tower coordinates) on type t;");
write(st, "\\\\ EL: [a, b, t] for every contact, 0-based points. The stress of largest smallest entry among the invariant");
write(st, "\\\\ stresses of total weight 1. Written by structure.gp.");
write(st, Str("OM = [", strjoin(vector(nT, t, strjoin(apply(z -> Str(z), Vec(tow(WMM[t]))), ", ")), "; "), "];"));
write(st, "EL = List();");
for (t = 1, nT, foreach(ORBS[JC[t]][5], q, write(st, Str("listput(EL, [", q[1], ",", q[2], ",", t, "]);"))));
out("stress_values_numeric", vector(nT, t, num(WMM[t])));
out("stress_values_tower", vector(nT, t, tow(WMM[t])~));
out("H_maxmin_tower", tow(HMM)~);
out("normF2_numeric", num(normF2));
out("H_maxmin_numeric", num(HMM));
\\ the max-min stress scaled so that its smallest entries are 1, and H for it; the ratio at the vertices
{
  my(tm = select(t -> abs(num(WMM[t]) - mm) < 10^-30, vector(nT, t, t))[1], m0 = WMM[tm]);
  out("stress_scaled_min1_tower", vector(nT, t, tow(WMM[t]/m0)~));
  out("stress_scaled_min1_numeric", vector(nT, t, num(WMM[t]/m0)));
  out("H_scaled_min1_tower", tow(HMM/m0)~);
  out("H_scaled_min1_numeric", num(HMM/m0));
  out("stress_scaled_total_weight_numeric", num(tot(WMM)/m0));
  for (u = 1, #V, if (V[u][tPP] != 0, out(Str("vertex_", u, "_ratio_tower"), tow(V[u][tQQ]/V[u][tPP])~)));
  out("maxmin_ratio_tower", tow(WMM[tQQ]/WMM[tPP])~);
  out("contact_orbit_names", vector(nT, t, Str(ORBS[JC[t]][4], ORBS[JC[t]][2], if (ORBS[JC[t]][3], "s", ""))));
  out("flex_orbit_names", vector(#JF, j, Str(ORBS[JF[j]][1], ":", ORBS[JF[j]][4], ORBS[JF[j]][2])));
}
out("seconds_total", (getwalltime() - T0)\1000);
