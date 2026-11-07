"""Certified distance between the exactly normalised stored code and the configuration with an exact Gram matrix
(paper, version 5). TRUSTED in ball arithmetic (FLINT arb); the rotation is chosen by untrusted code but is exactly
orthogonal (a product of Householder reflections with rational vectors).

For an exact Gram matrix G over K = Q(sqrt D) with G_BB positive definite and G = G_{.B} G_BB^{-1} G_{B.} (checked by the
caller), the configuration P* = G_{.B} L^{-T} (L the Cholesky factor of G_BB) has Gram matrix G; so has P* Q for every
orthogonal Q. The function returns upper bounds for the least-squares distance (sum over the points of |p*_a Q - x_a|^2)^(1/2)
and the largest coordinate distance between P* Q and the exactly normalised stored code x, for an exactly orthogonal
rational Q close to the least-squares rotation.
"""
import math, os, sys
from fractions import Fraction

sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/flintlib"))
import numpy as np
import flint

PREC = 300


def karb(x, Dsqrt):
    """A K element (kfield.K) as an arb ball."""
    return flint.arb(x.a) + flint.arb(x.b) * Dsqrt


def rational_orthogonal(Ot, d, bits=60):
    """An exactly orthogonal rational matrix close to the orthogonal matrix Ot: the product of the Householder
    reflections I - 2 v v^T / (v^T v) of a QR reduction of Ot, every v rounded to the grid 2^-bits (each factor is exactly
    orthogonal whatever the rounding), times the remaining diagonal of signs."""
    A = np.array(Ot, dtype=float)
    Hs = []
    for k in range(d):
        x = A[k:, k].copy()
        v = x.copy()
        v[0] += (1.0 if x[0] >= 0 else -1.0) * np.linalg.norm(x)
        vq = [Fraction(0)] * k + [Fraction(round(t * 2 ** bits), 2 ** bits) for t in v]
        vv = sum(t * t for t in vq)
        if vv == 0:
            continue
        Hs.append((vq, vv))
        vf = np.array([float(t) for t in vq])
        A = A - 2.0 * np.outer(vf, vf @ A) / float(vv)
    Q = flint.fmpq_mat(d, d, [int(i == j) for i in range(d) for j in range(d)])
    for vq, vv in Hs:
        ent = []
        for i in range(d):
            for j in range(d):
                f = Fraction(int(i == j)) - 2 * vq[i] * vq[j] / vv
                ent.append(flint.fmpq(f.numerator, f.denominator))
        Q = Q * flint.fmpq_mat(d, d, ent)
    s = [1 if A[i, i] >= 0 else -1 for i in range(d)]
    return Q * flint.fmpq_mat(d, d, [s[i] if i == j else 0 for i in range(d) for j in range(d)])


def config_from_gram(G, B, pts, D):
    """Ball enclosures of P* = G_{.B} L^{-T} for the points pts (lists of arb, one row per point)."""
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    Ds = flint.arb(D).sqrt()
    d = len(B)
    GBB = [[karb(G[a][b], Ds) for b in B] for a in B]
    L = [[flint.arb(0)] * d for _ in range(d)]
    for j in range(d):
        s = GBB[j][j] - sum((L[j][k] * L[j][k] for k in range(j)), flint.arb(0))
        if not s > 0:
            raise ValueError("Cholesky pivot not provably positive")
        L[j][j] = s.sqrt()
        for i in range(j + 1, d):
            L[i][j] = (GBB[i][j] - sum((L[i][k] * L[j][k] for k in range(j)), flint.arb(0))) / L[j][j]
    rows = []
    for a in pts:
        g = [karb(G[a][b], Ds) for b in B]
        # solve L y = g (forward substitution): p*_a = y (row), since P* L^T = G_{.B}
        y = []
        for i in range(d):
            y.append((g[i] - sum((L[i][k] * y[k] for k in range(i)), flint.arb(0))) / L[i][i])
        rows.append(y)
    flint.ctx.prec = old
    return rows


def distance(G, B, pts, dec, D):
    """Upper bounds of the least-squares and the largest coordinate distance between P* Q (points pts) and the exactly
    normalised stored code (decimal strings dec[a]), Q exactly orthogonal, close to the least-squares rotation."""
    rows = config_from_gram(G, B, pts, D)
    d = len(B)
    old = flint.ctx.prec
    flint.ctx.prec = PREC
    P = np.array([[float(x.mid()) for x in r] for r in rows])
    Xf = np.array([[float(t) for t in dec[a]] for a in pts])
    Xf /= np.linalg.norm(Xf, axis=1)[:, None]
    U, _, Vt = np.linalg.svd(P.T @ Xf)
    O = U @ Vt                                      # P O ~ Xf (untrusted)
    Q = rational_orthogonal(O, d)
    ok = Q.transpose() * Q == flint.fmpq_mat(d, d, [int(i == j) for i in range(d) for j in range(d)])
    Qa = [[flint.arb(Q[j, k]) for k in range(d)] for j in range(d)]
    tot = flint.arb(0)
    mx = 0.0
    for r, a in zip(rows, pts):
        q = [Fraction(t) for t in dec[a]]
        s2 = sum(x * x for x in q)
        nrm = flint.arb(flint.fmpq(s2.numerator, s2.denominator)).sqrt()
        for k in range(d):
            z = flint.arb(0)
            for j in range(d):
                z = z + r[j] * Qa[j][k]
            diff = z - flint.arb(flint.fmpq(q[k].numerator, q[k].denominator)) / nrm
            tot = tot + diff * diff
            mx = max(mx, math.nextafter(float(abs(diff).upper()), math.inf))
    D2 = tot.sqrt()
    flint.ctx.prec = old
    return {"rotation_exactly_orthogonal": bool(ok), "dist_ls_upper": math.nextafter(float(D2.upper()), math.inf),
            "dist_coord_upper": mx}
