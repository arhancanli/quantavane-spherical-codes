"""Census version 9: UNTRUSTED helpers for recognition at higher degree (Section 6 of paper_v9, item A1).

Nothing in this module enters a decision. It refines the zero P' of a square system to many digits and looks for exact
algebraic descriptions of its cosine and inner products by lattice reduction (FLINT's LLL). strict9.py checks every
proposal in exact arithmetic; a wrong proposal can only make that check fail.
"""
import math, os, sys

CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.environ.get("QV_FLINTLIB", os.path.expanduser("~/Library/Caches/quantavane-census/flintlib")))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, CERT)
import numpy as np
import flint


# ============================================================================================ refinement
def _scaled_float(num, e, sigma):
    """num * 2^(sigma - e) as a float (num an integer), truncated to 60 bits first (untrusted use only)."""
    b = num.bit_length() if num >= 0 else (-num).bit_length()
    if b <= 60:
        return math.ldexp(num, sigma - e)
    sh = b - 60
    return math.ldexp(num >> sh if num >= 0 else -((-num) >> sh), sigma - e + sh)


def residual_ints(sysm, Xs, Lx, Z, C, K):
    """Exact residual numerators of the square system at (Z, C) / 2^K: a list of (numerator, exponent) pairs, the value
    being numerator / 2^exponent. Xs / 2^Lx is the slice reference (exact dyadic integers)."""
    out = []
    one = 1 << (2 * K)
    for a in range(sysm.n):
        za = Z[a]
        out.append((sum(z * z for z in za) - one, 2 * K))
    CK = C << K
    for (a, b) in sysm.E:
        out.append((sum(x * y for x, y in zip(Z[a], Z[b])) - CK, 2 * K))
    for (s, t) in sysm.ST:
        acc = 0
        for a in range(sysm.n):
            acc += Xs[a][t] * Z[a][s] - Xs[a][s] * Z[a][t]
        out.append((acc, Lx + K))
    return out


def slice_ref(sysm):
    ratios = [[float(v).as_integer_ratio() for v in row] for row in sysm.xs]
    Lx = max(den.bit_length() - 1 for row in ratios for _, den in row)
    Xs = [[num * ((1 << Lx) // den) for num, den in row] for row in ratios]
    return Xs, Lx


def refine_hp(sysm, Z, C, K0, R, K, max_iter=400, log=None):
    """Chord steps z <- z - R F(z) with exact residuals, from the grid 2^-K0 to the grid 2^-K (untrusted). R is the
    double-precision approximate inverse of the census. Returns Z, C on the grid 2^-K and the history of the largest
    residual (as log2)."""
    sh = K - K0
    Z = [[z << sh for z in row] for row in Z]
    C = C << sh
    Xs, Lx = slice_ref(sysm)
    n, d = sysm.n, sysm.d
    hist = []
    stall = 0
    for it in range(max_iter):
        res = residual_ints(sysm, Xs, Lx, Z, C, K)
        lg = max(((abs(nm).bit_length() - e) if nm else -10 ** 9) for nm, e in res)
        hist.append(lg)
        if lg < -(K - 8):
            break
        if len(hist) > 1 and lg > hist[-2] - 1:
            stall += 1
            if stall >= 3:
                break
        else:
            stall = 0
        sigma = -lg                                    # F * 2^sigma has entries of size about 1
        Fs = np.array([_scaled_float(nm, e, sigma) for nm, e in res])
        dz = -(R @ Fs) if not hasattr(R, "dtype") or R.dtype == np.float64 else -(np.asarray(R, dtype=np.float64) @ Fs)
        mx = float(np.max(np.abs(dz)))
        if mx == 0.0 or not math.isfinite(mx):
            break
        E = math.frexp(mx)[1]
        scale = 60 - E
        q = np.rint(np.ldexp(dz, scale)).astype(np.int64)
        shift = K - sigma - scale                      # increment = q * 2^shift on the grid 2^-K
        inc = [int(x) << shift if shift >= 0 else ((int(x) + (1 << (-shift - 1))) >> -shift) for x in q]
        for a in range(n):
            row = Z[a]
            for k in range(d):
                row[k] += inc[a * d + k]
        C += inc[n * d]
    if log is not None:
        log["refine_log2_residuals"] = hist[:3] + ["..."] + hist[-3:] if len(hist) > 8 else hist
        log["refine_iterations"] = len(hist)
    return Z, C, hist[-1]


# ============================================================================================ lattice reduction
def powers_int(num, den_bits, deg, P):
    """round(x^k 2^P) for k = 0..deg, x = num / 2^den_bits, exactly."""
    out = []
    xk = 1
    ek = 0
    for k in range(deg + 1):
        # x^k = xk / 2^ek
        sh = ek - P
        out.append(xk >> sh if sh >= 0 else xk << -sh)
        xk *= num
        ek += den_bits
    return out


def rel_bits(v):
    return sum(math.log2(abs(int(a)) + 1) for a in v)


def lll_relation(vals, P):
    """A short integer relation among numbers given as integers round(x_i 2^P): LLL on the rows (e_i, x_i). Returns the
    coefficient vector and the size of its residual (log2), or None."""
    k = len(vals)
    rows = []
    for i, x in enumerate(vals):
        row = [0] * (k + 1)
        row[i] = 1
        row[k] = int(x)
        rows.append(row)
    L = flint.fmpz_mat(rows).lll()
    best = None
    for r in range(k):
        v = [int(L[r, j]) for j in range(k)]
        if not any(v):
            continue
        resid = int(L[r, k])
        if abs(resid) > (k + 2) * max(abs(a) for a in v) * 4:
            continue                                   # not a relation at this precision
        sz = rel_bits(v)
        if best is None or sz < best[1]:
            best = (v, sz, resid)
    return best


def accept_bits(P):
    """A relation found at precision 2^-P is accepted only if its coefficients use at most 3P/4 bits in all: LLL finds
    a vector of about P bits in a lattice without a relation, so a genuine relation must be much shorter."""
    return min(P - 64, (3 * P) // 4)


def minpoly_lll(num, den_bits, P, maxdeg, probe=False):
    """UNTRUSTED. The minimal polynomial (integer coefficients, constant first) of x = num / 2^den_bits, of degree at
    most maxdeg, found by LLL at precision 2^-P; accepted when its coefficients use at most accept_bits(P) bits in all
    and it vanishes at x to within the precision. Returns (coeffs, bits used) or None."""
    pw = powers_int(num, den_bits, maxdeg, P)
    # probe: one reduction at the largest degree (any relation of lower degree has a multiple there), then factor
    for deg in ([maxdeg] if probe else range(1, maxdeg + 1)):
        best = lll_relation(pw[:deg + 1], P)
        if best is None:
            continue
        v, sz, resid = best
        if sz > accept_bits(P):
            continue
        if abs(resid) > (deg + 2) * max(abs(a) for a in v) * 4:
            continue
        if v[-1] == 0 and not probe:
            continue                                   # found already at a lower degree
        # factor and keep the irreducible factor that vanishes at x
        f = flint.fmpz_poly(v)
        fac = f.factor()[1]
        bestf = None
        for g, _ in fac:
            co = [int(c) for c in g.coeffs()]
            val = sum(c * p for c, p in zip(co, pw))          # ~ g(x) 2^P
            if bestf is None or abs(val) < bestf[1]:
                bestf = (co, abs(val))
        co, val = bestf
        if val > (len(co) + 2) * max(abs(c) for c in co) * 4:
            continue
        if co[-1] < 0:
            co = [-c for c in co]
        return co, sz
    return None


def in_field_lll(gnum, den_bits, theta_pows, P):
    """UNTRUSTED. Rationals a_0..a_{D-1} with g = sum a_k theta^k, by LLL on (g, 1, theta, ..., theta^(D-1)) at
    precision 2^-P; theta_pows are round(theta^k 2^P). Returns a list of (numerator, denominator) or None."""
    sh = den_bits - P
    gi = gnum >> sh if sh >= 0 else gnum << -sh
    if abs(gi) <= 2:
        return [(0, 1)] * len(theta_pows)
    best = lll_relation([gi] + list(theta_pows), P)
    if best is None:
        return None
    v, sz, resid = best
    if v[0] == 0 or sz > accept_bits(P):
        return None
    if abs(resid) > (len(v) + 2) * max(abs(a) for a in v) * 4:
        return None
    m = v[0]
    return [(-a, m) for a in v[1:]]
