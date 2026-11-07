"""Replaces the two floating-point kernel calls of the strictness programs (census/strict.py, paper_v6/symstrict.py) by
the exact tests of exactcheck.py, for the stage that decides strictness (Theorem 3.13: (S') and the span on the box
around the exact configuration). Import this module before running strict.main() or symstrict.main().

  cq.transposed_enclosure  -> square_test(..., mode='box'): DF_E'(p) nonsingular for every p in the box, and the
                              enclosure of the combined stress, exactly (exactcheck.py);
  sphcert.span_check       -> span on the box, exactly (B^-1 in rationals).
Only calls whose caller is strict_entry (strict.py) or full_checks (symstrict.py) are replaced; the census re-run
inside these programs (cq.certify, used only for untrusted data) keeps the original functions. The data needed by the
exact test (the box radius, the extra contacts and their coefficients) are captured from sphcert.System.jac_radius and
cq.extra_rhs, which are called just before.
"""
import hashlib, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exactcheck as X      # sets the library paths
import numpy as np
import cq
import sphcert as S
from fractions import Fraction

SHA = {"exactpatch.py": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
       "exactcheck.py": X.SHA}
_orig_te = cq.transposed_enclosure
_orig_span = S.span_check
_orig_jr = S.System.jac_radius
_orig_rhs = cq.extra_rhs
LAST = {}


def _caller():
    return sys._getframe(2).f_code.co_name


def jac_radius(self, rad):
    self._rad = np.array(rad, dtype=float)
    return _orig_jr(self, rad)


def extra_rhs(extra, beta_f, Xf, radz, N, d):
    LAST["extra"] = (list(extra), np.array(beta_f, dtype=float), np.array(Xf, dtype=float))
    return _orig_rhs(extra, beta_f, Xf, radz, N, d)


def transposed_enclosure(sysm, Jc, Jr, R, b, brad, log, key):
    if _caller() not in ("strict_entry", "full_checks"):
        return _orig_te(sysm, Jc, Jr, R, b, brad, log, key)
    n, d = sysm.n, sysm.d
    Xf = np.asarray(sysm.xs, dtype=float)
    Lc, zi = X.dyadic_exp(Xf)
    Zc = [zi[a * d:(a + 1) * d] for a in range(n)]
    iso = bool(b[-1] == -1.0 and not np.any(b[:-1]))
    bspec = None
    if not iso:
        extra, lam, Xe = LAST["extra"]
        assert np.array_equal(Xe, Xf), "extra_rhs was called at another centre"
        bspec = (extra, lam)
    sq = X.Square(n, d, sysm.E, Xf, Zc, 0, Lc)
    inv = X.Inverse(np.asarray(R, dtype=np.float64))
    res = X.square_test(sq, inv, "box", rad=sysm._rad, bspec=bspec)
    X._ABS.clear()
    log[key + "_beta"] = X.up_float(res["beta_K"])
    log["exact_kernel"] = SHA
    if not res.get("ok_nonsingular"):
        log["exact_box_nonsingular"] = False
        return None, None
    lo = res["stress_lo"]
    log["exact_stress_min_lower"] = X.dn_float(min(lo))
    yt = np.zeros(sysm.nu)
    yt[n:n + sysm.m] = [X.dn_float(x) for x in lo]
    return yt, np.zeros(sysm.nu)


def span_check(Xf, rad, d):
    if _caller() not in ("strict_entry", "full_checks"):
        return _orig_span(Xf, rad, d)
    import scipy.linalg as sl
    import flint
    _, _, piv = sl.qr(np.asarray(Xf).T, pivoting=True, mode="economic")
    idx = sorted(int(i) for i in piv[:d])
    B = flint.fmpq_mat(d, d, [flint.fmpq(*float(Xf[i][k]).as_integer_ratio()) for i in idx for k in range(d)])
    try:
        Bi = B.inv()
    except Exception:
        return False, float("inf"), idx
    nrm = max(sum(abs(Fraction(int(Bi[i, k].p), int(Bi[i, k].q))) for k in range(d)) for i in range(d))
    rs = max(sum(Fraction(float(rad[i][k])) for k in range(d)) for i in idx)
    val = nrm * rs
    return bool(val < 1), X.up_float(val), idx


S.System.jac_radius = jac_radius
cq.extra_rhs = extra_rhs
cq.transposed_enclosure = transposed_enclosure
S.span_check = span_check
