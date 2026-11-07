"""The field K = Q(sqrt D) (D = 1: K = Q) with exact arithmetic over FLINT rationals, for the version 5 programs.

An element is a + b sqrt(D) with a, b in Q (flint.fmpq). Matrices over K are pairs (A0, A1) of flint.fmpq_mat, meaning
A0 + A1 sqrt(D). Everything here is exact; signs of elements are decided exactly.
"""
import os, sys

sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/flintlib"))
import flint

fq = flint.fmpq


class K:
    """a + b sqrt(D), D a squarefree positive integer fixed for the run (K.D); D = 1 means b = 0 always."""
    __slots__ = ("a", "b")
    D = 1

    def __init__(self, a, b=0):
        self.a = a if isinstance(a, fq) else fq(a)
        self.b = b if isinstance(b, fq) else fq(b)

    @staticmethod
    def of(x):
        return x if isinstance(x, K) else K(x)

    def __add__(s, o):
        o = K.of(o)
        return K(s.a + o.a, s.b + o.b)

    __radd__ = __add__

    def __sub__(s, o):
        o = K.of(o)
        return K(s.a - o.a, s.b - o.b)

    def __rsub__(s, o):
        return K.of(o) - s

    def __neg__(s):
        return K(-s.a, -s.b)

    def __mul__(s, o):
        if not isinstance(o, K):
            o = fq(o) if not isinstance(o, fq) else o
            return K(s.a * o, s.b * o)
        return K(s.a * o.a + K.D * s.b * o.b, s.a * o.b + s.b * o.a)

    __rmul__ = __mul__

    def inv(s):
        nrm = s.a * s.a - K.D * s.b * s.b
        if nrm == 0:
            raise ZeroDivisionError("zero element of K")
        return K(s.a / nrm, -s.b / nrm)

    def __truediv__(s, o):
        return s * K.of(o).inv()

    def __eq__(s, o):
        o = K.of(o)
        return s.a == o.a and s.b == o.b

    def __hash__(s):
        return hash((int(s.a.p), int(s.a.q), int(s.b.p), int(s.b.q)))

    def iszero(s):
        return s.a == 0 and s.b == 0

    def sign(s):
        """Exact sign of a + b sqrt(D)."""
        sa = (s.a > 0) - (s.a < 0)
        sb = (s.b > 0) - (s.b < 0)
        if sb == 0 or K.D == 1:
            return sa
        if sa == 0:
            return sb
        if sa == sb:
            return sa
        # opposite signs: compare a^2 with b^2 D
        t = s.a * s.a - K.D * s.b * s.b
        st = (t > 0) - (t < 0)
        return sa * st

    def __lt__(s, o):
        return (s - o).sign() < 0

    def __gt__(s, o):
        return (s - o).sign() > 0

    def __le__(s, o):
        return (s - o).sign() <= 0

    def __ge__(s, o):
        return (s - o).sign() >= 0

    def __float__(s):
        from fractions import Fraction
        a = float(Fraction(int(s.a.p), int(s.a.q)))
        b = float(Fraction(int(s.b.p), int(s.b.q)))
        return a + b * (K.D ** 0.5)

    def __repr__(s):
        if s.b == 0:
            return str(s.a)
        return f"{s.a} + {s.b}*sqrt({K.D})"

    def modp(s, p, r):
        """The image in F_p under sqrt(D) -> r (denominators must be invertible mod p)."""
        def m(q):
            return int(q.p) % p * pow(int(q.q) % p, -1, p) % p
        return (m(s.a) + m(s.b) * r) % p


ZERO, ONE = K(0), K(1)


def ksum(it):
    a, b = fq(0), fq(0)
    for x in it:
        a += x.a
        b += x.b
    return K(a, b)


def prime_for(D, start=(1 << 31) - 1):
    """The largest prime p < start in which D is a square (D = 1: any prime), with a square root r of D mod p."""
    p = start
    while True:
        p -= 1
        if p % 2 == 0:
            continue
        if not flint.fmpz(p).is_prime():
            continue
        if D == 1:
            return p, 1
        # Tonelli-Shanks is not needed for p = 3 mod 4
        if p % 4 == 3:
            r = pow(D, (p + 1) // 4, p)
            if r * r % p == D % p:
                return p, r


# ------------------------------------------------------------------ small dense linear algebra over K (lists of lists)
def mat_mul(A, B):
    n, m, k = len(A), len(B), len(B[0])
    return [[ksum(A[i][t] * B[t][j] for t in range(m) if not A[i][t].iszero() and not B[t][j].iszero())
             for j in range(k)] for i in range(n)]


def ldl_pivots(A):
    """The pivots of A = L D L^T (A symmetric over K), exactly; None if a zero pivot occurs before the end."""
    n = len(A)
    M = [row[:] for row in A]
    piv = []
    for k in range(n):
        if M[k][k].iszero():
            return piv, False
        piv.append(M[k][k])
        inv = M[k][k].inv()
        for i in range(k + 1, n):
            if M[i][k].iszero():
                continue
            f = M[i][k] * inv
            for j in range(k, n):
                if not M[k][j].iszero():
                    M[i][j] = M[i][j] - f * M[k][j]
    return piv, True


def is_pos_def(A):
    """Exact: True iff the symmetric matrix A over K is positive definite (every LDL^T pivot positive)."""
    piv, ok = ldl_pivots(A)
    return ok and all(x.sign() > 0 for x in piv), piv


def solve(A, B):
    """Exact solution X of A X = B over K (A square nonsingular), by Gauss-Jordan elimination."""
    n = len(A)
    M = [A[i][:] + B[i][:] for i in range(n)]
    w = len(M[0])
    for k in range(n):
        p = next((i for i in range(k, n) if not M[i][k].iszero()), None)
        if p is None:
            raise ZeroDivisionError("singular")
        M[k], M[p] = M[p], M[k]
        inv = M[k][k].inv()
        M[k] = [x * inv for x in M[k]]
        for i in range(n):
            if i != k and not M[i][k].iszero():
                f = M[i][k]
                M[i] = [x - f * y for x, y in zip(M[i], M[k])]
    return [row[n:] for row in M]
