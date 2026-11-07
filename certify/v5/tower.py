"""A tower of real quadratic fields over Q with exact arithmetic and one fixed real embedding (paper, version 5).

L = K0(s1, s2), K0 = Q(sqrt D0), s1^2 = Delta1 in K0, s2^2 = Delta2 in K0. An element of K0 is kfield.K (with K.D = D0);
an element of K1 = K0(s1) is a pair (a, b) meaning a + b s1; an element of L = K1(s2) is a pair (a, b) of K1 elements
meaning a + b s2. The representation is unique when [L:K0] = 4, which check_degree() proves (Delta1 and Delta2 are not
squares in K0, and Delta2/Delta1 is not a square in K0: norms to Q that are negative or not rational squares).
The real embedding fixes sqrt(D0) > 0 and the signs sigma1, sigma2 of s1 = sigma1 sqrt(Delta1), s2 = sigma2 sqrt(Delta2)
(Delta1, Delta2 must be positive in this embedding, which check_degree() also proves). Signs of elements are decided in
FLINT ball arithmetic at increasing precision; an element is zero exactly when its representation is zero.
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/flintlib"))
import flint
from kfield import K, fq

ZK = K(0)


class Tower:
    def __init__(self, D0, Delta1, Delta2, sigma1, sigma2):
        K.D = D0
        self.D0, self.Delta1, self.Delta2 = D0, Delta1, Delta2
        self.sigma1, self.sigma2 = sigma1, sigma2
        self._emb = {}

    # ---------------------------------------------------------------- elements: ((a0, a1), (b0, b1)) with K entries
    def el(self, a0=ZK, a1=ZK, b0=ZK, b1=ZK):
        """a0 + a1 s1 + (b0 + b1 s1) s2."""
        return ((a0, a1), (b0, b1))

    def const(self, q):
        return ((K(q), ZK), (ZK, ZK))

    @staticmethod
    def _add1(x, y):
        return (x[0] + y[0], x[1] + y[1])

    @staticmethod
    def _sub1(x, y):
        return (x[0] - y[0], x[1] - y[1])

    def _mul1(self, x, y):
        a, b = x
        c, d = y
        bz, dz = b.iszero(), d.iszero()
        if bz and dz:
            return (a * c, ZK)
        if bz:
            return (a * c, a * d)
        if dz:
            return (a * c, b * c)
        return (a * c + b * d * self.Delta1, a * d + b * c)

    def add(self, x, y):
        return (self._add1(x[0], y[0]), self._add1(x[1], y[1]))

    def sub(self, x, y):
        return (self._sub1(x[0], y[0]), self._sub1(x[1], y[1]))

    def neg(self, x):
        return ((-x[0][0], -x[0][1]), (-x[1][0], -x[1][1]))

    @staticmethod
    def _z1(x):
        return x[0].iszero() and x[1].iszero()

    def iszero(self, x):
        return self._z1(x[0]) and self._z1(x[1])

    def eq(self, x, y):
        return self.iszero(self.sub(x, y))

    def mul(self, x, y):
        a, b = x
        c, d = y
        bz, dz = self._z1(b), self._z1(d)
        if bz and dz:
            return (self._mul1(a, c), (ZK, ZK))
        if bz:
            return (self._mul1(a, c), self._mul1(a, d))
        if dz:
            return (self._mul1(a, c), self._mul1(b, c))
        bd = self._mul1(b, d)
        bdD = (bd[0] * self.Delta2, bd[1] * self.Delta2)
        return (self._add1(self._mul1(a, c), bdD), self._add1(self._mul1(a, d), self._mul1(b, c)))

    def _inv1(self, x):
        a, b = x
        nrm = a * a - b * b * self.Delta1          # in K0
        ni = nrm.inv()
        return (a * ni, -(b * ni))

    def inv(self, x):
        a, b = x
        # (a + b s2)^-1 = (a - b s2) / (a^2 - b^2 Delta2), the denominator in K1
        bb = self._mul1(b, b)
        den = self._sub1(self._mul1(a, a), (bb[0] * self.Delta2, bb[1] * self.Delta2))
        di = self._inv1(den)
        return (self._mul1(a, di), self._mul1((-b[0], -b[1]), di))

    def div(self, x, y):
        return self.mul(x, self.inv(y))

    # ---------------------------------------------------------------- embedding and signs
    def _consts(self, prec):
        if prec not in self._emb:
            old = flint.ctx.prec
            flint.ctx.prec = prec
            r0 = flint.arb(self.D0).sqrt()
            ka = lambda x: flint.arb(x.a) + flint.arb(x.b) * r0
            d1, d2 = ka(self.Delta1), ka(self.Delta2)
            s1 = d1.sqrt() * self.sigma1
            s2 = d2.sqrt() * self.sigma2
            self._emb[prec] = (r0, s1, s2)
            flint.ctx.prec = old
        return self._emb[prec]

    def arb(self, x, prec=200):
        r0, s1, s2 = self._consts(prec)
        old = flint.ctx.prec
        flint.ctx.prec = prec
        ka = lambda y: flint.arb(y.a) + flint.arb(y.b) * r0
        v = ka(x[0][0]) + ka(x[0][1]) * s1 + (ka(x[1][0]) + ka(x[1][1]) * s1) * s2
        flint.ctx.prec = old
        return v

    def sign(self, x):
        if self.iszero(x):
            return 0
        prec = 200
        while prec <= 6400:
            v = self.arb(x, prec)
            if v > 0:
                return 1
            if v < 0:
                return -1
            prec *= 2
        raise ArithmeticError("sign not decided")

    def lt(self, x, y):
        return self.sign(self.sub(x, y)) < 0

    def __float__(self):
        raise TypeError

    def to_float(self, x):
        return float(self.arb(x, 100).mid())

    # ---------------------------------------------------------------- reduction mod p
    def modp_setup(self, p):
        """Square roots r0, r1, r2 mod p of D0, Delta1(r0), Delta2(r0); None if one does not exist."""
        def sqrt_mod(a):
            a %= p
            if a == 0:
                return 0
            if pow(a, (p - 1) // 2, p) != 1:
                return None
            if p % 4 == 3:
                return pow(a, (p + 1) // 4, p)
            # Tonelli-Shanks
            q, s = p - 1, 0
            while q % 2 == 0:
                q //= 2
                s += 1
            z = 2
            while pow(z, (p - 1) // 2, p) != p - 1:
                z += 1
            m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
            while t != 1:
                i, t2 = 0, t
                while t2 != 1:
                    t2 = t2 * t2 % p
                    i += 1
                b = pow(c, 1 << (m - i - 1), p)
                m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
            return r
        r0 = sqrt_mod(self.D0)
        if r0 is None:
            return None
        d1 = self.Delta1.modp(p, r0)
        d2 = self.Delta2.modp(p, r0)
        r1, r2 = sqrt_mod(d1), sqrt_mod(d2)
        if r1 is None or r2 is None:
            return None
        return (p, r0, r1, r2)

    def modp(self, x, setup):
        p, r0, r1, r2 = setup
        m = lambda y: y.modp(p, r0)
        return (m(x[0][0]) + m(x[0][1]) * r1 + (m(x[1][0]) + m(x[1][1]) * r1) * r2) % p

    # ---------------------------------------------------------------- the degree
    def check_degree(self):
        """[L:Q] = 8 and the embedding is real: returns a dict of the exact facts used."""
        def norm(x):          # N_{K0/Q}(a + b sqrt D0) = a^2 - D0 b^2
            return x.a * x.a - self.D0 * x.b * x.b

        def is_rational_square(q):
            if q < 0:
                return False
            num, den = int(q.p), int(q.q)
            import math
            return math.isqrt(num) ** 2 == num and math.isqrt(den) ** 2 == den
        n1, n2 = norm(self.Delta1), norm(self.Delta2)
        n12 = n2 / n1
        facts = {"D0_squarefree_nonsquare": True,
                 "N(Delta1)": str(n1), "N(Delta2)": str(n2), "N(Delta2/Delta1)": str(n12),
                 "Delta1_not_square_in_K0": not is_rational_square(n1),
                 "Delta2_not_square_in_K0": not is_rational_square(n2),
                 "Delta2_over_Delta1_not_square_in_K0": not is_rational_square(n12),
                 "Delta1_positive": K(0) < self.Delta1, "Delta2_positive": K(0) < self.Delta2}
        import math
        facts["D0_squarefree_nonsquare"] = math.isqrt(self.D0) ** 2 != self.D0
        facts["degree_8"] = all(facts[k] for k in ("D0_squarefree_nonsquare", "Delta1_not_square_in_K0",
                                                   "Delta2_not_square_in_K0", "Delta2_over_Delta1_not_square_in_K0"))
        return facts
