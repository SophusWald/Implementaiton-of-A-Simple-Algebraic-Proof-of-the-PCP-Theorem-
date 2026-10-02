"""Shared algebra for the simple PCP prover and verifier.

Contents
--------
* F: elements of the field F_q = GF(2^t), and Vec: vectors in F_q^k.
  Both support the usual notation +, -, *, /, ** .
* The bit map rho : F_q -> F_2^t.
* F2Poly: polynomials over F_2 in the t bits of a field element. These are
  the "questions" P in P_3(t, F_2) that a Hadamard-encoded oracle answers.
* Poly: multivariate polynomials over F_q (low-degree extensions, the
  vanishing certificates, restriction to lines).
* The degree-reduction map Psi and the points Phi_lambda.
* Params: all parameters of one PCP instance.

Bits (elements of F_2) are represented by the field elements 0 and 1, so
they can be added and multiplied with + and * as well.
"""

import itertools
from math import comb

# ---------------------------------------------------------------------------
# The field F_q = GF(2^t)
# ---------------------------------------------------------------------------
# The paper takes t = 2c' * ceil(c_1 log(hm)). Here t is a small constant,
# chosen to keep the proof short: t = 2, so q = 4 and F_4 = {0, 1, omega, omega^2}.
# This allows only c = 2 (c + 1 must divide q - 1 = 3) and |H| <= 4.
#
# Completeness does not depend on q, but soundness needs q well above the
# degree D, so with q = 4 a cheating prover is rarely caught. For a better
# rejection rate (and a much longer proof) set T = 4, 6 or 8.
#
# F_q = F_2[x] / (modulus). The element b_0 + b_1 x + ... + b_{t-1} x^{t-1}
# is stored as the integer with binary digits b_{t-1} ... b_0. Python's
# operators ^, & and >> act on binary digits; they are used only in this
# section, to implement the field.

T = 2
Q = 2 ** T
_PRIMITIVE_MODULUS = {          # x generates F_q^* for these moduli
    2: 0b111,                   # x^2 + x + 1
    4: 0x13,                    # x^4 + x + 1
    6: 0x43,                    # x^6 + x + 1
    8: 0x11D,                   # x^8 + x^4 + x^3 + x^2 + 1
}
_MODULUS = _PRIMITIVE_MODULUS[T]

_EXP = [0] * (Q - 1)   # _EXP[i] = (integer of) x^i
_LOG = [0] * Q         # _LOG[x^i] = i
_e = 1
for _i in range(Q - 1):
    _EXP[_i] = _e
    _LOG[_e] = _i
    _e <<= 1                       # times x
    if _e & Q:                     # degree t appeared:
        _e ^= _MODULUS             #   reduce modulo the defining polynomial


class F:
    """An element of F_q, used like a number: x + y, x - y, x * y, x / y, x ** k."""

    __slots__ = ("value",)

    def __init__(self, value):
        assert 0 <= value < Q
        self.value = value

    def __add__(self, other):
        # coefficients of x^i are added mod 2: bitwise XOR of the integers
        return F(self.value ^ other.value)

    __sub__ = __add__              # characteristic 2: x - y = x + y

    def __neg__(self):
        return self

    def __mul__(self, other):
        if not isinstance(other, F):
            return NotImplemented  # e.g. scalar * Vec, handled by Vec
        if self.value == 0 or other.value == 0:
            return ZERO
        return F(_EXP[(_LOG[self.value] + _LOG[other.value]) % (Q - 1)])

    def __truediv__(self, other):
        assert other.value != 0, "division by zero"
        return self * F(_EXP[(-_LOG[other.value]) % (Q - 1)])

    def __pow__(self, k):
        """x ** k for an integer k >= 0, with the convention 0 ** 0 = 1."""
        if k == 0:
            return ONE
        if self.value == 0:
            return ZERO
        return F(_EXP[(_LOG[self.value] * k) % (Q - 1)])

    def __eq__(self, other):
        if isinstance(other, int):           # allows x == 0 and x == 1
            return self.value == other
        return isinstance(other, F) and self.value == other.value

    def __hash__(self):
        return hash(self.value)

    def __bool__(self):
        return self.value != 0

    def __repr__(self):
        return f"F({self.value})"


ZERO, ONE = F(0), F(1)


def element_of_order(k):
    """An element of multiplicative order exactly k (k must divide q - 1)."""
    assert (Q - 1) % k == 0
    return F(2) ** ((Q - 1) // k)        # F(2) is the generator x


def random_bit(rng):
    """A uniformly random bit, as the field element ZERO or ONE."""
    return F(rng.randrange(2))


def random_element(rng, nonzero=False):
    return F(rng.randrange(1 if nonzero else 0, Q))


def rho(x):
    """The F_2-linear bijection F_q -> F_2^t: the t coefficients of x."""
    return tuple((x.value >> i) & 1 for i in range(T))


def rho_inv(bits):
    return F(sum(bit << i for i, bit in enumerate(bits)))


# ---------------------------------------------------------------------------
# Vectors in F_q^k
# ---------------------------------------------------------------------------

class Vec(tuple):
    """A vector in F_q^k: u + v and u - v are coordinatewise, c * v scales."""

    def __add__(self, other):
        assert len(self) == len(other)
        return Vec(x + y for x, y in zip(self, other))

    __sub__ = __add__

    def __rmul__(self, c):
        return Vec(c * x for x in self)


def concat(*vectors):
    """(a, b, ...) as one vector, e.g. concat(a, b) in F_q^{2m}."""
    return Vec(x for v in vectors for x in v)


def zero_vector(k):
    return Vec([ZERO] * k)


def random_vector(k, rng):
    return Vec(random_element(rng) for _ in range(k))


def all_vectors(k):
    """Every point of F_q^k, in a fixed order."""
    return [Vec(p) for p in itertools.product([F(i) for i in range(Q)], repeat=k)]


# ---------------------------------------------------------------------------
# Polynomials over F_2 in the t bits z_0, ..., z_{t-1}
# ---------------------------------------------------------------------------

def _monomials_up_to(degree):
    return [frozenset(s) for k in range(degree + 1)
            for s in itertools.combinations(range(T), k)]


class F2Poly:
    """A multilinear polynomial over F_2 in z_0..z_{t-1}.

    Only the values on F_2^t matter, and z^2 = z there, so multilinear
    polynomials are enough. A monomial is a frozenset of variable indices
    (the empty set is the constant 1). The polynomial is the set of its
    monomials with coefficient 1.
    """

    def __init__(self, monomials=()):
        self.monomials = frozenset(monomials)

    def __add__(self, other):
        # coefficients add mod 2: keep the monomials occurring in exactly one
        return F2Poly(self.monomials.symmetric_difference(other.monomials))

    def __mul__(self, other):
        result = F2Poly()
        for m1 in self.monomials:
            for m2 in other.monomials:
                result = result + F2Poly([m1 | m2])   # z_i * z_i = z_i
        return result

    def __call__(self, bits):
        """Evaluate at a point of F_2^t. The answer is the bit ZERO or ONE."""
        return F(sum(all(bits[i] for i in mono) for mono in self.monomials) % 2)

    @staticmethod
    def constant(s):
        """The constant polynomial s, for a bit s."""
        return F2Poly([frozenset()] if s else [])

    @staticmethod
    def random(degree, rng):
        """Uniformly random element of P_degree(t, F_2)."""
        return F2Poly(m for m in _monomials_up_to(degree) if rng.randrange(2))

    @staticmethod
    def random_linear(rng):
        """Uniformly random L in P_1(t, F_2) with L(0) = 0."""
        return F2Poly(frozenset([i]) for i in range(T) if rng.randrange(2))

    @staticmethod
    def variable(i):
        """The polynomial z_i."""
        return F2Poly([frozenset([i])])

    @staticmethod
    def all(degree):
        """Every element of P_degree(t, F_2), in a fixed order.
        (Only feasible for small t: there are 2^(number of monomials).)"""
        monomials = _monomials_up_to(degree)
        return [F2Poly(mono for mono, bit in zip(monomials, bits) if bit)
                for bits in itertools.product((0, 1), repeat=len(monomials))]

    def __str__(self):
        if not self.monomials:
            return "0"
        ordered = sorted(self.monomials, key=lambda mono: (-len(mono), sorted(mono)))
        return " + ".join("*".join(f"z{i}" for i in sorted(mono)) or "1" for mono in ordered)

    @staticmethod
    def interpolate(func):
        """The multilinear polynomial agreeing with func : F_2^t -> F_2.

        Identify a monomial prod_{i in S} z_i with the set S, and S with its
        indicator vector 1_S in F_2^t. The coefficient of S is the sum of
        func(1_U) over all subsets U of S (Moebius inversion)."""
        subsets = _monomials_up_to(T)
        value = {U: func(tuple(1 if i in U else 0 for i in range(T))) for U in subsets}
        return F2Poly(S for S in subsets
                      if sum((value[U] for U in subsets if U <= S), ZERO) == 1)


def num_hadamard_questions(degree=3):
    """|P_degree(t, F_2)| = 2^(number of monomials of degree <= degree)."""
    return 2 ** sum(comb(T, k) for k in range(degree + 1))


# ---------------------------------------------------------------------------
# Univariate polynomials over F_q (coefficient lists, lowest degree first)
# ---------------------------------------------------------------------------

def uni_mul(p, r):
    out = [ZERO] * (len(p) + len(r) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(r):
            out[i + j] = out[i + j] + x * y
    return out


def vanishing_univariate(H):
    """Z_H(X) = prod_{gamma in H} (X - gamma)."""
    z = [ONE]
    for gamma in H:
        z = uni_mul(z, [-gamma, ONE])
    return z


def lagrange_basis(H, gamma):
    """The univariate polynomial that is 1 at gamma and 0 on H minus {gamma}:
    prod_{other != gamma} (X - other) / (gamma - other)."""
    p = [ONE]
    for other in H:
        if other != gamma:
            d = gamma - other
            p = uni_mul(p, [-other / d, ONE / d])
    return p


def Z_H(H, point):
    """The vector (Z_H(x_1), ..., Z_H(x_k)) for point = (x_1, ..., x_k)."""
    out = []
    for x in point:
        value = ONE
        for gamma in H:
            value = value * (x - gamma)
        out.append(value)
    return Vec(out)


# ---------------------------------------------------------------------------
# Multivariate polynomials over F_q
# ---------------------------------------------------------------------------

class Poly:
    """A polynomial over F_q in n variables, stored as
    {exponent tuple: nonzero coefficient}. Supports P + Q, P - Q, P * Q,
    P ** k, and adding a constant: P + c."""

    def __init__(self, n, terms=None):
        self.n = n
        self.terms = {e: c for e, c in (terms or {}).items() if c != 0}

    @staticmethod
    def constant(n, c):
        return Poly(n, {(0,) * n: c})

    @staticmethod
    def variable(n, j):
        """The polynomial X_j in n variables."""
        return Poly.univariate(n, j, [ZERO, ONE])

    @staticmethod
    def univariate(n, j, coeffs):
        """The polynomial sum_k coeffs[k] X_j^k in n variables."""
        terms = {}
        for k, c in enumerate(coeffs):
            e = [0] * n
            e[j] = k
            terms[tuple(e)] = c
        return Poly(n, terms)

    def __add__(self, other):
        if isinstance(other, F):
            other = Poly.constant(self.n, other)
        terms = dict(self.terms)
        for e, c in other.terms.items():
            terms[e] = terms.get(e, ZERO) + c
        return Poly(self.n, terms)

    __sub__ = __add__              # characteristic 2

    def __mul__(self, other):
        terms = {}
        for e1, c1 in self.terms.items():
            for e2, c2 in other.terms.items():
                e = tuple(x + y for x, y in zip(e1, e2))
                terms[e] = terms.get(e, ZERO) + c1 * c2
        return Poly(self.n, terms)

    def __pow__(self, k):
        result = Poly.constant(self.n, ONE)
        for _ in range(k):
            result = result * self
        return result

    def is_zero(self):
        return not self.terms

    def degree(self):
        return max((sum(e) for e in self.terms), default=0)

    def evaluate(self, point):
        assert len(point) == self.n
        total = ZERO
        for e, c in self.terms.items():
            for x, k in zip(point, e):
                c = c * x ** k
            total = total + c
        return total

    def rename(self, new_n, positions):
        """Move variable i to variable positions[i] in a ring with new_n
        variables (e.g. chi(X) -> chi(Y_2))."""
        terms = {}
        for e, c in self.terms.items():
            new_e = [0] * new_n
            for i, k in enumerate(e):
                new_e[positions[i]] = k
            terms[tuple(new_e)] = c
        return Poly(new_n, terms)

    def divide_by_vanishing(self, i, H):
        """PolyDivide: return (Q, R) with self = Q * Z_H(X_i) + R and
        deg_{X_i}(R) < |H|."""
        z = vanishing_univariate(H)          # monic of degree h
        h = len(H)
        remainder = dict(self.terms)
        quotient = {}
        top = max((e[i] for e in remainder), default=0)
        for d in range(top, h - 1, -1):      # cancel X_i^d, highest d first
            for e in [e for e in remainder if e[i] == d]:
                c = remainder.pop(e)
                q_exp = e[:i] + (d - h,) + e[i + 1:]
                quotient[q_exp] = quotient.get(q_exp, ZERO) + c
                for k in range(h):           # subtract c * X^(d-h) * (z - X^h)
                    r_exp = e[:i] + (d - h + k,) + e[i + 1:]
                    remainder[r_exp] = remainder.get(r_exp, ZERO) - c * z[k]
        return Poly(self.n, quotient), Poly(self.n, remainder)

    def restrict_to_line(self, a, b):
        """The coefficients of the univariate polynomial X -> self(a + X b)."""
        result = [ZERO] * (self.degree() + 1)
        for e, c in self.terms.items():
            p = [c]
            for aj, bj, k in zip(a, b, e):
                for _ in range(k):
                    p = uni_mul(p, [aj, bj])  # times (a_j + b_j X)
            for k, x in enumerate(p):
                result[k] = result[k] + x
        return result

    @staticmethod
    def lde(n, H, values):
        """LDE(f): the unique polynomial of individual degree < |H| that
        agrees with f on H^n, where values = {point of H^n: f(point)}
        (missing points mean f = 0)."""
        total = Poly(n)
        for point, value in values.items():
            if value == 0:
                continue
            term = Poly.constant(n, value)
            for j, gamma in enumerate(point):
                term = term * Poly.univariate(n, j, lagrange_basis(H, gamma))
            total = total + term
        return total


# ---------------------------------------------------------------------------
# Degree reduction (Psi) and the points Phi_lambda
# ---------------------------------------------------------------------------

def base_digits(k, base, num_digits):
    digits = []
    for _ in range(num_digits):
        digits.append(k % base)
        k //= base
    return digits


def psi_evaluate(coeffs, w, c, m1):
    """Psi(p)[w] = sum_k p_k * w_{0,k_0} * ... * w_{c-1,k_{c-1}}, where
    (k_0, ..., k_{c-1}) is k in base m1 and w_{i,j} = w[i*m1 + j]."""
    total = ZERO
    for k, pk in enumerate(coeffs):
        term = pk
        for i, digit in enumerate(base_digits(k, m1, c)):
            term = term * w[i * m1 + digit]
        total = total + term
    return total


def phi(lam, c, m1):
    """Phi_lambda = (1, lam^{m1^i}, ..., lam^{(m1-1) m1^i})_{0 <= i < c}, so
    that Psi(p)[Phi_lambda] = p(lambda). With 0^0 = 1, phi(0) = Phi_0."""
    return Vec(lam ** (j * m1 ** i) for i in range(c) for j in range(m1))


# ---------------------------------------------------------------------------
# Parameters and the graph <-> H^m dictionary
# ---------------------------------------------------------------------------

class Params:
    """Parameters of one PCP instance. Vertices are the points of H^m, so
    graphs have at most h^m vertices."""

    def __init__(self, h=3, m=2, c=2):
        assert (c + 1) % 2 == 1, "c + 1 must be odd (c + 1 = 2^{c'} +- 1)"
        self.h, self.m, self.c = h, m, c
        assert h <= Q, "H must be a subset of F_q"
        self.H = [F(i) for i in range(h)]    # h distinct field elements
        # Degree bound D for all four encoded polynomials. The largest is
        # B0, of degree <= deg(LDE(E)) + 3 deg(chi) = 2m(h-1) + 3m(h-1).
        self.D = 5 * m * (h - 1)
        self.m1 = 2
        while self.m1 ** c <= self.D:        # m1^c > D, so Psi is defined
            self.m1 += 1
        self.omega = element_of_order(3)     # colours are 1, omega, omega^2
        self.zeta = element_of_order(c + 1)

    def max_vertices(self):
        return self.h ** self.m

    def vertex_point(self, v):
        """Vertex v in {0, ..., h^m - 1} -> point of H^m (base-h digits)."""
        return Vec(self.H[d] for d in base_digits(v, self.h, self.m))


def edge_lde(params, edges):
    """LDE(E) on H^{2m}, where E(x, y) = 1 iff {x, y} is an edge."""
    values = {}
    for u, v in edges:
        pu, pv = params.vertex_point(u), params.vertex_point(v)
        values[concat(pu, pv)] = ONE
        values[concat(pv, pu)] = ONE
    return Poly.lde(2 * params.m, params.H, values)
