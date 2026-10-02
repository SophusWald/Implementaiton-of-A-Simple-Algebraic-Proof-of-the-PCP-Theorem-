"""The honest prover (Algorithm "Completeness Prover" in the appendix).

Input : a graph G = (V, E) with V = H^m, and a proper 3-colouring.
Output: the proof, i.e. eight oracles

    chi,       chi_lines         Had(LDE(Color))
    chi_prime, chi_prime_lines   Had(chi~),  chi~(Y1, Y2) = LDE(Color)(Y1) - LDE(Color)(Y2)
    A0,        A0_lines          Had(sum_i A^(i)(X) Y_i),  A^(i) from MultiDivide(LDE(Color)^3 - 1)
    B0,        B0_lines          Had(sum_i B^(i)(X) Y_i),  B^(i) from MultiDivide(LDE(E) (chi~^3 - 1))

For a polynomial f in k variables the two oracles are

    point oracle  f[a, P]          = P(rho(f(a)))                    a in F_q^k
    lines oracle  f_1[(a, b), w, P] = P(rho(Psi(f(a + X b))[w]))     a, b in F_q^k, w in F_q^{c m1}

with P in P_3(t, F_2). Each oracle is a huge table of bits (see
proof_length_bits), so the prover returns it as a function that computes
any requested entry. The verifier only ever queries entries.
"""

from functools import lru_cache

from algebra import Q, ONE, Poly, num_hadamard_questions, psi_evaluate, rho, edge_lde


def is_proper_coloring(n, edges, coloring):
    return len(coloring) == n and all(coloring[u] != coloring[v] for u, v in edges)


def multi_divide(P, H):
    """MultiDivide(P, Z_H(X_1), ..., Z_H(X_k)): the list Q_1, ..., Q_k with
    P = sum_i Q_i * Z_H(X_i). Asserts that P vanishes on H^k (zero remainder)."""
    quotients = []
    remainder = P
    for i in range(P.n):
        Qi, remainder = remainder.divide_by_vanishing(i, H)
        quotients.append(Qi)
    assert remainder.is_zero(), "P does not vanish on H^k"
    return quotients


def vanishing_certificate(P, H):
    """The polynomial sum_i P^(i)(X) * Y_i in 2k variables (X, Y)."""
    k = P.n
    certificate = Poly(2 * k)
    for i, Qi in enumerate(multi_divide(P, H)):
        Y_i = Poly.variable(2 * k, k + i)
        certificate = certificate + Qi.rename(2 * k, list(range(k))) * Y_i
    return certificate


def hadamard_oracle(f):
    """Had(f): (a, P) -> P(rho(f(a)))."""
    value_at = lru_cache(maxsize=None)(f.evaluate)

    def oracle(a, P):
        return P(rho(value_at(a)))
    return oracle


def lines_oracle(f, params):
    """Had(Psi*(L(f))): ((a, b), w, P) -> P(rho(Psi(f(a + X b))[w]))."""
    line = lru_cache(maxsize=None)(f.restrict_to_line)

    def oracle(a, b, w, P):
        coeffs = line(a, b)
        assert len(coeffs) <= params.D + 1
        return P(rho(psi_evaluate(coeffs, w, params.c, params.m1)))
    return oracle


def honest_prover(params, n, edges, coloring):
    """Build the proof for graph ([n], edges) from a proper 3-colouring
    (a list of colours in {0, 1, 2}, one per vertex)."""
    assert n <= params.max_vertices()
    assert is_proper_coloring(n, edges, coloring), "not a proper 3-colouring"
    m, H = params.m, params.H
    omega = params.omega
    colors = [ONE, omega, omega ** 2]

    # Color : H^m -> {1, omega, omega^2}. Grid points that are not vertices
    # of G are isolated vertices and get colour 1.
    color_values = {params.vertex_point(v): ONE for v in range(params.max_vertices())}
    for v in range(n):
        color_values[params.vertex_point(v)] = colors[coloring[v]]

    chi_hat = Poly.lde(m, H, color_values)                      # LDE(Color)
    A0 = vanishing_certificate(chi_hat ** 3 - ONE, H)

    chi_Y1 = chi_hat.rename(2 * m, list(range(m)))              # LDE(Color)(Y1)
    chi_Y2 = chi_hat.rename(2 * m, list(range(m, 2 * m)))       # LDE(Color)(Y2)
    chi_tilde = chi_Y1 - chi_Y2
    E_hat = edge_lde(params, edges)                             # LDE(E)
    B0 = vanishing_certificate(E_hat * (chi_tilde ** 3 - ONE), H)

    proof = {}
    for name, f in [("chi", chi_hat), ("chi_prime", chi_tilde), ("A0", A0), ("B0", B0)]:
        assert f.degree() <= params.D
        proof[name] = hadamard_oracle(f)
        proof[name + "_lines"] = lines_oracle(f, params)
    return proof


def proof_length_bits(params):
    """Total number of bits in the eight oracle tables.

    A point oracle over F_q^k has q^k * |P_3| entries; a lines oracle has
    q^{2k} * q^{c m1} * |P_3| entries. The four polynomials live in
    k = m, 2m, 2m, 4m variables."""
    m, c, m1 = params.m, params.c, params.m1
    questions = num_hadamard_questions(3)
    total = 0
    for k in (m, 2 * m, 2 * m, 4 * m):
        total += Q ** k * questions                      # point oracle
        total += Q ** (2 * k) * Q ** (c * m1) * questions  # lines oracle
    return total
