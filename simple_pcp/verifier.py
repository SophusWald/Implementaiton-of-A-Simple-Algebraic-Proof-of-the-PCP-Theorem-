"""The PCP verifier for 3-COLOR (Algorithm "The PCP Verifier" in the appendix).

The verifier gets the graph and oracle access to the proof:

    proof["chi"](a, P)                 chi[a, P]                 a in F_q^m
    proof["chi_lines"](a, b, w, P)     chi_1[(a, b), w, P]       a, b in F_q^m, w in F_q^{c m1}
    proof["chi_prime"], proof["chi_prime_lines"]    over F_q^{2m}
    proof["A0"],        proof["A0_lines"]           over F_q^{2m}
    proof["B0"],        proof["B0_lines"]           over F_q^{4m}

where P is a polynomial in P_3(t, F_2) and every answer is a bit. It reads
only a handful of entries and outputs accept or reject.

Notation follows the appendix: field elements, vectors and bits are added
and multiplied with + and *, and x ** k is the power x^k. Since we are in
characteristic 2, "-" is the same as "+".
"""

from algebra import (ZERO, ONE, F2Poly, rho, rho_inv, random_bit, random_element, random_vector,
                     zero_vector, concat, phi, Z_H, edge_lde)


def verify(params, n, edges, proof, rng):
    """One run of the verifier. Returns (accept, list of failed tests)."""
    m, c, m1, H, zeta = params.m, params.c, params.m1, params.H, params.zeta
    assert n <= params.max_vertices()

    # ---- the randomness (red in the paper), shared by all tests -----------
    a, b = random_vector(m, rng), random_vector(m, rng)
    alpha, beta = random_vector(2 * m, rng), random_vector(4 * m, rng)
    u, v = random_vector(c * m1, rng), random_vector(c * m1, rng)
    lam = random_element(rng, nonzero=True)
    P = [F2Poly.random(i, rng) for i in range(4)]       # P_i in P_i(t, F_2)
    R = F2Poly.random(3, rng)
    s = random_bit(rng)
    L = F2Poly.random_linear(rng)                       # L(0) = 0

    # ---- subroutines ------------------------------------------------------
    def SC(g, poly):
        """Self-correction: g[poly + R] - g[R]."""
        return g(poly + R) - g(R)

    def LC(theta, u0):
        """Line correction: - sum_{i=1}^{c+1} theta[u0 + zeta^i v, L]."""
        total = ZERO
        for i in range(1, c + 2):
            total = total - theta(u0 + zeta ** i * v, L)
        return total

    def Lambda(gamma):
        """Lambda_gamma(z) = L(rho(gamma * y^3 - gamma)) where y = rho^{-1}(z)."""
        return F2Poly.interpolate(lambda z: L(rho(gamma * rho_inv(z) ** 3 - gamma)))

    s_poly = F2Poly.constant(s)
    zeros_m = zero_vector(m)
    Z_a = Z_H(H, a)
    Z_ab = Z_H(H, concat(a, b))
    Phi_lam, Phi_0 = phi(lam, c, m1), phi(ZERO, c, m1)

    # theta[w, P]: the lines-table entry of A0 / B0 for the line
    # through `start` in direction `direction`
    def lines_entry(oracle, start, direction):
        return lambda w, poly: oracle(start, direction, w, poly)

    theta_A1 = lines_entry(proof["A0_lines"], concat(a, zeros_m), alpha)
    theta_A2 = lines_entry(proof["A0_lines"], concat(a, Z_a), alpha)
    theta_B1 = lines_entry(proof["B0_lines"], concat(a, b, zeros_m, zeros_m), beta)
    theta_B2 = lines_entry(proof["B0_lines"], concat(a, b, Z_ab), beta)

    failed = []

    def check(name, ok):
        if not ok:
            failed.append(name)

    # ---- low-degree tests -------------------------------------------------
    for name, k in [("chi", m), ("chi_prime", 2 * m), ("A0", 2 * m), ("B0", 4 * m)]:
        Pi, Pi_lines = proof[name], proof[name + "_lines"]
        a2, b2 = random_vector(k, rng), random_vector(k, rng)   # fresh (blue)

        f = lambda poly: Pi(a2, poly)                      # Pi[a', .]
        f1 = lambda w, poly: Pi_lines(a2, b2, w, poly)     # Pi_1[(a', b'), ., .]
        f1_at_u = lambda poly: f1(u, poly)                 # Pi_1[(a', b'), u, .]

        for g_name, g in [("point", f), ("lines", f1_at_u)]:
            # Hadamard test: affine-ness and multiplicativity of the encoding
            for i in range(4):
                check(f"LDT.{name}.{g_name}.affine{i}",
                      SC(g, P[i] + s_poly) == g(P[i]) + s)
            for i in (1, 2):
                check(f"LDT.{name}.{g_name}.mult{i}",
                      SC(g, L * P[i]) == SC(g, L) * SC(g, P[i]))

        # degree-c test:  Pi_1[u, L] + sum_i Pi_1[u + zeta^i v, L] = 0
        check(f"LDT.{name}.degree_c", f1(u, L) - LC(f1, u) == 0)
        # lines table agrees with the point table at a' + lambda b'
        check(f"LDT.{name}.line_vs_point", LC(f1, Phi_lam) == Pi(a2 + lam * b2, L))

    # ---- 3-colouring tests ------------------------------------------------
    for name, theta in [("A1", theta_A1), ("A2", theta_A2),
                        ("B1", theta_B1), ("B2", theta_B2)]:
        g = lambda poly: theta(u, poly)
        for i in (0, 1):
            check(f"ZERO.{name}.affine{i}", SC(g, P[i] + s_poly) == g(P[i]) + s)
        check(f"ZERO.{name}.degree_c", LC(theta, u) == theta(u, L))

    chi, chi_prime, A0, B0 = proof["chi"], proof["chi_prime"], proof["A0"], proof["B0"]

    # A0 is 0 at (a, 0) and equals chi(a)^3 - 1 at (a, Z_H(a))
    check("ZERO.A1.line_vs_point",
          LC(theta_A1, Phi_lam) == A0(concat(a, zeros_m) + lam * alpha, L))
    check("ZERO.A1.is_zero", LC(theta_A1, Phi_0) == 0)
    check("ZERO.A2.line_vs_point",
          LC(theta_A2, Phi_lam) == A0(concat(a, Z_a) + lam * alpha, L))
    check("ZERO.A2.equals_val",
          LC(theta_A2, Phi_0) == SC(lambda poly: chi(a, poly), Lambda(ONE)))

    # chi'(a, b) = chi(a) - chi(b)
    check("CONS.chi_prime", chi_prime(concat(a, b), L) == chi(a, L) - chi(b, L))

    # B0 is 0 at (a, b, 0) and equals LDE(E)(a, b) (chi'(a, b)^3 - 1) at (a, b, Z_H(a, b))
    E_ab = edge_lde(params, edges).evaluate(concat(a, b))
    check("ZERO.B1.line_vs_point",
          LC(theta_B1, Phi_lam) == B0(concat(a, b, zeros_m, zeros_m) + lam * beta, L))
    check("ZERO.B1.is_zero", LC(theta_B1, Phi_0) == 0)
    check("ZERO.B2.line_vs_point",
          LC(theta_B2, Phi_lam) == B0(concat(a, b, Z_ab) + lam * beta, L))
    check("ZERO.B2.equals_prop",
          LC(theta_B2, Phi_0) == SC(lambda poly: chi_prime(concat(a, b), poly), Lambda(E_ab)))

    return (not failed), failed
