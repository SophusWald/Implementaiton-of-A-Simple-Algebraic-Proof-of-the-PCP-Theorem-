"""Simulate the PCP: random 3-colourable graph -> honest prover -> verifier.

    python simulate.py                 # 9 vertices (h = 3, m = 2), 5 verifier runs
    python simulate.py --runs 20 --seed 1 --edge-prob 0.7

An honest proof for a properly coloured graph should be accepted every time.
The oracles chi and A0 are also written out in full to a text file
(default proof_tables.txt).
"""

import argparse
import math
import random

from algebra import Params, T, F2Poly, all_vectors, rho_inv
from prover import honest_prover, proof_length_bits
from verifier import verify


def random_3colorable_graph(n, edge_prob, rng):
    """Plant a random colouring, then add each edge between two differently
    coloured vertices with probability edge_prob."""
    coloring = [rng.randrange(3) for _ in range(n)]
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)
             if coloring[u] != coloring[v] and rng.random() < edge_prob]
    return edges, coloring


def element_name(x):
    """Field elements of F_4 as 0, 1, w, w^2 (other fields: the stored integer)."""
    if T == 2:
        return ["0", "1", "w", "w^2"][x.value]   # x = w is the generator
    return str(x.value)


def write_oracle_table(out, name, oracle, k):
    """Write the point oracle `name` : F_q^k x P_3(t, F_2) -> F_2 in full:
    one row per point a, one answer bit per question P."""
    questions = F2Poly.all(3)
    out.write(f"Point oracle {name} : F_{2 ** T}^{k} x P_3({T}, F_2) -> F_2\n")
    out.write(f"{len(all_vectors(k))} points x {len(questions)} questions "
              f"= {len(all_vectors(k)) * len(questions)} bits\n\n")
    out.write("Question j (column j of the answers):\n")
    for j, P in enumerate(questions):
        out.write(f"  P_{j:<3} = {P}\n")
    out.write("\n")
    width = len(", ".join(["w^2"] * k)) + 2
    out.write(f"{'point a':<{width}}  answers {name}[a, P_0], ..., {name}[a, P_{len(questions) - 1}]"
              f"   value (decoded, not part of the proof)\n")
    for a in all_vectors(k):
        answers = [oracle(a, P) for P in questions]
        bits = "".join(str(bit.value) for bit in answers)
        # the answers to the questions z_0, ..., z_{t-1} are the bits rho(value)
        value = rho_inv([oracle(a, F2Poly.variable(i)).value for i in range(T)])
        point = "(" + ", ".join(element_name(x) for x in a) + ")"
        out.write(f"{point:<{width}}  {bits}   {element_name(value)}\n")
    out.write("\n\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h", type=int, default=3, help="|H|")
    parser.add_argument("--m", type=int, default=2, help="vertices are H^m")
    parser.add_argument("--c", type=int, default=2, help="degree-reduction parameter")
    parser.add_argument("--edge-prob", type=float, default=0.5)
    parser.add_argument("--runs", type=int, default=5, help="verifier runs")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--tables", default="proof_tables.txt",
                        help="file to write the oracles chi and A0 to")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    params = Params(h=args.h, m=args.m, c=args.c)
    n = params.max_vertices()

    edges, coloring = random_3colorable_graph(n, args.edge_prob, rng)
    print(f"graph: {n} vertices, {len(edges)} edges")
    print(f"colouring: {coloring}")
    print(f"parameters: q = 2^{T}, h = {params.h}, m = {params.m}, c = {params.c}, "
          f"D = {params.D}, m1 = {params.m1}")
    print(f"proof length: about 2^{math.log2(proof_length_bits(params)):.0f} bits "
          f"(served as oracles, never written out)")

    proof = honest_prover(params, n, edges, coloring)

    with open(args.tables, "w") as out:
        out.write(f"graph: {n} vertices, edges {edges}\ncolouring: {coloring}\n\n")
        write_oracle_table(out, "chi", proof["chi"], params.m)
        write_oracle_table(out, "A0", proof["A0"], 2 * params.m)
    print(f"wrote the oracles chi and A0 to {args.tables}")

    accepted = 0
    for run in range(args.runs):
        accept, failed = verify(params, n, edges, proof, rng)
        accepted += accept
        print(f"run {run + 1}: {'ACCEPT' if accept else 'REJECT  failed: ' + ', '.join(failed)}")
    print(f"accepted {accepted}/{args.runs}")


if __name__ == "__main__":
    main()
