"""Check that the JavaScript port computes exactly what simple_pcp computes.

    python tests/crosscheck.py            # needs node on the PATH

For several fields, graphs and colourings, both implementations build the
proof and run the verifier with the same stream of random numbers. The test
compares the four polynomials chi, chi_prime, A0, B0 term by term, and the
verifier's verdict and list of failed tests on every run.

simple_pcp fixes the field in algebra.py (T = 2). For other fields the test
copies simple_pcp to a temporary directory and changes that one line.
The cheating proofs are built here from simple_pcp's own functions, the same
way prover.js builds them: run MultiDivide and drop the remainder.
"""

import importlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SIMPLE_PCP = os.path.join(HERE, "..", "..", "simple_pcp")
RUN_JS = os.path.join(HERE, "run_js.js")
MASK = 0xFFFFFFFF


class Mulberry32:
    """The random source of algebra.js (class Random), with Python's
    randrange interface, so both implementations draw the same numbers."""

    def __init__(self, seed):
        self.state = seed & MASK

    def next32(self):
        self.state = (self.state + 0x6D2B79F5) & MASK
        t = self.state
        t = ((t ^ (t >> 15)) * (t | 1)) & MASK
        t = (t ^ (t + ((t ^ (t >> 7)) * (t | 61)))) & MASK
        return (t ^ (t >> 14)) & MASK

    def randrange(self, start, stop=None):
        if stop is None:
            start, stop = 0, start
        return start + self.next32() % (stop - start)


def load_simple_pcp(t, workdir):
    """Import algebra, prover, verifier from a copy of simple_pcp with T = t."""
    target = os.path.join(workdir, f"t{t}")
    os.makedirs(target)
    for name in ("algebra.py", "prover.py", "verifier.py"):
        shutil.copy(os.path.join(SIMPLE_PCP, name), target)
    algebra_path = os.path.join(target, "algebra.py")
    with open(algebra_path) as f:
        source = f.read()
    assert source.count("\nT = 2\n") == 1
    with open(algebra_path, "w") as f:
        f.write(source.replace("\nT = 2\n", f"\nT = {t}\n"))
    for name in ("algebra", "prover", "verifier"):
        sys.modules.pop(name, None)
    sys.path.insert(0, target)
    try:
        return [importlib.import_module(name) for name in ("algebra", "prover", "verifier")]
    finally:
        sys.path.remove(target)


def python_proof(algebra, prover, params, n, edges, coloring, strict):
    if strict:
        return prover.honest_prover(params, n, edges, coloring), None
    # The cheating prover: honest_prover without the two assertions.
    Poly, ONE, ZERO = algebra.Poly, algebra.ONE, algebra.ZERO
    m, H = params.m, params.H

    def certificate(P):
        k = P.n
        out, remainder = Poly(2 * k), P
        for i in range(k):
            Qi, remainder = remainder.divide_by_vanishing(i, H)
            out = out + Qi.rename(2 * k, list(range(k))) * Poly.variable(2 * k, k + i)
        return out

    colours = [ONE, params.omega, params.omega ** 2, ZERO]
    values = {params.vertex_point(v): ONE for v in range(params.max_vertices())}
    for v in range(n):
        values[params.vertex_point(v)] = colours[coloring[v]]
    chi = Poly.lde(m, H, values)
    A0 = certificate(chi ** 3 - ONE)
    chi_tilde = chi.rename(2 * m, list(range(m))) - chi.rename(2 * m, list(range(m, 2 * m)))
    B0 = certificate(algebra.edge_lde(params, edges) * (chi_tilde ** 3 - ONE))
    polys = {"chi": chi, "chi_prime": chi_tilde, "A0": A0, "B0": B0}
    proof = {}
    for name, f in polys.items():
        proof[name] = prover.hadamard_oracle(f)
        proof[name + "_lines"] = prover.lines_oracle(f, params)
    return proof, polys


def honest_polys(algebra, prover, params, n, edges, coloring):
    """The polynomials inside honest_prover (which only returns oracles)."""
    _, polys = python_proof(algebra, prover, params, n, edges, coloring, strict=False)
    return polys


def normalise(terms):
    """{exponent tuple: coefficient value} from either implementation."""
    return {tuple(e): c for e, c in terms}


# K4 on vertices 0..3 plus a few edges; colourings as in the demo.
K4_EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
SCENARIOS = [
    # t, h, label, n, edges, colouring, strict, runs
    (2, 3, "honest", 9, [(0, 4), (1, 5), (2, 7), (4, 8), (5, 6)],
     [0, 1, 2, 0, 1, 2, 0, 1, 2], True, 15),
    (2, 3, "closest", 9, K4_EDGES + [(4, 5), (6, 8)], [0, 1, 2, 0, 0, 1, 2, 0, 1], False, 15),
    (4, 2, "four colours", 4, K4_EDGES, [0, 1, 2, 3], False, 15),
    (6, 3, "honest", 9, [(0, 4), (1, 5), (2, 7), (4, 8), (5, 6)],
     [0, 1, 2, 0, 1, 2, 0, 1, 2], True, 6),
    (6, 3, "closest", 9, K4_EDGES + [(4, 5), (6, 8)], [0, 1, 2, 0, 0, 1, 2, 0, 1], False, 8),
    (6, 3, "four colours", 9, K4_EDGES + [(4, 5), (6, 8)], [0, 1, 2, 3, 0, 1, 2, 0, 1], False, 8),
]


def main():
    failures = 0
    with tempfile.TemporaryDirectory() as workdir:
        modules = {}
        for t, h, label, n, edges, coloring, strict, runs in SCENARIOS:
            if t not in modules:
                modules[t] = load_simple_pcp(t, workdir)
            algebra, prover, verifier = modules[t]
            params = algebra.Params(h=h, m=2, c=2)
            seed = 1000 * t + h

            proof, polys = python_proof(algebra, prover, params, n, edges, coloring, strict)
            if polys is None:
                polys = honest_polys(algebra, prover, params, n, edges, coloring)
            py_polys = {name: normalise((e, c.value) for e, c in f.terms.items())
                        for name, f in polys.items()}
            rng = Mulberry32(seed)
            py_runs = [list(verifier.verify(params, n, edges, proof, rng)) for _ in range(runs)]

            spec = dict(t=t, h=h, m=2, c=2, n=n, edges=edges, coloring=coloring,
                        strict=strict, seed=seed, runs=runs)
            out = subprocess.run(["node", RUN_JS], input=json.dumps(spec), capture_output=True,
                                 text=True, check=True)
            js = json.loads(out.stdout)
            js_polys = {name: normalise(terms) for name, terms in js["polynomials"].items()}

            same_polys = js_polys == py_polys
            same_runs = js["runs"] == py_runs
            rejected = sum(not accept for accept, _ in py_runs)
            status = "ok" if same_polys and same_runs else "MISMATCH"
            failures += status != "ok"
            print(f"q = {2 ** t:<3} h = {h}  {label:<13} polynomials "
                  f"{'agree' if same_polys else 'DIFFER'}, {runs} verifier runs "
                  f"{'agree' if same_runs else 'DIFFER'} ({rejected} rejected)  {status}")
            if not same_runs:
                for i, (p, j) in enumerate(zip(py_runs, js["runs"])):
                    if p != j:
                        print(f"   run {i + 1}: python {p}\n          js     {j}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
