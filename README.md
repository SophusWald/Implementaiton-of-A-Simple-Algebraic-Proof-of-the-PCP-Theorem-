# A simple algebraic PCP for 3-COLOR, in Python
The paper
[*A Simple Algebraic Proof of the PCP Theorem*](https://drive.google.com/file/d/1pWzJvBc48sUURmrZ0LVTzeWnPD2oInfn/view?usp=sharing)
(Amireddy, Behera, Srinivasan, Sudan, Willumsgaard),
presents a simple PCP-protocol for 3-coloring, which does not rely on the rather abstract notion of composition of PCPs,
but instead on the simpler cod-concatenation.
As a result of this, we are able to implement the full algorithm in simple readable python code,
which can both simulate the honest prover and verifier.

We note that this implementation is not optimized for efficiency, but rather for simplicity.

The implementation was primarily written with Claude Code.

## Running it

The implementation is written in Python 3. There are no dependencies.
The verifier.py and prover.py script correspond to the verifier and prover algorithms, 
while the simulate.py script will run both of these instances on a randomly selected graph.

```
cd simple_pcp
python simulate.py                       # 9 vertices, 5 verifier runs
python simulate.py --runs 20 --seed 1 --edge-prob 0.7
python simulate.py --h 2 --m 3           # 8 vertices
```

The simulation:

1. generates a random 3-colourable graph together with a proper 3-colouring.
2. runs the honest prover on it to build the proof.
3. runs the verifier several times. Each run prints
   `ACCEPT`, or `REJECT` together with the tests that failed. An honest
   proof is always accepted.

It also writes the oracles `chi` and `A0` in full to `proof_tables.txt`
(see below).

## Files

| File | Contents |
|---|---|
| `simple_pcp/verifier.py` | `verify(params, n, edges, proof, rng)`: one run of the verifier. Returns accept/reject and the list of failed tests. |
| `simple_pcp/prover.py` | `honest_prover(params, n, edges, coloring)`: builds the proof from a proper 3-colouring. |
| `simple_pcp/simulate.py` | Generates a random 3-colorable graph, and hands it to prover.py which generates a certificate.
verifier.py then checks the certificate and print out `chi` and `A0` in a seperate textfile. |
| `simple_pcp/algebra.py` | Shared algebra: the field `F_q`, vectors, polynomials over `F_2` and over `F_q`, the map `Psi`, the points `Phi_lambda`, and the parameters. |

## The proof

To gain a full overview of how the protocol works, we refer to the paper 
[*A Simple Algebraic Proof of the PCP Theorem*](https://drive.google.com/file/d/1pWzJvBc48sUURmrZ0LVTzeWnPD2oInfn/view?usp=sharing).
We give a short description here. 

Vertices corresponds to points in a grid `H^m`, and the colours correspond to `1, ω, ω²` in `F_q`, where `ω` is the primitive third root of unity, and `q` is a power of 2.
The proof consists of eight oracles: a point oracle and a lines oracle for each of four polynomials.

| Polynomial | Variables | Meaning |
|---|---|---|
| `chi` = LDE(Color) | m | the colouring |
| `chi_prime` = chi(Y1) − chi(Y2) | 2m | differences of colours |
| `A0` = Σ A⁽ⁱ⁾(X) Yᵢ | 2m | certificate that chi³ − 1 vanishes on `H^m` ie. the coloring uses only 3 colours.|
| `B0` = Σ B⁽ⁱ⁾(X) Yᵢ | 4m | certificate that LDE(E)·(chi_prime³ − 1) vanishes on `H^{2m}` ie. every edge has different coloured verticies. |

Each polynomials are encoded in a point oracle and a lines oracle. 
The point oracle is encoded using the Hadamard encoding, 
while the Lines oracle is encoded first with the Reed-Muller encoding and then Hadamard.

Each oracle is a huge table of bits, larger than can be reasonably stored,
so the prover hands it to the verifier as a function that computes any requested entry.
The verifier only ever queries individual entries.

### Proof length (default parameters: q = 4, h = 3, m = 2, c = 2)

| Oracle | Bits |
|---|---|
| `chi` | 256 |
| `A0` | 4,096 |
| `B0` | 2²⁰ |
| lines oracles | 2³² (`chi`) up to 2⁵⁶ (`B0`) |

The total is about 2⁵⁶ bits, almost all of it the lines oracle of `B0`.

### Proof length (default parameters: q = 4, h = 3, m = 2, c = 2)


## Parameters and simplifications

- **The field is fixed at `q = 4`** (`T = 2` in `algebra.py`), to keep the proof short.
  This forces `c = 2` and `|H| ≤ 4`. Completeness holds for every `q`. Soundness,
  however, needs `q` well above the degree bound `D`, so at `q = 4` a cheating
  prover is rarely caught. Setting `T = 4`, `6` or `8` gives a better rejection
  rate and a much longer proof (2¹¹⁹, 2¹⁹⁸ and 2³⁰¹ bits).
- **The degree bound** is `D = 5m(h − 1)`, the exact bound for `B0`, instead of
  the paper's `c₂·hm`.
- **Graphs have exactly `h^m` vertices.** Unused grid points act as isolated
  vertices of colour 1.
