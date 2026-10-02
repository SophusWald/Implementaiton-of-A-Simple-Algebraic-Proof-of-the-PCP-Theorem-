# A simple algebraic PCP for 3-COLOR, in Python
The paper
[*A Simple Algebraic Proof of the PCP Theorem*](https://drive.google.com/file/d/1pWzJvBc48sUURmrZ0LVTzeWnPD2oInfn/view?usp=sharing)
(Amireddy, Behera, Srinivasan, Sudan),
presents a simple PCP-protocol which does not rely on the rather abstract notion of composition of PCPs.
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
| `simple_pcp/simulate.py` | Random graph → prover → verifier, plus the printout of `chi` and `A0`. |
| `simple_pcp/algebra.py` | Shared algebra: the field `F_q`, vectors, polynomials over `F_2` and over `F_q`, the map `Psi`, the points `Phi_lambda`, and the parameters. |

## The proof

To gain a full overview of how the protocol works, we refer to the paper 
[*A Simple Algebraic Proof of the PCP Theorem*](https://drive.google.com/file/d/1pWzJvBc48sUURmrZ0LVTzeWnPD2oInfn/view?usp=sharing).
We give a short description here. 

Vertices corresponds to points in a grid `H^m`, and the colours correspond to `1, ω, ω²` in `F_q`, where `ω` is the primitive third root of unity.
The proof consists of eight oracles: a point oracle and a lines oracle for each of four polynomials.

| Polynomial | Variables | Meaning |
|---|---|---|
| `chi` = LDE(Color) | m | the colouring |
| `chi_prime` = chi(Y1) − chi(Y2) | 2m | differences of colours |
| `A0` = Σ A⁽ⁱ⁾(X) Yᵢ | 2m | certificate that chi³ − 1 vanishes on `H^m` |
| `B0` = Σ B⁽ⁱ⁾(X) Yᵢ | 4m | certificate that LDE(E)·(chi_prime³ − 1) vanishes on `H^{2m}` |

For a polynomial `f` in `k` variables:

- the **point oracle** is `f[a, P] = P(ρ(f(a)))` for `a ∈ F_q^k` and `P ∈ P_3(t, F_2)`;
- the **lines oracle** is `f_1[(a, b), w, P] = P(ρ(Ψ(f(a + Xb))[w]))`.

Each oracle is a huge table of bits, so the prover hands it to the
verifier as a function that computes any requested entry. The verifier
only ever queries individual entries.

### Proof length (default parameters: q = 4, h = 3, m = 2, c = 2)

| Oracle | Bits |
|---|---|
| `chi` | 256 |
| `chi_prime`, `A0` | 4,096 each |
| `B0` | 2²⁰ |
| lines oracles | 2³² (`chi`) up to 2⁵⁶ (`B0`) |

The total is about 2⁵⁶ bits, almost all of it the lines oracle of `B0`.

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

## The table printout

`proof_tables.txt` lists, for `chi` and `A0`, the 16 questions `P_0, …, P_15`
(all polynomials in two bits). It then gives one row per point `a`, with the
16 answer bits. The last column is the decoded value `f(a)`, written with
`w = ω`; it is only there for readability and is not part of the proof.
Vertex `v` is the point `(v mod h, v div h)` of `H^2`, with `H = {0, 1, w}`.
