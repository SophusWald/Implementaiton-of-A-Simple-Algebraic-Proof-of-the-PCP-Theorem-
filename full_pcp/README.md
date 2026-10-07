# The PCP for 3-COLOR in the browser

The protocol of `simple_pcp`, ported line by line to JavaScript, with an
interactive page: graph drawing, cheating provers, and 100 verifier runs.

Open `index.html` in a browser. There is no build step. MathJax is loaded
from a CDN, and the rest works offline. To put the page on a website, copy
`index.html` together with the three `.js` files.

## Files

| File | Contents |
|---|---|
| `algebra.js` | Port of `simple_pcp/algebra.py`. The field is chosen at run time with `setField(t)`, for q = 4, 16, 64, 256. |
| `prover.js` | Port of `simple_pcp/prover.py` (`honestProver`), plus `cheatingProver`, `fakeHadamardEncoding` and `fakeZeroSlices`. |
| `verifier.js` | Port of `simple_pcp/verifier.py` (`verify`): the same tests, under the same names. |
| `index.html` | The page. |
| `tests/crosscheck.py` | Checks the port against `simple_pcp` (needs `node`). |

## Cheating provers

A cheating prover builds its proof like the honest prover, but from a false
colouring. Where the honest prover requires MultiDivide to leave remainder
zero, the cheating prover drops the nonzero remainder.

- **closest 3-colouring**: exactly one monochromatic edge. The remainder of the B0 division is dropped, and `ZERO.B2.equals_prop` catches it.
- **closest 4-colouring**: a proper colouring that uses the field element 0 as a fourth colour. The remainder of the A0 division is dropped, and `ZERO.A2.equals_val` catches it.
- **fake Hadamard encoding**: the closest-3-colouring proof, but the table `chi_prime` answers quadratic questions as if the colouring were proper. It stays linear in the question, and its answers to linear questions are honest. The zero test passes, and only `LDT.chi_prime.point.mult*` catches it.
- **functions that are not polynomials**: the closest-3-colouring proof, but A0 and B0 are changed on the slice Y = Z_H(X) to satisfy the zero-test identity there, with the lines table moved to match. The zero test passes, and `ZERO.B2.line_vs_point` catches it.

At q = 64, a single verifier run rejects each of these provers with
probability about 1/2 (the fake Hadamard encoding: about 2/3). At q = 4 it
rejects them much less often, because q is smaller than the degree bound D.
q must be a power of 4: the colours 1, ω, ω² need 3 | q − 1.

## Testing

```
python tests/crosscheck.py
```

For q = 4, 16 and 64, and for honest and cheating colourings, the script runs
both implementations on the same stream of random numbers. It then checks
that `chi`, `chi_prime`, `A0` and `B0` agree term by term, and that every
verifier run gives the same verdict and the same list of failed tests. For
q ≠ 4 it uses a temporary copy of `simple_pcp` with `T` changed.
`simple_pcp` itself is not modified.

## Differences from simple_pcp

- The field size is a parameter, not a constant.
- Polynomials store monomials as packed integers. This limits them to 8
  variables (m ≤ 2) and total degree below 64. The page offers h = 2, 3, 4
  (4, 9, 16 vertices).
- `verify` accepts a precomputed LDE(E) and also returns the names of all
  tests run. Its verdicts are unchanged.
