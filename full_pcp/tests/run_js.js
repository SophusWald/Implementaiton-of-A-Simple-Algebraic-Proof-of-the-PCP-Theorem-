// Helper for crosscheck.py: reads a scenario as JSON on stdin, builds the
// proof and runs the verifier with the JavaScript port, prints JSON.

"use strict";

const path = require("path");
for (const file of ["algebra.js", "prover.js", "verifier.js"]) {
  require(path.join(__dirname, "..", file));
}
const PCP = globalThis.PCP;

const spec = JSON.parse(require("fs").readFileSync(0, "utf8"));
PCP.setField(spec.t);
const params = new PCP.Params(spec.h, spec.m, spec.c);
const built = spec.strict
  ? PCP.honestProver(params, spec.n, spec.edges, spec.coloring)
  : PCP.cheatingProver(params, spec.n, spec.edges, spec.coloring);

const polynomials = {};
for (const [name, f] of Object.entries(built.polynomials)) polynomials[name] = f.termList();

const rng = new PCP.Random(spec.seed);
const runs = [];
for (let i = 0; i < spec.runs; i += 1) {
  const { accept, failed } = PCP.verify(params, spec.n, spec.edges, built.proof, rng);
  runs.push([accept, failed]);
}
process.stdout.write(JSON.stringify({ polynomials, runs }));
