# Independent exact audit of We-Solve/research issue 17

This directory independently checks the finite reduced domain at dilation
`k = 5` described in Sections 8.2 and 9.2 of
<https://github.com/We-Solve/research/issues/17>.

It contains a new implementation and new certificates. It does not use the
submitter's unavailable programs, witnesses, or execution logs. No NetworkX,
graph atlas, downloaded graph list, or external numerical solver is used.

## Run

From this directory:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -pedantic audit.cpp -o audit
./audit
python verify.py
```

Do not add `-DNDEBUG` if assertion checking is desired. The recorded run used
g++ 13.3.0 and Python 3.12.14, with C++ assertions enabled.

The first program generates graph representatives and all certificates. The
second program uses only Python's standard library and reads the saved files;
it imports no code from the first program.

## Results from the executed run

| Vertices | Core graphs | Candidate weight vectors | Feasible inputs / checked certificates |
| ---: | ---: | ---: | ---: |
| 5 | 3 | 107 | 107 |
| 6 | 15 | 1,953 | 1,953 |
| 7 | 114 | 40,941 | 40,237 |
| Total | 132 | 43,001 | 42,297 |

Every feasible input has a spanning tree with a positive integer maximum
removal coefficient. The new witnesses have coefficients 1, 2, and 3 in
40,884, 1,410, and 3 cases respectively. These counts describe the first
witness chosen by this implementation. They need not match a different
implementation's histogram, since graph labels and tree enumeration order
affect the first witness.

The C++ program examined 515,910 input/tree pairs before finding those
witnesses. Its recorded runtime was approximately 0.35 seconds; the separate
Python certificate checker took approximately 3.76 seconds. Runtimes are
environment dependent and include the work described in their JSON logs.

## Why the graph generator is exhaustive

Start with the unique zero-vertex graph. At stage `n`, take each generated
isomorphism representative on `n-1` vertices, add a new vertex, and consider
every subset of its old vertices as the new vertex's neighborhood.

Every graph on `n` vertices occurs up to isomorphism in this process: remove
one vertex; the remainder has a representative at the preceding stage by
induction, and the removed vertex's neighborhood appears in the extension
loop. Nothing is filtered out before this generation step.

To remove isomorphic duplicates, first place vertices in nondecreasing
degree order, and then try every permutation within each equal-degree block.
The minimum resulting adjacency bit mask is the canonical representation.
This is exact: any isomorphism preserves degrees, and all remaining
permutations are enumerated. Two graphs have the same canonical mask if and
only if they are isomorphic.

The resulting numbers of unlabeled simple graphs on 0 through 7 vertices
were `1, 1, 2, 4, 11, 34, 156, 1044`. These counts are program outputs; they
are not assumed or hardcoded in the algorithm. Every generated canonical
representative is included in `graph_catalogue.json`.

The core filter keeps exactly the graphs with 5 through 7 vertices,
minimum degree at least 3, connectivity after deleting any one vertex, and

```text
2m <= 5(n-1) <= 3m.
```

The complete resulting edge lists are in `graphs.json`.

## Exact feasibility and maximum checking

For a core with `n` vertices and `m` edges, every reduced input has edge
weights 2 or 3. Exactly `q = 5(n-1)-2m` edges must have weight 3. Both
programs enumerate every choice of those `q` edges.

Membership in `5 P_G` is checked by the total-sum equation, nonnegative
coordinates, and the complete vertex-subset inequalities

```text
w(E(S)) <= 5(|S|-1),  2 <= |S| < n.
```

The C++ search enumerates all spanning trees and minimizes the exact
selected-edge nonnegativity and vertex-subset upper bounds for removal.
Ratios are compared by integer cross multiplication. No floating-point
feasibility or maximality decisions are made.

Each accepted witness is additionally checked directly: the tree is valid,
the coefficient is a positive integer, and the residual satisfies every
coordinate, total-sum, and vertex-subset condition at scale `5-coefficient`.
A blocking constraint then proves maximality. A blocking constraint is
one of:

1. The coefficient has reached the scale bound 5.
2. A selected edge has zero residual, preventing further subtraction.
3. A subset inequality is tight and its bound decreases faster under
   further removal than the subset's residual weight.

The Python checker repeats the domain enumeration, verifies all saved
certificates from direct residual inequalities and blocking constraints,
and checks that there is exactly one certificate for every feasible input
and no certificate for any infeasible or out-of-domain input. It does not
compute or import the optimum-ratio formula.

## Files and encoding

- `audit.cpp`: complete generator, exact search, and direct witness validation.
- `verify.py`: separate standard-library saved-certificate checker.
- `results.json`: generator and search summary.
- `verification.json`: separate checker summary.
- `graph_catalogue.json`: all unlabeled graph representatives through order 7.
- `graphs.json`: explicit core graph edge lists and per-graph counts.
- `witnesses.jsonl`: one JSON record for each of the 42,297 feasible inputs.

Vertices are labeled `0` through `n-1`. Within each core graph, edges are
ordered lexicographically. Bit `i` in a weight mask or tree mask refers to
edge `i` in that order. `canonical_mask` instead uses the complete graph's
lexicographic vertex-pair list, retaining zero bits for absent edges.

Witness records contain the graph identifier, mask of weight-three edges,
tree mask, integer coefficient, and blocking-constraint description.
`blocking_kind` is 0 for the scale bound, 1 for an edge, and 2 for a vertex
subset. For kind 1, `blocking_index` is the zero-based edge index; for kind
2 it is the bit mask of vertices in the blocking subset.

## Scope and limitations

Together with the mathematically proved reductions in the submission,
these finite checks verify the claimed conclusion at `k=5` for every
connected loopless multigraph with at most seven vertices. The code itself
enumerates the reduced simple core domain, rather than infinitely many
multigraphs.

The standalone generator removes dependence on the NetworkX graph atlas.
The separate Python checker establishes certificate correctness and full
weight-domain coverage for the generated core list; graph-isomorphism
completeness is provided by the C++ generator and the induction argument
above. Both new programs were written as part of this review, so the two
implementations are not independent human or formal-proof audits.

This is an ordinary exact computer-assisted verification, not a
proof-assistant verification. It does not resolve arbitrary graphs and
dilations or establish novelty. It also does not authenticate the
submitter's historical execution claims or reproduce its particular
witness ordering, random samples, logs, or file hashes.
