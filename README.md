# Isolated and sporadic points on modular curves

These are the computations accompanying Filip Najman's paper *Isolated and sporadic points on X₁(N) and X₀(N)*. We consider noncuspidal closed points over the rational numbers and include CM points.

| Curve and question | Levels |
| --- | --- |
| An isolated point on X₁(N) | 17, 19 and every N ≥ 21 |
| An isolated point on X₀(N) | Every N except 1–10, 12, 13, 16, 18, 20, 24, 25, 32, 36 and 49 |
| A sporadic point on X₁(N) | Every N except 1–20, 22–27, 32 and the two open levels 83 and 131. The first 27 exceptions have no such point. |

The [level table](x1-sporadic/levels.csv) gives the noncuspidal degree bounds, the lower bounds for the minimum density degree, and their sources for 1 ≤ N ≤ 720. Clark, Genao, Pollack and Saia prove existence for every larger N. The field for a closed point and the proofs of isolation at the two open levels can be checked independently of the classification table.

## Reproducing the computations

Run these commands from the repository root. The table uses Python 3.9 or later and SymPy. The isolation computations use SageMath (tested with 10.7), Magma, a C++17 compiler and FLINT development headers accessible to `pkg-config`.

```sh
python3 x1-sporadic/levels.py
python3 x1-isolated/verify.py
sage -python x0-isolated/verify.py
sage -python x0-isolated/residual_degrees.py
```

The first command rebuilds the CSV. It gives 691 positive cases, 27 negative cases and the two open cases. The second checks the four cuspidal exhaustions, nine explicit small-level points, and the winding, irreducibility and cusp-jet calculations for the degree-143 point at level 131. Add `--rebuild-jets` to regenerate its expansion matrix or `--exact-points` to repeat the nine Riemann–Roch calculations in characteristic zero. The last two commands compute the cuspidal class counts, check the rank-one cases and recover the five residual-point polynomials on X₀(N).

The [sporadic-level instructions](x1-sporadic/README.md) give the individual tests behind the level table. For the two isolated diamond orbits of degree 57 at level 83, use the [level-83 instructions](level-83/README.md). Their sporadicity is open. The degree-143 point at level 131 is nonsporadic since the curve maps with degree 130 to a positive-rank elliptic curve.

## Sources

The finite inputs come from the [minimum CM degree computations](https://github.com/fsaia/least-cm-degree/tree/7eef6c9fd9d6d0bdea56eb7f339ffb1da7ee1c3c), the [point census](https://github.com/x1mn-sporadic/x1mn-sporadic.github.io/tree/dca0c196ae2bc3ffc09319b03171bba7323c1d16), [van Hoeij's modular-curve data](https://www.math.fsu.edu/~hoeij/files/X1N/), [Sutherland's models](https://math.mit.edu/~drew/X1_altcurves.html), and the [LMFDB](https://www.lmfdb.org/). The X₀(N) computation uses the classification in [Derickx and Najman](https://arxiv.org/abs/2511.09015v2).

OpenAI's GPT-6.2 Astra carried out the mathematical work on the new results, prepared the first draft and organized the computations in this repository.
