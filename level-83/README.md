# Two degree-57 points on X₁(83)

The two points are the residual divisors of the fibres at 1 of E₃₅E₄₀/(E₃E₃₄) and at −1 of E₂E₉/(E₁₈E₃₃). The [polynomials](records/minima-all.json) are in ascending coefficient order. They are irreducible and have different squarefree factorizations modulo 2, so their residue fields are different.

From this directory run

```sh
sage -python replay.py /tmp/x183-isolated
```

This rebuilds the rational coefficient strata through degree 113, the coordinate polynomials through degree 89, the jet check for the second point and the isogeny-class witnesses. The first point also has a [residue computation](../x1-sporadic/certificates/siegel_residue83.py) and a [tangent-space check](../x1-sporadic/certificates/tangent83.py).

Add `--rebuild-jets` to regenerate the cusp expansions, `--hecke` to verify modular correspondences through 31, or `--full --threads 8` to repeat the larger fibre searches. The larger runs use FLINT 3, GMP, OpenMP and a C++17 compiler. The [degree-52 calculations](degree-52) concern necessary conditions for a moving family. Whether there is any sporadic point at level 83 remains open.
