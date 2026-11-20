# Sporadic points on X₁(N)

`python3 levels.py` rebuilds [levels.csv](levels.csv) from the frozen inputs in [data](data). The column `noncuspidal_degree_upper` bounds the degree of some point. A strict inequality with `density_lower` proves existence. Ten rows also use the cuspidal fibre argument. A witness bound in one of those rows need not be the degree of a known sporadic point.

The code in [import_inputs.py](import_inputs.py) reads checkouts of the minimum-CM-degree and point-census repositories identified in the [main README](../README.md). Pass their directories in that order. The form labels, exact Hecke coefficients, characters and elliptic modular degrees are already included. [lmfdb_queries.sql](lmfdb_queries.sql) lists the upstream fields. The table checks the genus against the total dimension, including oldform multiplicities, and uses rational arithmetic for the spectral comparisons.

## Individual calculations

Run the commands below from the [certificates](certificates) directory. The large modular-unit and cusp-jet tests at 90, 94, 96, 108, 120 and 150 have separate [commands and parameters](certificates/composite-commands.txt). These require Magma or SageMath with NumPy and, for 94, a C++17 compiler. Retain the Python and C++ assertions.

| Levels | Commands |
| --- | --- |
| 54, 59, 60, 71 | Run `magma -b ex54.m`, and analogously `ex59.m`, `ex60.m`, `ex71.m`. Then run `magma -b jets60.m`, `magma -b N:=59 primejets.m`, `magma -b N:=71 primejets.m`. |
| 81, 87, 95 | `magma -b verify81.m`, `magma -b verify87.m`, `magma -b verify95.m` |
| 99, 107 | `python3 siegel_contacts.py`, `python3 siegel_prime.py` |
| 117 | `sage -python degree117.py`, `magma -b differential117.m` |
| 121 | `python3 quotient121.py`, `sage -python factors121.py`. The first also needs NumPy and C++17. |
| 123 | `magma -b tate123.m` verifies the field and marking. `sage -python siegel_quotient123.py` reconstructs the fibre. |
| 138 | `sage -python nodes138.py`, `magma -b brandt138.m`, `sage -python factors138.py`, `sage -python fibre138.py` |
| 180 | `magma -b base36.m`, `magma -b point180.m` |
| 192 | `sage periods192.sage` |
| 210, 420 | `python3 levels210_420.py` |
| 231, 303 | `sage -python igusa231_303.py --verify-curves` checks the numerical inequalities, the elliptic curves, traces and CM markings. |

For the bound at 83 run `sage -python density83/series.py`. Other searches at 83 and 131 are retained with the corresponding scripts and finite input records. The [two degree-57 points](../level-83) and the [degree-143 point](../x1-isolated) have their own isolation checks. A search with a specified set of modular units only tests the indicated fibres.
