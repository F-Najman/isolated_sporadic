Van Hoeij searches on X_1(83)

The proof and search scope are in note/weighted83.tex in the archive.
From this directory run

    python3 replay.py /tmp/x183-check

This reconstructs the 435,204 fibres from the explicit modular-unit
catalogue, checks their cusp orders and all proper-factor exclusions,
and compares the normalized factor patterns with the recorded hashes.
It needs Python 3.11+ and NumPy. To recompute the elliptic norm
polynomials and the finite-field factorizations, run

    sage -python replay.py /tmp/x183-replay --factor --threads 8

Add --spaces to also regenerate the canonical cusp jets modulo 167
and verify the two- and three-cusp Riemann--Roch tests. The arithmetic
replay requires Sage with NumPy, SymPy and its Cremona database, a
C++17 compiler with OpenMP, FLINT 3, GMP and pkg-config. Keep assertions
enabled. The input exponents are fixed in catalogue.txt.gz. A repeat
of the heuristic selection can be made with

    sage -python weighted.py /tmp/x183-lll

The selection used 768 parameter choices, changing the lattice weights
at the cusps, and retained units of degree at most 175. The LLL output
can depend on the fpylll/fplll version. Use the archived catalogue for
the exact search reported in the note.

The 22 targets are all nonzero reduced a/b with |a| and b at most four
and b positive. The catalogue consists of 7,960 balanced units with at
most four factors on each side, 11,219 balanced units with five through
fourteen factors on each side, and 603 unbalanced units. For unbalanced
units the exponent sum is positive and the first norm at the nonsplit
cusp is independently known to be (-a)^41.

Of 433,002 positive remainders, 137 survived the first prime group and
none survived the second. Every tested proper point degree below 57
is excluded. The smallest positive remainder has degree 57. The
sibling siegel_residue83.py proves irreducibility of the known example.

For the Riemann--Roch test let P be the degree-57 point. All 820
pairs of rational cusps with multiplicity 32 and all 10,660 triples
with multiplicity 12 give l(P+A)=1. The source of the liftable cusp
jets is ../eisenjets.py. See the note for the deduction for smaller
multiplicities and for effective cuspidal subtractions.
