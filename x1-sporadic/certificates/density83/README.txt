X_1(83): the density bound and the signed modular-unit fibres

From the directory containing this file, run

    sage -python series.py

This checks the norm of j modulo 2 by its Laurent expansion at the only
pole, the multiplicities of its five elliptic-curve values, and the
orders of E(F_2) and E(F_4). The bound b >= 52 follows in the note.

To repeat the complete modular-unit enumeration, choose a writable
working directory with about 1 GB free. These commands require SageMath
with its Cremona curve database, NumPy, SymPy, g++, OpenMP, GMP, FLINT
3 and pkg-config. Keep assertions enabled.

    sage -python series.py --units ./run83
    sage -python cases.py ./run83
    python3 sieve.py ./run83 --replay --jobs 4 --threads 4

The jobs and threads arguments control the number of concurrent chunks
and the threads per norm calculation. Reducing either does not change
the arithmetic or the recorded hashes. The manifest records a SHA-256
digest of the normalized prime and factor degrees for each input block
and stage. The replay verifies those digests, all multiplicity totals,
the recorded subset-sum exclusions and the final list of survivors.
It produces the following counts.

      factors    positive fibres     after primes 5--29     after 31--61
         2               181                    0                   0
         3            35,769                    2                   0
         4         3,750,521                  491                   0

The least noncuspidal degree in these signed fibres is 57. For a
separate proof of attainment, run the sibling certificate

    sage -python ../siegel_residue83.py

It computes the characteristic polynomial of a second modular unit on
the degree-57 remainder by exact Newton sums and proves it irreducible
modulo 11. The paper describes the isolation of this point and the
remaining density degrees 52 through 57 relevant to its sporadicity.
