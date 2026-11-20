"""Recompute the degree and irreducibility of the single level-131 fibre."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / 'fibre'))
from vanhoeij131 import CUSPS, fibre, sieve

u = (19, 54, 27, 37)
auxiliary = (10, 65, 11, 55)
degree, residual, orders = fibre(u, 1)
assert (degree, residual) == (231, 143)
assert [(h, m) for h, m in zip(CUSPS, orders) if m > 0] == [(1, 19), (16, 39), (31, 29), (46, 1)]
_, _, d, possible, witnesses = sieve((u, 1, auxiliary, 24))
assert d == 143 and possible == []
for p, factors in witnesses:
    print('prime', p, 'factor degrees with multiplicity', factors)
print('The residual divisor is one closed point of degree 143.')
