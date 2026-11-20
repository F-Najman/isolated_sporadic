"""Verify the upper bound on l((T-H+S)_+) at level 83."""
import numpy as np
from sage.all import GF, matrix

n = 83
exponents = {35: 1, 40: 1, 3: -1, 34: -1}
contact = {10: 8, 18: 27, 26: 3}
data = np.load('eis83.npz')
assert int(data['prime']) == 167
jets = matrix(GF(167), data['rows'])
selected = []
pole_degree = 0
for c, (orbit_size, precision, start) in enumerate(data['profile']):
    orbit_size, precision, start = map(int, (orbit_size, precision, start))
    # The units are balanced. At the unique nonrational closed cusp the order is zero.
    if c == 0:
        order = 0
        assert orbit_size == 41
    else:
        assert orbit_size == 1
        numerator = sum(e * (i*c % n) * ((i*c % n) - n) for i, e in exponents.items())
        assert numerator % (2*n) == 0
        order = numerator // (2*n)
    poles = max(0, -order)
    pole_degree += poles*orbit_size
    multiplicity = max(0, poles + 1 - contact.get(c, 0))
    assert multiplicity <= precision
    for j in range(orbit_size):
        first = start + j*precision
        selected.extend(range(first, first + multiplicity))
assert pole_degree == 95 and len(selected) == 174
rank = jets.matrix_from_rows(selected).rank()
assert rank == 173
print('Effective degree', len(selected), 'jet rank', rank, 'Riemann-Roch upper bound', 1+len(selected)-rank)
