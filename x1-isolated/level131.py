"""Evaluate canonical differentials at B=(T-H+S)_+ on X_1(131)."""
import sys

import numpy as np
from sage.all import GF, matrix


data = np.load(sys.argv[1])
assert int(data['prime']) == 263
forms = matrix(GF(263), data['rows'])
exponents = {19: 1, 54: 1, 27: -1, 37: -1}
contacts = {}
pole_degree = 0
rows = []

for h, profile in enumerate(data['profile']):
    orbit_size, precision, start = map(int, profile)
    if h == 0:
        # The nonsplit orbit has neither a pole nor contact with u=1.
        assert orbit_size == 65
        multiplicity = 1
    else:
        assert orbit_size == 1
        c = h if h % 2 else h + 131
        residues = {i: c*i % 131 for i in exponents}
        numerator = sum(exponents[i]*r*(r-131) for i, r in residues.items())
        order, remainder = divmod(numerator, 262)
        assert remainder == 0
        parity = sum(e*(c*i//131) for i, e in exponents.items()) % 2
        if order == 0 and parity == 0:
            contacts[h] = min(min(r, 131-r) for r in residues.values())
        poles = max(0, -order)
        pole_degree += poles
        multiplicity = max(0, 1 + poles - contacts.get(h, 0))
    assert multiplicity <= precision
    for j in range(orbit_size):
        first = start + j*precision
        rows.extend(range(first, first + multiplicity))

rank = int(forms.matrix_from_rows(rows).rank())
assert pole_degree == 231 and contacts == {1: 19, 16: 39, 31: 29, 46: 1}
assert len(rows) == 357 and rank == 356
print('pole degree', pole_degree, 'contacts', contacts, 'effective degree', len(rows), 'rank', rank)
