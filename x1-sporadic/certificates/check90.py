"""Verify every enlarged jet system that contains a candidate cuspidal divisor S+E at level 90. Requires eis90.npz in the working directory."""

import numpy as np, itertools
from sage.all import matrix, GF
s = np.load('eis90.npz')
k = GF(int(s['prime']))
R = matrix(k, s['rows'])
O = [(int(d), [list(range(int(start) + j * int(m), int(start) + (j + 1) * int(m))) for j in range(int(d))]) for d, m, start in s['profile']]
lo = {1: 3, 2: 2}
base = [r for d, q in O for p in q for r in p[:lo.get(d, 0)]]
B = R.matrix_from_rows(base)
print('base', B.rank(), len(base), flush=True)
assert B.rank() == len(base)
T = R * B.right_kernel_matrix().transpose()
hi1 = [[r for p in q for r in p[3:]] for d, q in O if d == 1]
hi2 = [[r for p in q for r in p[2:]] for d, q in O if d == 2]
other = [[r for p in q for r in p] for d, q in O if d in (4, 6)]
cases = 0
for a, b in itertools.combinations(hi1, 2):
    for c in hi2:
        for e in other:
            ix = a + b + c + e
            assert T.matrix_from_rows(ix).rank() == len(ix)
            cases += 1
print('PASS 90', cases)
