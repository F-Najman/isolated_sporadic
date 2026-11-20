"""Verify every cuspidal divisor S+E with bounded excess. Bucket budgets give a finite covering by larger divisors. For dependent supersets, the relation-space recursion rules out every relation whose support fits into an allowed divisor. Run under Sage Python with N, excess, step size, and bucket size as arguments."""

import numpy as np, random, sys, itertools, math
from sage.all import matrix, GF, identity_matrix
N, k, s, size = map(int, sys.argv[1:5])
data = np.load(f'eis{N}.npz')
G = GF(int(data['prime']))
R = matrix(G, data['rows'])
O = [(int(d), [list(range(int(start) + j * int(m), int(start) + (j + 1) * int(m))) for j in range(int(d))]) for d, m, start in data['profile']]
seed = 20261118
random.Random(seed).shuffle(O)
F = len(O)
parts = [O[i:i + size] for i in range(0, F, size)]

def inds(o, u):
    d, lists = o
    m = (d == 1) + u // d
    return [i for L in lists for i in L[:m]]
base = [j for o in O for j in inds(o, s - 1)]
B = R.matrix_from_rows(base)
print('base', B.rank(), len(base), flush=True)
assert B.rank() == len(base)
T = R * B.right_kernel_matrix().transpose()
Bset = set(base)
ext = B.stack(identity_matrix(GF(int(data['prime'])), R.ncols()))
framing = ext.matrix_from_rows(ext.transpose().pivots())
Z = R * framing.inverse()
Z0 = Z.matrix_from_columns(list(range(len(base))))
lookup = {v: (i, degree, ell + 1 - (degree == 1)) for i, (degree, pts) in enumerate(O) for seq in pts for ell, v in enumerate(seq)}

def loww(z, ids):
    choices = {}
    for j, y in enumerate(z):
        if y:
            i, degree, extra = lookup[ids[j]]
            weight = degree * extra
            if i not in choices or weight > choices[i][0]:
                choices[i] = (weight, j)
    return sorted(choices.values(), reverse=True)

def exclude(K, ids):
    h = K.nrows()
    if not h:
        return True
    chosen = []
    total = 0
    for x, j in loww(K.row(0), ids):
        chosen.append(j)
        total += x
        if total > k:
            break
    if total <= k:
        return False
    if h == 1:
        return True
    gone = set()
    for j in chosen:
        col = K.column(j)
        position = next((i for i, y in enumerate(col) if y))
        key = tuple(col / col[position])
        if key in gone:
            continue
        gone.add(key)
        sub = matrix(G, [col]).right_kernel_matrix() * K
        if not exclude(sub, ids):
            return False
    return True
levels = [[tuple((v for o in bucket for v in inds(o, min(k, s * (i + 1) - 1)) if v not in Bset)) for i in range(k // s + 1)] for bucket in parts]
print('start', len(parts), math.comb(k // s + len(parts) - 1, len(parts) - 1), flush=True)
bad = []
good = 0
for alloc in itertools.combinations_with_replacement(range(len(parts)), k // s):
    if len(bad) >= 2:
        break
    cnt = [alloc.count(j) for j in range(len(parts))]
    U = sum((levels[j][c] for j, c in enumerate(cnt) if c), ())
    M = T.matrix_from_rows(U)
    if M.rank() != len(U):
        K = M.left_kernel_matrix()
        K = K.augment(-K * Z0.matrix_from_rows(U))
        if not exclude(K, list(U) + base):
            bad.append(cnt)
            print('possible dependency', cnt, flush=True)
            continue
    good += 1
    if good % 10000 == 0:
        print(good, flush=True)
print('DONE', N, good, len(bad), bad, flush=True)
assert not bad and good == math.comb(k // s + len(parts) - 1, len(parts) - 1)
