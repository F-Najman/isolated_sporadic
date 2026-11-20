"""Prepare certified Gram--Schmidt approximations for the level-180 lattice."""
from ast import literal_eval
from pathlib import Path
from sage.all import QQ, ZZ, diagonal_matrix, gcd, matrix

level, bound = 180, 117
valuation_path = Path(__file__).with_name('valuations180.txt')
valuations = matrix(ZZ, literal_eval(valuation_path.read_text()))
weights = valuations.right_kernel_matrix()[0]
weights *= 1 if min(weights) > 0 else -1
assert min(weights) == 1
weighted = valuations * diagonal_matrix(weights)
basis = weighted.BKZ(block_size=25)
assert basis.row_module() == weighted.row_module()
orthogonal, coefficients = basis.change_ring(QQ).gram_schmidt()
assert coefficients * orthogonal == basis
norms = [row.dot_product(row) for row in orthogonal]
rank, dimension = basis.dimensions()
scale, budget_scale, norm_scale = 2**15, 2**15, 2**12

units = [a for a in range(level) if gcd(a, level) == 1]
def fold(j):
    return min(j % level, (-j) % level)

for a in units:
    permutation = [fold(a * j) for j in range(dimension)]
    assert weighted.row_module() == weighted.matrix_from_columns(permutation).row_module()
pairs = []
seen = set()
for i in range(dimension):
    for j in range(i + 1, dimension):
        if (i, j) in seen:
            continue
        pairs.append((i, j))
        seen.update(tuple(sorted((fold(a*i), fold(a*j)))) for a in units)
assert len(seen) == dimension*(dimension-1)//2 and len(pairs) == 326

with open('sphere180.in', 'w') as stream:
    print(rank, dimension, bound, scale, budget_scale, norm_scale, len(pairs), file=stream)
    for row in basis:
        print(*row, file=stream)
    for row in coefficients:
        print(*[(x*scale).round() for x in row], file=stream)
    print(*[(r*norm_scale).floor() for r in norms], file=stream)
    for i, j in pairs:
        target = [((QQ(bound) * (g[i]-g[j])*scale) / (2*r)).round()
                  for g, r in zip(orthogonal, norms)]
        print(i, j, *target, file=stream)
