"""Verify the finite calculations for isolated points on X_0(N) over Q.

Run with `sage -python verify.py`.
"""
from math import prod

from sage.all import (BinaryQF_reduced_representatives, CuspForms, Gamma0,
                      Infinity, J0, Matrix, ModularSymbols, ZZ, divisors,
                      euler_phi, gcd, identity_matrix, lcm)


def cusp_coordinates(n):
    """Return degrees and group coordinates for closed cusps other than infinity."""
    symbols = ModularSymbols(n, 2, sign=0)
    periods = symbols.cuspidal_submodule().integral_period_mapping()
    cusps = symbols.cusps()
    weights, rows = [], []
    for d in divisors(n)[:-1]:
        orbit = [a for a in cusps if gcd(a.denominator(), n) == d]
        assert len(orbit) == euler_phi(gcd(d, n//d))
        weights.append(len(orbit))
        rows.append(sum(periods(symbols([Infinity, a])) for a in orbit))
    assert sum(weights)+1 == len(cusps)

    # The period lattice is Z^(2g). The relation lattice consists of the
    # integer row vectors whose product with the period matrix is integral.
    periods_matrix = Matrix(rows)
    scale = lcm([r.denominator() for r in periods_matrix.list()])
    numerators = (scale*periods_matrix).change_ring(ZZ)
    block = numerators.stack(scale*identity_matrix(ZZ, numerators.ncols()))
    relations = block.left_kernel_matrix()[:, :len(rows)]
    smith, left, right = relations.smith_form()
    assert left*relations*right == smith
    active = [i for i in range(len(rows)) if smith[i,i] > 1]
    moduli = tuple(int(smith[i,i]) for i in active)
    generators = [tuple(int(right[j,i] % smith[i,i]) for i in active) for j in range(len(rows))]
    return weights, moduli, generators


def effective_classes(g, weights, moduli, generators):
    """Group elements represented using closed cusps of total degree at most g."""
    zero = (0,)*len(moduli)
    by_degree = [set() for _ in range(g+1)]
    by_degree[0].add(zero)
    for weight, v in zip(weights, generators):
        for d in range(weight, g+1):
            for old in by_degree[d-weight]:
                by_degree[d].add(tuple((a+b) % m for a, b, m in zip(old, v, moduli)))
    return len(set.union(*by_degree))


# Theorems 1.1 and 1.2 of Derickx--Najman.
no_sporadic_cm = (set(range(1,11)) | {12,13} | set(range(15,19)) | set(range(20,27))
                  | set(range(28,34)) | set(range(35,38)) | set(range(39,42))
                  | set(range(46,51)) | {53,59,60,61,65,70,71,72,79,80,83,87,89,94,96,101,131,144})
no_sporadic = no_sporadic_cm - {15,17,21,37}
low_genus = {1,2,3,4,5,6,7,8,9,10,12,13,16,18,20,24,25,32,36,49}
bielliptic = {53,61,65,79,83,89,101,131}
remaining = no_sporadic - low_genus - bielliptic
assert {n for n in no_sporadic if Gamma0(n).genus() <= 1} == low_genus
assert len(no_sporadic) == 54 and len(remaining) == 26

print('N genus cuspidal_subgroup represented_at_degree_g', flush=True)
for n in sorted(remaining):
    g = int(Gamma0(n).genus())
    weights, moduli, generators = cusp_coordinates(n)
    order = prod(moduli)
    assert order <= J0(n).rational_cuspidal_subgroup().order()
    covered = effective_classes(g, weights, moduli, generators)
    print(n, g, order, covered, flush=True)
    assert covered < order

print('N genus q_p_coefficient_of_last_echelon_form', flush=True)
for p in [53,61,79,89,101]:
    space = CuspForms(Gamma0(p), 2)
    g = int(space.dimension())
    basis = space.echelon_basis()
    assert [f.q_expansion(g+1).valuation() for f in basis] == list(range(1,g+1))
    f = basis[-1].q_expansion(p+1)
    assert f[g] == 1 and f[p] != 0
    print(p, g, f[p], flush=True)

print('N genus degree_of_CM_ramification_point', flush=True)
for p, expected in [(83,3),(131,5)]:
    g = int(Gamma0(p).genus())
    h = len(BinaryQF_reduced_representatives(-p, primitive_only=True))
    print(p, g, h, flush=True)
    assert h == expected and h < g-1
print('All finite checks passed.')
