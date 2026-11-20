"""Compute the winding ranks for all character orbits at level 131."""
from sage.all import DirichletGroup, EllipticCurve, Infinity, ModularSymbols, matrix


generator = DirichletGroup(131).gen()
for order, dimension, target, last in [(1, 11, 10, 24), (5, 10, 10, 11), (13, 10, 10, 10), (65, 10, 10, 10)]:
    character = (generator**(130//order)).minimize_base_ring()
    ambient = ModularSymbols(character, 2, sign=1)
    cusp = ambient.cuspidal_submodule()
    assert character.is_even() and character.order() == order and cusp.dimension() == dimension
    period = cusp.rational_period_mapping()
    winding = ambient([0, Infinity])
    orbit = [period(ambient.hecke_operator(i)(winding)) for i in range(1, last+1)]
    rank = matrix(orbit).rank()
    print('character order', order, 'dimension', dimension, 'winding rank', rank, flush=True)
    assert rank == target

E = EllipticCurve([0, -1, 1, 1, 0])
assert E.conductor() == 131 and E.rank(proof=True) == 1
print('The elliptic quotient at level 131 has rank one.')
