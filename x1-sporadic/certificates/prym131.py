"""Rational two-torsion and winding rank of the degree-13 diamond part."""
from sage.all import DirichletGroup, Infinity, ModularSymbols, ZZ, identity_matrix, matrix

chi = (DirichletGroup(131).gen()**10).minimize_base_ring()
assert chi.order() == 13 and chi.is_even()
ambient = ModularSymbols(chi, 2, sign=1)
space = ambient.cuspidal_submodule()
assert space.dimension() == 10
at_three = (1+3*chi(3))*identity_matrix(chi.base_ring(),10) - space.hecke_matrix(3)
order = ZZ(at_three.determinant().norm())
assert order > 0 and order % 2 == 1
print("Prym order at 3", order, flush=True)
period = space.rational_period_mapping()
winding = ambient([0, Infinity])
vectors = []
for n in range(1,35):
    vectors.append(period(ambient.hecke_operator(n)(winding)))
    if n >= 10 and matrix(vectors).rank() == 10:
        break
assert matrix(vectors).rank() == 10
print("winding rank", len(vectors), matrix(vectors).rank(), flush=True)
