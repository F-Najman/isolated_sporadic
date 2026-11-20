"""Count the new-character part of J_Y(F_2) by Eichler–Shimura."""
from sage.all import DirichletGroup, ModularSymbols, ZZ, identity_matrix

chi = (DirichletGroup(131).gen()**10).minimize_base_ring()
space = ModularSymbols(chi, 2, sign=1).cuspidal_submodule()
assert chi.order() == 13 and space.dimension() == 10
fac = (1 + 2*chi(2))*identity_matrix(chi.base_ring(), 10) - space.hecke_matrix(2)
order = ZZ(fac.det().norm())
assert order == 773727027031129171413080804164258787077
assert order % 13 == 1 and order % 3 != 0
print("Prym order at two", order)
