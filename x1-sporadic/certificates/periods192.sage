# sage periods192.sage
# The period lattice of the 192.a rational newform on Gamma_1(192).
from sage.all import *
N = 192
M = ModularSymbols(N, 2)
ap = {2: 0, 3: -1, 5: -2, 7: -4, 11: -4, 13: 2, 17: -6, 19: 4}
spaces = [A for A in M.cuspidal_submodule().new_subspace().decomposition()
          if A.dimension() == 2 and all(A.hecke_matrix(p) == a*identity_matrix(2)
                                       for p, a in ap.items())]
assert len(spaces) == 1
period = spaces[0].integral_period_mapping()
lattice_vectors = []
for a, c, d in [(1153, 1344, 577), (385, 1728, 193)]:
    assert (a*d-1) % c == 0 and c % N == 0 and (a-1) % N == (d-1) % N == 0
    # g = [a, (ad-1)/c; c, d] belongs to Gamma_1(N).
    # This loop joins infinity to g(infinity) in the upper half-plane.
    lattice_vectors.append(period(M([Infinity, QQ(a)/c])))
print(matrix(lattice_vectors))
assert abs(matrix(lattice_vectors).det()) == 1
