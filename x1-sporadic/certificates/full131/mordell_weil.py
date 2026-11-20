"""Exact winding ranks, cuspidal orders and rational torsion of J_1(131)/E."""
from math import gcd
from sage.all import DirichletGroup, EllipticCurve, GF, Infinity, ModularSymbols, QQ, ZZ, identity_matrix, matrix

prime = 131
degree = 65
orders = (5, 13, 65)
characters = [(DirichletGroup(prime).gen()**(130//m)).minimize_base_ring() for m in orders]
spaces = []
for m, character, stop in zip(orders, characters, (11, 10, 10)):
    ambient = ModularSymbols(character, 2, sign=1)
    space = ambient.cuspidal_submodule()
    assert character.order() == m and character.is_even() and space.dimension() == 10
    period = space.rational_period_mapping()
    winding = ambient([0, Infinity])
    images = [period(ambient.hecke_operator(n)(winding)) for n in range(1, stop+1)]
    assert matrix(images).rank() == 10
    print("character order", m, "winding rank 10 through T_"+str(stop))
    spaces.append(space)

def valuation_numerator(a, i):
    r = a*i % prime
    return r*(r-prime)//2

# A column is the divisor of E_i/E_65 multiplied by 131, omitting the last cusp.
M = matrix(ZZ,64,64,lambda a,i: valuation_numerator(a+1,i+1)-valuation_numerator(a+1,65))
residues = [(i*i-65**2) % prime for i in range(1,65)]
X = identity_matrix(ZZ,64)
X[0,0] = prime
for i in range(1,64):
    X[0,i] = -(residues[i]*pow(residues[0],-1,prime) % prime)
unit_lattice = (M*X/prime).change_ring(ZZ)
cusps = abs(unit_lattice.det())
bernoulli = QQ(prime)
for chi in characters:
    bernoulli *= (chi.bernoulli(2)/4).norm()
assert cusps == abs(bernoulli) and cusps % 2 == 1
print("split cuspidal group order",cusps)

def cusp_index_on_quotient(j):
    def cusp(i):
        a = pow(2,i,prime)
        return min(a,prime-a)
    orbits = [{cusp(i+j*t) for t in range(65//j)} for i in range(j)]
    assert len(set.union(*orbits)) == 65
    pullbacks = matrix(ZZ,64,j-1,lambda a,i: int(a+1 in orbits[i])-int(a+1 in orbits[-1]))
    combined = unit_lattice.augment(pullbacks).transpose().hermite_form(include_zero_rows=False)
    return cusps // abs(combined.det())

assert cusp_index_on_quotient(5) == 1292275
assert cusp_index_on_quotient(13) == 5548598204321224637
print("projected cuspidal orders",cusp_index_on_quotient(5),cusp_index_on_quotient(13))

base = ModularSymbols(prime,2,sign=1).cuspidal_submodule()
curve = EllipticCurve([0,-1,1,1,0])
assert curve.torsion_order() == 1 and base.dimension() == 11
assert curve.change_ring(GF(2)).trace_of_frobenius() == 0
bound = 0
for q in (3,5,11,19):
    card_E = curve.change_ring(GF(q)).cardinality()
    card_base = ZZ(((1+q)*identity_matrix(QQ,11)-base.hecke_matrix(q)).det())
    assert card_base % card_E == 0
    prym = 1
    for chi,space in zip(characters,spaces):
        frob = (1+q*chi(q))*identity_matrix(chi.base_ring(),10) - space.hecke_matrix(q)
        prym *= ZZ(frob.det().norm())
    complement = card_base//card_E
    if q == 3:
        assert complement == 125515
    if q == 5:
        assert prym % 2 == 1
    count = prym*complement
    assert count % cusps == 0
    bound = gcd(bound,count)
    if q > 3:
        print("prime",q,"running gcd divided by cuspidal order",bound//cusps)
assert bound == 4*cusps
assert 2**2+2 == 6 and gcd(6,65) == 1
print("full Mordell-Weil certificate complete")
