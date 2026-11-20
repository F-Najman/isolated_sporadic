"""Run with sage -python factors121.py. Verify the arithmetic inputs for level 121."""
from sage.all import (DirichletGroup, EllipticCurve, Gamma0, Gamma1, GammaH, Infinity,
                      ModularSymbols, PolynomialRing, PowerSeriesRing, QQ, euler_phi,
                      identity_matrix, matrix, prod)
from sage.databases.cremona import CremonaDatabase

level = 121
H = GammaH(level, [-1,12])
assert (Gamma1(level).genus(), H.genus(), Gamma0(level).genus()) == (526,26,6)
assert (2*526-2 - 11*(2*26-2))//10 == 50
assert max(120*r-50*r*(r-1)//2 for r in range(1,12)) == 210
assert (121**2-11**2)//24 == 605
assert EllipticCurve([1,-1,0,-2,-1]).cm_discriminant() == -7

# The even characters have orders 1,5,11,55. At each nontrivial order
# the rational cuspidal periods of T_n{0,infinity} span the entire
# plus space for one character; Galois conjugation covers the others.
generator = DirichletGroup(level).gen(0)
total_dimension = Gamma0(level).genus()
for order,dimension,bound in [(5,5,6),(11,10,10),(55,10,10)]:
    character = (generator**(110//order)).minimize_base_ring()
    assert character.order() == order and character.is_even()
    ambient = ModularSymbols(character,2,sign=1)
    cusp = ambient.cuspidal_submodule()
    assert cusp.dimension() == dimension
    period = cusp.rational_period_mapping()
    winding = ambient([0,Infinity])
    periods = [period(ambient.hecke_operator(n)(winding)) for n in range(1,bound+1)]
    assert matrix(periods).rank() == dimension
    total_dimension += euler_phi(order)*dimension
    print('Character order, cuspidal dimension and winding rank', order, dimension, dimension, flush=True)
assert total_dimension == 526

curves = list(CremonaDatabase().iter_optimal([11,121]))
elliptic = [(c.cremona_label(), int(c.ap(2)), int(c.rank())) for c in curves]
assert elliptic == [('11a1',-2,0),('121a1',-1,0),('121b1',0,1),('121c1',1,0),('121d1',2,0)]
print('Cremona labels, traces and ranks', elliptic, flush=True)

ring = PolynomialRing(QQ,'x')
x = ring.gen()
projector = prod(x**2-a*x+2 for _,a,rank in elliptic if not rank)
assert projector == x**8+3*x**6+8*x**4+12*x**2+16
assert projector % (x**2+2) == 16
assert list(projector % (x**7-1)) == [16,1,12,0,8,0,3]

# A count from the modular curve, independent of the enumeration of points
# on Fisher's equations. Eichler-Shimura gives the zeta numerator.
S = ModularSymbols(H,2,sign=1).cuspidal_submodule()
A = S.hecke_matrix(2)
diamond = S.diamond_bracket_matrix(2)
numerator = (identity_matrix(ring,26)-A*x+2*diamond*x**2).det()
series = PowerSeriesRing(QQ,'u',20)
u = series.gen()
zeta = series(numerator)/((1-u)*(1-2*u))
counts = [int(zeta[i]) for i in range(16,20)]
assert counts == [319800,639480,1278840,2557560]
print('Effective divisors in degrees 16 through 19', counts, flush=True)
