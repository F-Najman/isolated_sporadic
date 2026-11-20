"""Run with sage -python factors138.py. Exact winding spaces at level 138."""
from sage.all import DirichletGroup, Gamma0, Gamma1, Infinity, ModularSymbols, euler_phi, matrix
from sage.databases.cremona import CremonaDatabase

level = 138
assert (Gamma1(level).genus(), Gamma1(69).genus(), Gamma0(level).genus()) == (441, 133, 21)
characters = {}
for character in DirichletGroup(level):
    if character.is_even():
        characters.setdefault(character.order(), character.minimize_base_ring())
assert sorted(characters) == [1, 2, 11, 22]

total = 0
for order, character in sorted(characters.items()):
    symbols = ModularSymbols(character, 2, sign=1)
    cusp = symbols.cuspidal_submodule()
    dimension = cusp.dimension()
    target = dimension - int(order == 1)
    period = cusp.rational_period_mapping()
    winding = symbols([0, Infinity])
    vectors = []
    for n in range(1, 49):
        vectors.append(period(symbols.hecke_operator(n)(winding)))
        if n >= target and matrix(vectors).rank() == target:
            break
    else:
        raise ArithmeticError(f'Insufficient winding rank at character order {order}')
    total += euler_phi(order) * dimension
    print('order, dimension, winding rank, last operator', order, dimension, target, n, flush=True)
assert total == 441

curves = list(CremonaDatabase().iter_optimal([level]))
matches = [curve for curve in curves if curve.rank(proof=True) > 0]
assert len(matches) == 1
curve = matches[0]
assert curve.conductor() == 138 and curve.rank(proof=True) == 1
assert curve.root_number() == -1 and curve.ap(5) == -2
assert curve.has_multiplicative_reduction(2)
print('positive-rank curve, model, rank, a_5', curve.cremona_label(), curve.a_invariants(), 1, -2)
