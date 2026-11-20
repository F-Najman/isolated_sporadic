"""Compute T_5 directly on cyclic order-69 subgroups of the supersingular curve."""
from collections import Counter
from itertools import product
from sympy import Matrix, eye


def act(quaternion, line, p):
    # A line is spanned by (line,1), with p denoting the span of (1,0).
    # The matrices of i and j are ((0,-1),(1,0)) and ((a,b),(b,-a)).
    a, b = {3: (1, 1), 5: (0, 2), 23: (2, 8)}[p]
    assert (a*a + b*b + 1) % p == 0
    r, s, t, u = quaternion
    x, y = (1, 0) if line == p else (line, 1)
    X = ((r + t*a - u*b)*x + (-s + t*b + u*a)*y) % p
    Y = ((s + t*b + u*a)*x + (r - t*a + u*b)*y) % p
    if X == Y == 0:
        return None
    return X * pow(Y, -1, p) % p if Y else p


def pair_action(q, pair):
    return (act(q, pair[0], 3), act(q, pair[1], 23))


# Projectively a quaternion and any nonzero scalar multiple agree.
# The first four units represent 1,i,j,k. The other eight represent
# (1 +/- i +/- j +/- k)/2 after discarding the scalar denominator.
units = [tuple(int(k == j) for k in range(4)) for j in range(4)]
units += [(1,) + signs for signs in product((-1, 1), repeat=3)]
orbits = {}
representatives = []
for pair in product(range(4), range(24)):
    if pair in orbits:
        continue
    orbit = {pair_action(q, pair) for q in units}
    assert len(orbit) == 12 and not orbit.intersection(orbits)
    for image in orbit:
        orbits[image] = len(representatives)
    representatives.append(pair)
assert len(representatives) == 8

# Six endomorphisms of norm 5, one for each order-5 kernel on A.
isogenies = []
for j in (1, 2, 3):
    for sign in (-1, 1):
        q = [1, 0, 0, 0]
        q[j] = 2*sign
        isogenies.append(q)
kernels = [{L for L in range(6) if act(q, L, 5) is None} for q in isogenies]
assert all(len(kernel) == 1 for kernel in kernels) and set.union(*kernels) == set(range(6))

def hecke_row(pair):
    return Counter(orbits[pair_action(q, pair)] for q in isogenies)

rows = [hecke_row(pair) for pair in representatives]
assert all(hecke_row(pair) == rows[i] for pair, i in orbits.items())
T = Matrix([[row[i] for i in range(8)] for row in rows])
assert T == T.T
line = (T + 2*eye(8)).nullspace()
assert len(line) == 1
v = list(line[0])
scale = next(x for x in v if x)
v = [x/scale for x in v]
assert sorted(v) == [-1,-1,0,0,0,0,1,1]
print('Orbit representatives in P^1(F_3) x P^1(F_23)', representatives)
print('T_5 on supersingular classes')
for row in T.tolist():
    print(*row)
print('Eigenvector for a_5=-2', v)
print('Multiplicity of lifted coordinates', {c: 22*v.count(c) for c in (-1,0,1)})
