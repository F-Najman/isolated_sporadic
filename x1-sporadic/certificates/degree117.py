"""Run with sage -python degree117.py inside the unpacked computation archive."""
import json
import sys
from fractions import Fraction
from pathlib import Path

from sage.all import Gamma1, Infinity, ModularSymbols, QQ, ZZ, divisors, matrix

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
from hecke import values

raw = json.loads((root / 'data/form_metadata.json').read_text())
forms = [f for f in raw if 117 % f['level'] == 0]
assert sum(f['dim']*len(divisors(117//f['level'])) for f in forms) == Gamma1(117).genus() == 409
elliptic = {f['label'] for f in forms if f['dim'] == 1}
assert elliptic == {'39.2.a.a','117.2.a.a'}
possible = {f['label']:f for f in forms if not (f['analytic_rank_proved'] and f['analytic_rank'] == 0)}
checks = {'117.2.e.a':(2,4), '117.2.r.a':(2,3), '117.2.t.a':(5,12), '117.2.z.a':(2,6)}
assert set(possible) == set(checks) | {'117.2.a.a'}

J2 = 117**2 * (1-Fraction(1,9)) * (1-Fraction(1,169))
S = lambda p: (p-1)*J2/24
assert S(2) == 504 and 100 < S(2)/5 <= 101
raw_hecke = {r['label']:r for r in json.loads((root / 'data/hecke_spectral.json').read_text())}
for label,(p,b) in checks.items():
    assert values(possible[label], raw_hecke)[p][2] == b
    bound = S(p)/max(p*p+1, (p+1)**2-b)
    assert bound > 72
    print(label, p, b, bound)

# Check every possible elliptic quotient independently in the trivial-character
# new modular-symbol spaces. An elliptic curve over Q has trivial nebentypus.
classes = {}
for level in divisors(117):
    M = ModularSymbols(level,2)
    for A in M.cuspidal_submodule().new_submodule().decomposition():
        if A.dimension() == 2:
            assert level not in classes
            classes[level] = (M,A)
assert set(classes) == {39,117}
assert classes[39][1].hecke_matrix(2) == 1
M,A = classes[117]
assert A.hecke_matrix(2) == -1 and A.abelian_variety().modular_degree() == 8
period = A.integral_period_mapping()
vectors = []
for a,c,d in [(469,585,469), (469,1287,118)]:
    assert (a*d-1) % c == 0 and c % 117 == 0 and (a-1) % 117 == (d-1) % 117 == 0
    vectors.append(period(M([Infinity,QQ(a)/c])))
assert matrix(vectors) == matrix(ZZ, [[2,-1],[3,-1]])
assert abs(matrix(vectors).det()) == 1
m = 36*8
for d in (71,72):
    ks = [k for k in range(1,20) if S(2)*k <= 5*d*k+3*m
          and (m <= d if k == 1 else 101 <= 2*(d*k-m))]
    assert ks == ([5] if d == 71 else [5,6])
    assert 5*d-m < 101 and 5*d-m < m
assert (1-(-1)) % 5 != 0  # the conductor-39 and conductor-117 5-torsion differ
print('X1(117): m = 288, geometric gonality >= 101, density degree >= 72.')
