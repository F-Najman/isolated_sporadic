"""Integral incidence bounds for new elliptic quotients of X_1(N)."""
from fractions import Fraction
from functools import cache
import json
from math import ceil, floor, isqrt
from pathlib import Path

DATA = Path(__file__).resolve().parent / 'data'
DEGREES = json.loads((DATA / 'elliptic_degrees.json').read_text())


@cache
def possible_optimal_degrees(label, diamond_degree):
    """Possible degrees of the optimal map from X_1(N) in this new isogeny class."""
    info = DEGREES[label]
    M = diamond_degree * info['optimal_X0_degree']
    if label == '192.2.a.a':
        # The two loops in certificates/periods192.sage span the optimal period lattice.
        return (M,)
    degrees = set()
    for curve in info['class']:
        torsion = curve['torsion']
        generator_degree = curve['isogeny_to_X0_optimal']
        for a in range(1, isqrt(torsion // generator_degree) + 1):
            e = a * a * generator_degree
            if torsion % e == 0 and diamond_degree % e == 0:
                degrees.add(M // e)
    assert M in degrees
    return tuple(sorted(degrees))


def elliptic_excluded(n, form, d, counts, phi_n, gonality):
    """Bound the integer k using rational gonality and all good primes."""
    label = form['label']
    if form['level'] != n or form['dim'] != 1 or label not in DEGREES:
        return False
    if all(S <= (p*p+1)*d for p,S in counts.items()):
        return False  # We need uniqueness on the characteristic-zero translate.
    # Bounds for u=m/k in the simultaneous-prime incidence inequalities.
    lower_u, upper_u = Fraction(0), Fraction(d)
    for p, S in counts.items():
        t = form['initial_traces'][p-1]**2 - 2*p
        slack = (p*p+1)*d - S
        if t > 0:
            upper_u = min(upper_u, slack/t)
        elif t < 0:
            lower_u = max(lower_u, slack/t)
        elif slack < 0:
            return True
    if lower_u > upper_u or upper_u <= 0:
        return True
    for m in possible_optimal_degrees(label, phi_n//2):
        if lower_u <= m <= upper_u:
            return False  # k=1 remains possible.
        # For k>=2, intersections with translates give gonality <= 2(dk-m).
        least_k = max(2, ceil(m/upper_u), ceil(Fraction(2*m+gonality,2*d)))
        if lower_u == 0 or least_k <= floor(m/lower_u):
            return False
    return True


def old_equality_excluded(n, form, d, counts, diamond_order):
    """At a scalar elliptic endpoint the family is an optimal degree-d map."""
    label = form['label']
    if not form['level'] < n < 778 or form['dim'] != 1 or label not in DEGREES:
        return False
    traces = form['initial_traces']
    if not any(S > (p*p+1)*d and S == ((p+1)**2-traces[p-1]**2)*d for p,S in counts.items()):
        return False
    data = DEGREES[label]
    # After quotienting by the translations of the diamonds, the map
    # descends to X0(n). Derickx--Najman, Lemma 3.3, gives the divisibility.
    modulus = diamond_order * data['optimal_X0_degree']
    for curve in data['class']:
        for e in range(1, curve['torsion']+1):
            if curve['torsion']%e == 0 and diamond_order%e == 0 and (d*e)%modulus == 0:
                return False
    return True


def level_255_endpoint():
    """The complete table of rational elliptic isogeny classes at divisors of 255."""
    data = json.loads((DATA / 'elliptic_255.json').read_text())
    assert all(255 % conductor == 0 and rank == 0 for label, conductor, rank in data)
    assert {row[0] for row in data} == {'15.a1', '17.a1', '51.a1', '85.a1'}
    return 257


@cache
def _hecke_orders():
    return {r['label']:r for r in json.loads((DATA / 'hecke_spectral.json').read_text())}


def scalar_equality_excluded(form, d, counts, beta):
    """Exclude exact negative-scalar endpoints for higher-dimensional simple factors."""
    if form['dim'] == 1:
        return False
    candidates = [p for p, S in counts.items()
                  if beta[p] < 2*p and S == d*((p+1)**2-beta[p])]
    if not candidates:
        return False
    from hecke import values
    extrema = values(form, _hecke_orders())
    # A rational lower bound for the smallest conjugate that equals the
    # exact trace mean forces every conjugate to equal that mean.
    return any(extrema[p][0] == extrema[p][1] == beta[p] for p in candidates)


def combined_elliptic_excluded(form, d, counts):
    """Compare the same mixed-class norm to the incidence degree at all primes."""
    if form['dim'] != 1:
        return False
    lower, upper = Fraction(0), Fraction(d)
    for p, S in counts.items():
        t = form['initial_traces'][p - 1] ** 2 - 2*p
        gap = (p*p + 1)*d - S
        if t < 0:
            lower = max(lower, gap/t)
        elif t > 0:
            upper = min(upper, gap/t)
        elif gap < 0:
            return True
    return lower > upper or upper <= 0
