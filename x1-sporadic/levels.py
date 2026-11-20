"""Rebuild the level table from the source data in this directory."""
import csv
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from pathlib import Path
from certificates.levels210_420 import small_prime_gonality
from certificates.igusa231_303 import density_bound as igusa_elliptic_bound
from joint import joint_excluded
from top import top_bound
from hecke import eigenvalue_lower_bounds
from elliptic import elliptic_excluded, old_equality_excluded, scalar_equality_excluded, level_255_endpoint

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
inputs = json.loads((DATA / 'level_inputs.json').read_text())
forms = json.loads((DATA / 'form_metadata.json').read_text())
cm_values = [0] + inputs['cm_minimum']
no_cm = set(inputs['cm_no'])
unknown_cm = set(inputs['cm_undecided'])
assert len(no_cm) == 67 and len(unknown_cm) == 227 and no_cm.isdisjoint(unknown_cm)
rank_zero = (set(range(1, 37)) | set(range(38, 43)) | set(range(44, 53)) |
             {54,55,56,59,60,62,64,66,68,69,70,71,72,75,76,78,81,84,87,90,94,
              96,98,100,108,110,119,120,132,140,150,168,180})
negative = set(range(1, 21)) | set(range(22, 28)) | {32}
unit_cases = {54,59,60,71,90,94,96,108,120,150}
primes = (2,3,5,7,11,13,17,19)
form_levels = defaultdict(list)
for f in forms:
    form_levels[f['level']].append(f)

@cache
def divisors(n):
    return [m for m in range(1, n + 1) if n % m == 0]

@cache
def phi(n):
    return sum(math.gcd(a, n) == 1 for a in range(1, n + 1))

@cache
def jordan(n):
    value = n * n
    for p in divisors(n):
        if len(divisors(p)) == 2:
            value = value // (p*p) * (p*p-1)
    return value

def index(n):
    return {1:1, 2:3}.get(n, jordan(n)//2)

def occurs_in_positive_rank(f):
    return not (f['analytic_rank_proved'] and f['analytic_rank'] == 0)

# The file contains the raw Hecke data. Isolate real roots over Q each run.
beta = eigenvalue_lower_bounds(forms, DATA)

def spectral_bound(n, good, counts, gonality):
    forms_n = [f for m in divisors(n) for f in form_levels[m]]
    dimension = sum(f['dim'] * len(divisors(n//f['level'])) for f in forms_n)
    genus = 1 + Fraction(jordan(n),24) - Fraction(sum(phi(m)*phi(n//m) for m in divisors(n)),4)
    if dimension != genus:
        return None  # The downloaded subset does not cover this level.
    possible = [f for f in forms_n if occurs_in_positive_rank(f)]
    if not possible:
        return None  # Handled by the published rank-zero list.
    bounds = {f['label']:max((max(counts[p] / max(p*p+1, (p+1)**2-beta[f['label']][p]), top_bound(f, p, counts[p])),
                              p, f['label']) for p in good) for f in possible}
    basic, p, label = min(bounds.values())
    lo = math.ceil(basic)
    hi = gonality
    while lo < hi:
        d = (lo + hi) // 2
        if all(bounds[f['label']][0] > d or elliptic_excluded(n, f, d, counts, phi(n), gonality)
               or old_equality_excluded(n, f, d, counts, phi(n)//2)
               or scalar_equality_excluded(f, d, counts, beta[f['label']])
               or joint_excluded(f, d, counts) for f in possible):
            lo = d + 1
        else:
            hi = d
    reason = 'integral or equality refinement on newform factors' if lo > math.ceil(basic) else f'newform comparison (limiting form {label}, p={p})'
    return lo, reason

def build_table():
    limits = range(1, len(cm_values))
    direct = {n:(cm_values[n], 'CGPS minimum CM degree') for n in limits}
    assert direct[121][0] == 110
    direct[121] = (110, 'CM discriminant -7, generator of the cyclic kernel of (2+sqrt(-7))^2')
    for n, degree in ((210, 96), (420, 384)):
        assert cm_values[n] == degree
        direct[n] = (degree, 'CM discriminant -24, marking (sqrt(-6)-48)/N, certificates/levels210_420.py')
    for n, degree in ((231, 240), (303, 400)):
        assert cm_values[n] == degree
        direct[n] = (degree, 'CM discriminant -24, marking (sqrt(-6)-216)/N, certificates/igusa231_303.py')
    for ns, candidate in inputs['noncuspidal'].items():
        n = int(ns)
        if candidate[0] < direct[n][0]:
            direct[n] = tuple(candidate)
    for ns, candidate in json.loads((DATA / 'extra_points.json').read_text()).items():
        n = int(ns)
        if candidate[0] < direct[n][0]:
            direct[n] = tuple(candidate)
    # Project to gcd, then lift. A fiber of a finite map has a point of
    # degree at most the degree of the map times the degree downstairs.
    witnesses = {}
    for n in limits:
        degree, origin = direct[n]
        for m,(old_degree,old_origin) in direct.items():
            h = math.gcd(n,m)
            candidate = old_degree * (index(n)//index(h))
            if candidate < degree:
                degree = candidate
                origin = f'{old_origin} (level {m} to {h} to {n})'
        witnesses[n] = (degree,origin)
    bounds = {}
    rows = []
    for n in limits:
        g = inputs['census_gonality_low'].get(str(n),1)
        reason = 'gonality or Frey bound'
        density = math.ceil(Fraction(g, 1 if n in rank_zero else 2))
        if n >= 5:
            good = [p for p in primes if n % p]
            counts = {p:Fraction((p-1)*jordan(n),24) for p in good}
            g = max(g, math.ceil(Fraction(975*index(n),98304)),
                    max(math.ceil(counts[p]/(p*p+1)) for p in good))
            if n in (138, 210, 420):
                red_prime = 3 if n == 420 else 2
                nodes = (red_prime-1)*jordan(n//red_prime)//24
                g = max(g, small_prime_gonality(red_prime, nodes))
            density = math.ceil(Fraction(g, 1 if n in rank_zero else 2))
            reason = 'rank zero and gonality' if n in rank_zero else 'Frey and gonality'
            if n in (210, 420):
                reason = f'Frey and two-component reduction at {red_prime}'
            if n == 180:
                nodes = 4*jordan(36)//24
                assert nodes == 144 and n in rank_zero
                # A noninvariant function has degree at least 3*nodes/4.
                # An invariant function descends to X_1(36) of F_5-gonality 8.
                # The constant and the two-invariant cases cost at least nodes.
                igusa_bound = min(nodes, 3*nodes//4 + 4*8)
                assert igusa_bound == 140
                g = max(g, igusa_bound)
                density, reason = g, 'rank zero and Igusa components at 5'
            basic,p = max((counts[p]/(p+1)**2,p) for p in good)
            if math.ceil(basic)>density:
                density,reason=math.ceil(basic),f'supersingular incidence at p={p}'
            refined=spectral_bound(n,good,counts,g)
            if refined and math.ceil(refined[0])>density:
                density,reason=refined
            if n == 83:
                # density83/series.py counts the five supersingular points on
                # X_0(83) with both Verschiebung lifts rational over F_4.
                split = 5 * (phi(n)//2)
                bound = math.ceil(Fraction(split,4))
                assert split == 205 and bound == 52 and bound <= g
                if density < bound:
                    density,reason = bound, 'supersingular points on an etale double cover at 2'
            if n in (117,123):
                # Projective bound for the last elliptic incidence factor.
                # The full period lattice at 117 is replayed by degree117.py.
                target, label, m, t, h = {
                    117: (72, '117.2.a.a', 288, -3, 5),
                    123: (80, '123.2.a.b', 160, -4, 3)
                }[n]
                remaining = [f for f in forms if n % f['level'] == 0 and occurs_in_positive_rank(f)]
                for f in remaining:
                    if f['label'] != label:
                        assert max(max(top_bound(f, p, S), S / max(p*p+1, (p+1)**2-beta[f['label']][p]))
                                   for p, S in counts.items()) > target
                d = target - 1
                assert counts[2] > 5*d and counts[2]*(h+1) > 5*d*(h+1)-t*m
                assert m > d and 2*((h-1)*d-m) < g and h*d-m < g
                density,reason = max(density,target), 'projective bound for elliptic incidence'
            if n == 117:
                # differential117.m gives six geometric supersingular zeros on X0.
                # Each has 36 lifts. Choosing the divisor through another point
                # forces one more coincidence, beyond the degree 10*72.
                assert counts[2] + 6 * (phi(n)//2) == 10*72
                density,reason = max(density,73), 'supersingular zeros of the 117a1 differential at 2'
            if n == 121:
                # The determinant bound descends a family to W_19 on X(11).
                # quotient121.py rules out the elliptic translate there.
                assert counts[2] == 605 and counts[2] // 5 == 121
                assert max(120*r - 50*math.comb(r,2) for r in range(1,12)) // 11 == 19
                density,reason = max(density,121), 'ramified descent to X(11) and degree-19 divisor sieve'
            if n == 138:
                # nodes138.py gives the -2 direction on the eight cyclic classes.
                largest_class = 4 * (phi(69)//2)
                assert g >= 108 and nodes - largest_class == 88
                if density < nodes - largest_class:
                    density,reason = nodes - largest_class, 'elliptic direction on the reduction graph at 2'
            if n in (231, 303):
                nodes = jordan(n // 3) // 12
                bound = igusa_elliptic_bound(n, nodes, counts[2])
                if density < bound:
                    density, reason = bound, 'moving elliptic projections on the Igusa components at 3'
            if n == 255 and density < level_255_endpoint():
                assert counts[2] == 256 * 9 and counts[2] > 256 * 5
                density,reason=level_255_endpoint(),'higher-dimensional equality obstruction'
        if n not in no_cm|unknown_cm and density<=cm_values[n]:
            density,reason=cm_values[n]+1,'CGPS sporadic CM point'
        if str(n) in inputs['proved_sporadic']:
            degree,source=inputs['proved_sporadic'][str(n)]
            if degree>=density:
                density,reason=degree+1,source+' sporadic'
        for m in divisors(n)[:-1]:
            if bounds[m]>density:
                density,reason=bounds[m],f'finite map to level {m}'
        bounds[n]=density
        degree,origin=witnesses[n]
        if n in negative:
            assert degree>=density
            status,why='no','small-degree torsion classifications and gonality'
        elif degree<density:
            status,why='yes','point below density bound'
        elif n in unit_cases:
            status,why='yes','cuspidal fibers and modular units'
        else:
            status,why='unresolved',''
        rows.append(dict(n=n,status=status,reason=why,CM_status_CGPS='no' if n in no_cm else 'undecided' if n in unknown_cm else 'yes',least_CM_degree=cm_values[n],noncuspidal_degree_upper=degree,point_source=origin,gonality_lower=g,density_lower=density,density_source=reason))
    return rows

if __name__=='__main__':
    rows=build_table()
    with (ROOT/'levels.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(Counter(r['status'] for r in rows))
    print('Unresolved:',*[r['n'] for r in rows if r['status']=='unresolved'])
