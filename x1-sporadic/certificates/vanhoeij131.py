"""Search the two-over-two Siegel-unit fibres U=+1,-1 on X_1(131).

Run with Sage Python. Build the accompanying unit_traces.cpp as unit_traces
in this directory, or set X1_UNIT_TRACES to the path of the executable.
"""
import argparse
import os
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, islice
from pathlib import Path
from subprocess import check_output

from cypari2 import Pari
from sympy import primerange
from siegel_contacts import _leading, _zero_at_zeta

gp = Pari()
ENGINE = os.environ.get('X1_UNIT_TRACES', str(Path(__file__).resolve().with_name('unit_traces')))

N = 131
CUSPS = tuple(range(1, 66))


def folded(i):
    return min(i % N, -i % N)


def representatives():
    groups = defaultdict(list)
    for a, b in combinations_with_replacement(CUSPS, 2):
        groups[(a*a + b*b) % N].append((a, b))
    seen = set()
    result = []
    for group in groups.values():
        for first, second in combinations(group, 2):
            assert set(first).isdisjoint(second)
            if (first, second) in seen:
                continue
            for k in CUSPS:
                pair = [tuple(sorted(folded(k*t) for t in side)) for side in (first, second)]
                seen.add(tuple(sorted(pair)))
            result.append(first + second)
    return result


def exponents(unit):
    assert len(unit) % 2 == 0
    half = len(unit)//2
    e = Counter(unit[:half])
    e.subtract(unit[half:])
    return dict(e)


@lru_cache(maxsize=32)
def local_data(unit):
    e = exponents(unit)
    # The nonrational cusps form one orbit. Check at one geometric cusp
    # that neither sign is attained, using its value in Q(zeta_262).
    s2 = sum(i*i*t for i,t in e.items())
    top, bottom = _leading(N, N, 1, e, s2)
    for sign in (1, -1):
        delta = Counter(top)
        delta.subtract({i: sign*x for i,x in bottom.items()})
        assert not _zero_at_zeta(2*N, delta)
    out = []
    for c in CUSPS:
        a = c if c % 2 else c+N
        numer = 0
        parity = 0
        powers = [0]*N
        for i, exp in e.items():
            quotient, r = divmod(a*i, N)
            numer += exp*r*(r-N)
            parity += exp*quotient
            powers[r] += exp
            powers[N-r] += exp
        v, rem = divmod(numer, 2*N)
        assert rem == 0
        out.append((v, (-1)**(parity % 2), powers))
    return out


def fibre(unit, sign):
    data = local_data(unit)
    pole_degree = sum(max(-v, 0) for v, _, _ in data)
    orders = []
    for v, leading, powers in data:
        if v == 0 and leading == sign:
            # The least folded index with nonzero exponent occurs once.
            orders.append(next(k for k in CUSPS if powers[k]))
        else:
            orders.append(min(0, v))
    return pole_degree, -sum(orders), orders


def prepare(unit, sign, auxiliary):
    _, degree, orders = fibre(unit, sign)
    constants = [0]*degree
    rows = []
    lower = degree + 1
    for (v, sf, ef), (a, sh, eh), w in zip(local_data(unit), local_data(auxiliary), orders):
        if a == 0:
            for j in range(degree):
                constants[j] -= w * sh**(j + 1)
        if a < 0:
            lower = max(lower, w - min(v, 0) - a*degree + 4)
            rows.append((v, sf*sign, w, -a, sh, *ef, *eh))
    suffix = '\n'.join(' '.join(map(str, row)) for row in rows) + '\n'
    return degree, constants, len(rows), suffix, lower


def sieve(args):
    unit, sign, auxiliary, prime_count = args
    degree, constants, nrows, suffix, lower = prepare(unit, sign, auxiliary)
    # Every proper factorization has a monic factor of degree at most d/2.
    cutoff = degree//2
    possible = (1 << (cutoff + 1)) - 2
    witnesses = []
    for p in islice(primerange(lower + 1, 1000000), prime_count):
        raw = (f'{N} {degree} {p} {nrows}\n' +
               ' '.join(str(v % p) for v in constants) + '\n' + suffix)
        out = check_output([ENGINE], input=raw.encode())
        poly = gp.Pol([int(v) for v in out.split()])
        assert poly.poldegree() == degree
        factors = poly.factormod(p, 1)
        pattern = [(int(d), int(m)) for d,m in zip(factors[0], factors[1])]
        assert sum(d*m for d,m in pattern) == degree
        subsets = 1
        for d, m in pattern:
            for _ in range(m):
                subsets |= subsets << d
        new = possible & subsets
        if new != possible:
            witnesses.append((p, pattern))
        possible = new
        if not possible:
            break
    left = [d for d in range(1, cutoff+1) if possible & (1 << d)]
    return unit, sign, degree, left, witnesses


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--primes', type=int, default=24)
    parser.add_argument('--max-residual', type=int)
    parser.add_argument('--summary', action='store_true')
    options = parser.parse_args()
    units = representatives()
    aux = min(units, key=lambda u: (max(abs(v) for v,_,_ in local_data(u)),
                                    sum(abs(v) for v,_,_ in local_data(u))))
    cases = []
    cuspidal = 0
    for u in units:
        for sign in (1, -1):
            _, d, _ = fibre(u, sign)
            assert d >= 0
            if d == 0:
                cuspidal += 1
                continue
            if options.max_residual is None or d <= options.max_residual:
                shift = 3 if u == aux else 1
                h = tuple(folded(i*shift) for i in aux)
                cases.append((d, (u, sign, h, options.primes)))
    cases.sort()
    assert len(units) == 256 and cuspidal == 31
    if options.max_residual is None:
        assert len(cases) == 481 and (cases[0][0], cases[-1][0]) == (143,328)
    if not cases:
        print('unit classes', len(units), 'cuspidal fibres', cuspidal,
              'no positive residual fibre in the selected range')
        return
    print('unit classes', len(units), 'cuspidal fibres', cuspidal, 'selected fibres', len(cases),
          'range', (cases[0][0], cases[-1][0]), 'auxiliary', aux, flush=True)
    u, sign, _, _ = cases[0][1]
    degree, residual, orders = fibre(u, sign)
    print('first fibre', u, sign, 'degree', degree, 'residual', residual,
          'cusp contacts', [(c,w) for c,w in zip(CUSPS,orders) if w > 0], flush=True)
    if options.summary:
        print('counts by upper cutoff', [(n, sum(d<=n for d,_ in cases))
                                          for n in (167,180,200,225,260,350,500)], flush=True)
        return
    args = [case for _,case in cases]
    irreducible = 0
    with ProcessPoolExecutor(max_workers=options.workers) as pool:
        for result in pool.map(sieve, args, chunksize=1):
            print(result, flush=True)
            irreducible += not result[3]
    print('proved irreducible', irreducible, 'out of', len(cases))
    assert irreducible == len(cases)


if __name__ == '__main__':
    main()
