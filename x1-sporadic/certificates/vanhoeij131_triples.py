"""Enumerate the signed three-over-three Siegel-unit fibres of X_1(131)."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations_with_replacement
from math import comb
from pathlib import Path

import numpy as np
from sympy import isprime
from vanhoeij131 import CUSPS, N, folded, sieve


def buckets():
    groups = [[] for _ in range(N)]
    for triple in combinations_with_replacement(CUSPS, 3):
        groups[sum(i*i for i in triple) % N].append(triple)
    return groups


def zero_representative(a, b):
    # A lexicographically minimal translate has an index 1. Only translate
    # by inverses of indices of the input to obtain every such translate.
    for i in set(a+b):
        k = pow(i, -1, N)
        first, second = sorted(tuple(sorted(folded(k*t) for t in side)) for side in (a,b))
        if (first, second) < (a,b):
            return False
    return True


def split_prime():
    q = 131*1000000+1
    while not isprime(q):
        q += 131
    for x in range(2, q):
        zeta = pow(x, (q-1)//N, q)
        if zeta != 1:
            assert pow(zeta, N, q) == 1
            return q, zeta


def enumerate_fibres():
    groups = buckets()
    prime, zeta = split_prime()
    sin2 = [(2 - pow(zeta,i,prime) - pow(zeta,-i,prime)) % prime for i in range(66)]
    cusps = np.array([c if c%2 else c+N for c in CUSPS])
    # The two nonzero square classes have representatives 1 and 2.
    assert pow(2,65,N) == N-1
    count = 0
    for norm in (0,1,2):
        members = groups[norm]
        triples = np.array(members)
        numer = triples[:,:,None]*cusps
        rs = numer % N
        valuations = (rs*(rs-N)).sum(axis=1)
        parities = (numer//N).sum(axis=1) % 2
        first_order = np.minimum(rs,N-rs).min(axis=1)
        high = np.array([sin2[i]*sin2[j]*sin2[k] % prime for i,j,k in members])
        for j, a in enumerate(members):
            if norm == 0 and a[0] != 1:
                break
            disjoint = np.all(~np.isin(triples[j+1:], a),axis=1)
            ids = np.flatnonzero(disjoint) + j+1
            if norm == 0:
                ids = np.array([k for k in ids if zero_representative(a,members[k])], dtype=int)
            if not len(ids):
                continue
            count += len(ids)
            # The values of |product E_i|^2 at a nonrational cusp are
            # products of 2-zeta^i-zeta^(-i), up to a common normalization.
            # Unequal reductions exclude both signs in characteristic zero.
            assert np.all(high[ids] != high[j])
            vals = valuations[j]-valuations[ids]
            assert not np.any(vals % (2*N))
            vals //= 2*N
            poles = np.maximum(0,-vals).sum(axis=1)
            zero = (vals == 0)
            same = (parities[j] == parities[ids])
            order = np.minimum(first_order[j],first_order[ids])
            plus = poles-(zero*same*order).sum(axis=1)
            minus = poles-(zero*~same*order).sum(axis=1)
            for k, pole, pos, neg in zip(ids,poles,plus,minus):
                unit = a+members[k]
                yield unit, 1, int(pos), int(pole)
                yield unit, -1, int(neg), int(pole)
    # With a common index, cancellation gives a nonconstant two-over-two
    # unit and a unique cancelled index. Every orbit has size 65.
    raw = sum(comb(len(group),2) for group in groups)
    assert raw-65*65*256 == 65*count


def read_cases(path):
    with Path(path).open() as f:
        for line in f:
            a,b,c,i,j,k,s,d,pole = map(int,line.split())
            yield (a,b,c,i,j,k),s,d,pole


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-cases', type=Path)
    parser.add_argument('--read-cases', type=Path)
    parser.add_argument('--max-residual', type=int)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    src = read_cases(args.read_cases) if args.read_cases else enumerate_fibres()
    selected = []
    counts = Counter()
    out = args.write_cases.open('w') if args.write_cases else None
    for unit, sign, degree, pole in src:
        assert degree >= 0
        counts[degree] += 1
        if out:
            print(*unit,sign,degree,pole,file=out)
        if args.max_residual and 0 < degree <= args.max_residual:
            selected.append((degree,unit,sign))
    if out:
        out.close()
    print('signed fibres', sum(counts.values()), 'cuspidal', counts[0],
          'positive degree range', min(d for d in counts if d), max(counts),flush=True)
    print('counts through each cutoff', [(b,sum(n for d,n in counts.items() if 0<d<=b))
                           for b in (129,150,160,180,200,220,250,260,300,350,400,450,500)],flush=True)
    if selected:
        selected.sort()
        h = (10,65,11,55)
        tasks = [(u,s,h,24) for d,u,s in selected]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for (d,u,s), result in zip(selected,pool.map(sieve,tasks,chunksize=1)):
                assert (u,s,d) == result[:3]
                print(result,flush=True)


if __name__ == '__main__':
    main()
