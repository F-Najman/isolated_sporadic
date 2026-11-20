"""Weighted LLL on products of E_i supported at the rational cusps of X_1(131)."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from random import Random

import numpy as np
from fpylll import IntegerMatrix, LLL
from vanhoeij131 import CUSPS, N, fibre, folded, sieve

indices = np.arange(1, 66)
odd = np.where(indices % 2, indices, indices+N)
products = indices[:, None] * odd
residues = products % N
orders = residues*(residues-N)
parities = (products//N) % 2
contacts = np.minimum(residues, N-residues)


def congruence_basis(allowed):
    # Sum e_i = 0 and sum i^2 e_i = 0 mod 131.
    b, a, *rest = allowed
    inv = pow(a*a-b*b, -1, N)
    rows = []
    first = np.zeros(65, dtype=np.int64)
    first[a-1], first[b-1] = N, -N
    rows.append(first)
    for i in rest:
        t = (i*i-b*b)*inv % N
        if t > 65:
            t -= N
        row = np.zeros(65, dtype=np.int64)
        row[i-1], row[a-1], row[b-1] = 1, -t, t-1
        rows.append(row)
    return np.array(rows)


def canonical(unit):
    half = len(unit)//2
    forms = []
    for k in CUSPS:
        pos = tuple(sorted(folded(k*i) for i in unit[:half]))
        neg = tuple(sorted(folded(k*i) for i in unit[half:]))
        forms.append(min(pos,neg) + max(pos,neg))
    return min(forms)


def search(runs, seed):
    rng = Random(seed)
    pool = set()
    examined = 0
    for run in range(runs):
        nforced = run % 5
        forced = ([1] + rng.sample(range(2,66), nforced-1)) if nforced else []
        cutoffs = [rng.choice((1,5,9,13,17,21,25,29,33)) for c in forced]
        allowed = [i for i in CUSPS if all(folded(i*c) >= m for c,m in zip(forced,cutoffs))]
        if len(allowed) < 3:
            continue
        basis = congruence_basis(allowed)
        divs = basis @ orders
        assert not np.any(divs % (2*N))
        divs //= 2*N
        weights = np.array([rng.choice((2,3,4,5,7,9)) for c in CUSPS])
        for c in forced:
            weights[c-1] = 2048
        lattice = IntegerMatrix.from_matrix(np.hstack((divs*weights, basis)).tolist())
        LLL.reduction(lattice, delta=.97)
        out = np.empty((len(basis),130), dtype=np.int64)
        for j in range(len(basis)):
            out[j] = list(lattice[j])
        for row in out:
            e = row[65:]
            if np.sum(abs(e)) <= 4 or np.gcd.reduce(abs(e)) != 1:
                continue
            v = e @ orders // (2*N)
            assert np.array_equal(row[:65], v*weights)
            if any(v[c-1] for c in forced):
                continue
            examined += 1
            degrees = int(np.maximum(0,-v).sum())
            s = (e @ parities) % 2
            hit = np.min(contacts[e!=0], axis=0)
            least = min(degrees - int(hit[(v==0) & (s==k)].sum()) for k in (0,1))
            if least <= 190:
                pos = tuple(i for i,t in enumerate(e,1) for _ in range(max(t,0)))
                neg = tuple(i for i,t in enumerate(e,1) for _ in range(max(-t,0)))
                pool.add(canonical(pos + neg))
    cases = []
    for u in sorted(pool):
        for sign in (1,-1):
            degree, residual, _ = fibre(u,sign)
            if 0 < residual <= 190:
                cases.append((residual,u,sign,degree))
    return examined, sorted(cases)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', type=int, default=512)
    parser.add_argument('--seed', type=int, default=131)
    parser.add_argument('--workers', type=int, default=4)
    opts = parser.parse_args()
    examined, cases = search(opts.runs, opts.seed)
    print('weighted LLL runs', opts.runs, 'retained rows with repetitions', examined,
          'distinct positive fibres at most 190', len(cases), flush=True)
    print('residual degrees', sorted(Counter(c[0] for c in cases).items()), flush=True)
    aux = (10,65,11,55)
    irreducible = 0
    with ProcessPoolExecutor(max_workers=opts.workers) as pool:
        for result in pool.map(sieve, [(u,s,aux,24) for d,u,s,degree in cases]):
            print(result, flush=True)
            irreducible += not result[3]
    print('proved irreducible', irreducible, 'out of', len(cases))
