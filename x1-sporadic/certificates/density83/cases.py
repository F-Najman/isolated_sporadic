"""Generate the level-83 modular-unit input files and build the norm programs."""
import argparse
import subprocess
from collections import Counter
from itertools import combinations_with_replacement
from math import comb, prod
from pathlib import Path

import numpy as np
from sympy import isprime

P = 83
H = 41
SRC = Path(__file__).resolve().parent
SPLIT = 83000831
ZETA = 1495617


def least_pair(a,b):
    """In the zero square-sum class the least translate contains index 1."""
    for j in a+b:
        m = pow(j,-1,P)
        def image(indices):
            return tuple(sorted(min(m*y%P,-m*y%P) for y in indices))
        if tuple(sorted((image(a),image(b)))) < (a,b):
            return False
    return True


def small_cases(root, k, previous):
    cusps = np.array([a if a%2 else a+P for a in range(1,H+1)])
    values = [(2-pow(ZETA,i,SPLIT)-pow(ZETA,-i,SPLIT))%SPLIT for i in range(H+1)]
    buckets = [[] for _ in range(P)]
    for indices in combinations_with_replacement(range(1,H+1),k):
        buckets[sum(x*x for x in indices)%P].append(indices)
    units = 0
    hist = Counter()
    with (root/f'cases{k}.dat').open('w') as out:
        for s in (0,1,2):
            members = buckets[s]
            if not members:
                continue
            groups = np.asarray(members)
            products = groups[:,:,None]*cusps
            r = products % P
            orders = (r*(r-P)).sum(axis=1)
            signs = (products//P).sum(axis=1)%2
            contact = np.minimum(r,P-r).min(axis=1)
            leading = [prod(values[i] for i in t)%SPLIT for t in members]
            for i,a in enumerate(members):
                if s == 0 and a[0] != 1:
                    break
                ids = np.flatnonzero(np.all(~np.isin(groups[i+1:],a),axis=1))+i+1
                if s == 0:
                    ids = np.asarray([j for j in ids if least_pair(a,members[j])],dtype=int)
                if not len(ids):
                    continue
                assert all(leading[i]!=leading[j] for j in ids)
                val = orders[i]-orders[ids]
                assert not np.any(val%(2*P))
                val //= 2*P
                pole = np.maximum(0,-val).sum(axis=1)
                agree = signs[i]==signs[ids]
                cancel = (val==0)*np.minimum(contact[i],contact[ids])
                units += len(ids)
                for sign, matching in ((1,agree),(-1,~agree)):
                    deg = pole - (cancel*matching).sum(axis=1)
                    for j,d,t in zip(ids,deg,pole):
                        print(*a,*members[j],sign,d,t,file=out)
                        hist[int(d)] += 1
    assert sum(comb(len(b),2) for b in buckets) == H*(units+H*previous)
    assert (units,hist[0],min(d for d in hist if d)) == {
        2:(100,19,57), 3:(18132,495,57)}[k]
    return units


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workdir', type=Path)
    args = parser.parse_args()
    root = args.workdir.resolve()
    root.mkdir(exist_ok=True, parents=True)
    assert isprime(SPLIT) and pow(ZETA,P,SPLIT)==1 and ZETA!=1
    two = small_cases(root,2,0)
    small_cases(root,3,two)
    subprocess.run(['g++','-O3','-std=c++17',str(SRC/'four.cpp'),'-o',str(root/'four')],check=True)
    with (root/'four-hist.log').open('w') as out:
        subprocess.run([str(root/'four'),'389',str(root/'cases4.dat')],stdout=out,check=True)
    lines = (root/'four-hist.log').read_text().splitlines()
    assert lines[0]=='units 1876505 signed 3753010 prime-value collisions resolved 0 raw 110946697'
    hist = dict(tuple(map(int,line.split())) for line in lines[1:])
    assert hist[0]==2489 and min(d for d in hist if d)==57 and max(hist)==267
    assert sum(hist.values())==3753010
    libs = subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    for k in (2,3,4):
        subprocess.run(['g++','-O3','-std=c++17','-fopenmp',f'-DUNIT_FACTORS={k}',
                        str(SRC/'norm83.cpp'),'-o',str(root/f'norm{k}'),*libs,'-lgmpxx','-lgmp'],check=True)
    print('Generated the three searches in', root)


if __name__ == '__main__':
    main()
