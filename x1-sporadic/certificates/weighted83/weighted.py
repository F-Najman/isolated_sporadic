"""Find modular units by weighted LLL on cusp valuations at level 83."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from math import gcd
from pathlib import Path
from random import Random

import numpy as np

N = 83
H = 41
INDICES = np.arange(1,H+1)
PRODUCT = INDICES[:,None] * (INDICES + N*(1-INDICES%2))
R = PRODUCT % N
ORDER = 12*R*(R-N) + 2*N*N
DEN = 24*N
FOLD = np.minimum(R,N-R)
PARITY = PRODUCT//N % 2
PLACE = 83000831
ZETA = 1495617


def congruences(allowed, balanced):
    i,j,*tail = map(int,allowed)
    rows = []
    def row(values):
        a = np.zeros(H,dtype=np.int64)
        for index, value in values.items():
            a[index-1]=value
        rows.append(a)
    if balanced:
        row({i:N,j:-N})
        inverse = pow(j*j-i*i,-1,N)
        for k in tail:
            t = (k*k-i*i)*inverse%N
            row({k:1,j:-t,i:t-1})
    else:
        row({i:12*N})
        for k in [j,*tail]:
            t = -k*k*pow(i*i,-1,N)%N
            t += N*((-1-t)*pow(N,-1,12)%12)
            row({k:1,i:t})
    return np.array(rows)


def representative(e):
    # The diamonds permute exponents, changing the unit by at most a sign.
    # Positive total places zeros at the nonsplit cusp; include reciprocal
    # targets below so that inverting a balanced unit also preserves the search.
    images = e[(FOLD-1).T]
    if e.sum()<0:
        images = -images
    if e.sum()==0:
        images = np.concatenate((images,-images))
    # Exponents here fit signed bytes; lexicographic comparison is only an
    # ordering convention, not a numerical test.
    assert abs(images).max()<127
    return min(row.tobytes() for row in images.astype(np.int8))


def block(spec):
    from fpylll import IntegerMatrix, LLL
    seed, count, cutoff = spec
    rng = Random(seed)
    retained = set()
    selected = 0
    seen = set()
    for run in range(count):
        balanced = run%3 == 0
        number = run%5
        forced = ([0]+rng.sample(range(1,H),number-1)) if number else []
        lower = [rng.choice((1,1,3,5,8,11,15,19,23,27)) for _ in forced]
        allowed = [i for i in range(1,H+1)
                   if all(FOLD[i-1,c]>=m for c,m in zip(forced,lower))]
        if len(allowed)<3:
            continue
        basis = congruences(allowed,balanced)
        nums = basis@ORDER
        assert not np.any(nums%DEN)
        valuations = nums//DEN
        weights = np.array([rng.choice((1,2,3,4,6,8,11)) for _ in range(H)])
        weights[forced] = 2048
        weighted = np.column_stack((valuations*weights,
                                    (basis.sum(axis=1)//12)*rng.choice((8,16,32,64)),
                                    basis))
        lattice = IntegerMatrix.from_matrix(weighted.tolist())
        LLL.reduction(lattice,delta=.97)
        exponents = np.asarray([list(row)[-H:] for row in lattice],dtype=np.int64)
        # Neighbouring reduced basis vectors can have overlapping poles.
        m = min(len(exponents),6)
        rows = [*exponents]
        for i in range(m):
            for j in range(i):
                rows.extend((exponents[i]+exponents[j],exponents[i]-exponents[j]))
        for e in rows:
            key=e.tobytes()
            if key in seen:
                continue
            seen.add(key)
            if np.gcd.reduce(abs(e))!=1:
                continue
            total = int(e.sum())
            nums = e@ORDER
            assert total%12==0 and not np.any(nums%DEN)
            vals=nums//DEN
            if any(vals[c] for c in forced):
                continue
            pole = int(np.maximum(-vals,0).sum()) + H*max(0,-total//12)
            if not pole or pole>cutoff:
                continue
            selected+=1
            retained.add(representative(e))
    return selected,retained


@lru_cache(maxsize=None)
def cyclotomic(c,i):
    a = i*c%N
    return (2-pow(ZETA,a,PLACE)-pow(ZETA,-a,PLACE))%PLACE


def nonsplit_check(e,targets):
    # A contact at a nonsplit cusp with value a/b would force the square of
    # its absolute value, at every cyclotomic conjugate, to equal (a/b)^2.
    bad=set(targets)
    for c in range(1,H+1):
        value=1
        for i in np.flatnonzero(e):
            value=value*pow(cyclotomic(c,int(i)+1),int(e[i]),PLACE)%PLACE
        bad={t for t in bad if (value*t[1]**2-t[0]**2)%PLACE==0}
        if not bad:
            return
    raise ArithmeticError(('unresolved nonsplit cusp',list(e),bad))


def write_units(root,pool,height):
    targets=[(a,b) for a in range(-height,height+1) for b in range(1,height+1)
             if a and gcd(a,b)==1]
    histogram=Counter()
    with (root/'units.txt').open('w') as units, (root/'fibres.txt').open('w') as fibres:
        for label,rep in enumerate(sorted(pool)):
            e=np.frombuffer(rep,dtype=np.int8).astype(np.int64)
            total=int(e.sum())
            assert total>=0 and total%12==0 and (e@(INDICES**2))%N==0
            assert not np.any((e@ORDER)%DEN)
            if total==0:
                nonsplit_check(e,targets)
            vals=e@ORDER//DEN
            pole=int(np.maximum(-vals,0).sum())
            parity=e@PARITY%2
            contact=FOLD[e!=0].min(axis=0)
            print(label,*e,file=units)
            family = ('unbalanced' if total else
                      ('balanced-small' if abs(e).sum()<=8 else 'balanced-long'))
            histogram['units',family]+=1
            for a,b in targets:
                cancel=(int(contact[(vals==0)&(parity==int(a<0))].sum())
                        if abs(a)==b else 0)
                d=pole-cancel
                assert d>=0
                print(label,a,b,pole,d,file=fibres)
                histogram['fibres',family]+=1
                if d:
                    histogram['positive',family]+=1
                    if d<57:
                        print('residual of degree below 57:',label,a,b,d,flush=True)
    print('counts',sorted(histogram.items()),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out',type=Path)
    parser.add_argument('--runs',type=int,default=768)
    parser.add_argument('--jobs',type=int,default=8)
    parser.add_argument('--height',type=int,default=4)
    parser.add_argument('--cutoff',type=int,default=175)
    parser.add_argument('--from-units',type=Path,
                        help='use an explicit saved catalogue in place of performing LLL')
    args=parser.parse_args()
    args.out.mkdir(exist_ok=True,parents=True)
    if args.from_units:
        pool=set()
        for label,line in enumerate(args.from_units.read_text().splitlines()):
            i,*values=map(int,line.split())
            e=np.array(values,dtype=np.int64)
            assert i==label and len(e)==H and max(abs(e))<127
            rep=representative(e)
            assert rep==e.astype(np.int8).tobytes() and rep not in pool
            pool.add(rep)
        write_units(args.out,pool,args.height)
        raise SystemExit
    batches=[(830000+i,min(48,args.runs-i*48),args.cutoff)
             for i in range((args.runs+47)//48)]
    pool=set()
    selected=0
    with ProcessPoolExecutor(args.jobs) as executor:
        for n,units in executor.map(block,batches):
            pool.update(units)
            selected+=n
    # Include the known degree-57 point as a positive control.
    control=np.zeros(H,dtype=np.int64)
    control[[34,39]]=1
    control[[2,33]]=-1
    pool.add(representative(control))
    print('LLL attempts',args.runs,'selected rows including repeats',selected,
          'diamond and inversion representatives',len(pool),flush=True)
    write_units(args.out,pool,args.height)
