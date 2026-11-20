"""Cusp-contact search over Q in the full spans of equal-weight two-factor products."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import reduce
from itertools import combinations, combinations_with_replacement
import json
from math import gcd, lcm
from pathlib import Path


def primitive(row):
    g = reduce(gcd, row, 0)
    if not g:
        return tuple(row)
    if next(v for v in row if v) < 0:
        g = -g
    return tuple(v//g for v in row)


def eliminate(v, r, p):
    if not v[p]:
        return v
    return primitive([r[p]*a-v[p]*b for a, b in zip(v, r)])


def extend(space, row):
    """Return RREF over Q with primitive integer rows. None forces a support zero."""
    for p, r in space:
        row = eliminate(row, r, p)
    if not any(row):
        return space
    row = primitive(row)
    p = next(i for i, v in enumerate(row) if v)
    if sum(v!=0 for v in row)==1:
        return None
    result = []
    for q, r in space:
        r = eliminate(r, row, p)
        if sum(v!=0 for v in r)==1:
            return None
        result.append((q,r))
    result.append((p,row))
    return tuple(sorted(result))


def null_vector(space, n):
    piv = {p for p, _ in space}
    frees = [j for j in range(n) if j not in piv]
    assert len(frees)==1
    q = frees[0]
    w = [Fraction(0)]*n
    w[q] = Fraction(1)
    for p, r in space:
        w[p] = Fraction(-r[q],r[p])
    m = lcm(*(v.denominator for v in w))
    assert all(w)
    return primitive([int(m*v) for v in w])


def expansion(pair, h, length):
    c = h if h%2 else h+83
    residues = [c*i%83 for i in pair]
    val = -(sum(r*(83-r) for r in residues)//166)
    coeff = [0]*length
    coeff[0] = (-1)**sum(c*i//83 for i in pair)
    for s in (v for r in residues for v in (r,83-r)):
        for k in range(s,length,83):
            for j in range(length-1,k-1,-1):
                coeff[j] -= coeff[j-k]
    return val, coeff


def flags_for(support, expansions):
    n = len(support)
    base, flags = 0, []
    for h, column in enumerate(expansions):
        local = [column[m] for m in support]
        shifts = [v for v, _ in local]
        low = min(shifts)
        base += shifts[0]-low
        if shifts.count(low)==1:
            continue
        shifts = [v-low for v in shifts]
        rows = []
        space = ()
        for e in range(len(local[0][1])):
            row = tuple(0 if x>e else f[e-x] for x, (_, f) in zip(shifts,local))
            if not any(row):
                continue
            next_space = extend(space,row)
            if next_space != space:
                rows.append((e,row))
                if next_space is None:
                    break
                space = next_space
        else:
            raise ValueError("insufficient cusp precision")
        assert rows[0][0]==0 and len(rows)<=n
        if rows[-1][0]>0:
            flags.append((h,rows))
    flags.sort(key=lambda z: -z[1][-1][0])
    return base, flags


def search(pairs, support, expansions, cutoff, answer, families, counts):
    base, flags = flags_for(support, expansions)
    n = len(support)
    maximum = [0]*(len(flags)+1)
    for t in range(len(flags)-1,-1,-1):
        maximum[t] = maximum[t+1]+flags[t][1][-1][0]
    if base-maximum[0]>cutoff:
        counts["excluded_by_individual_cusps"]+=1
        return
    seen = set()
    def recurse(t, space, need):
        if maximum[t]<need:
            return
        if len(space)==n-1:
            if space in seen:
                return
            seen.add(space)
            w = null_vector(space,n)
            contact = sum(next(e for e,row in rows if sum(a*b for a,b in zip(row,w)))
                          for _,rows in flags)
            d = base-contact
            assert d>=0
            if d<=cutoff:
                answer.append((d,[pairs[j] for j in support],w))
            return
        if t==len(flags):
            # A generic member realizes exactly the first jet outside the space
            # at every cusp, even when earlier branches also contain this space.
            if space in seen:
                return
            seen.add(space)
            total = sum(next(e for e,row in rows if extend(space,row)!=space)
                        for _,rows in flags)
            d = base-total
            assert 0<=d<=cutoff
            families.append((d,[pairs[j] for j in support],space))
            return
        rows = flags[t][1]
        for e,row in rows:
            next_space = extend(space,row)
            if next_space==space:
                continue
            if e+maximum[t+1]>=need:
                recurse(t+1,space,need-e)
            if next_space is None:
                break
            space = next_space
    recurse(0,(),base-cutoff)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("cutoff",type=int)
    cli.add_argument("output",type=Path)
    a = cli.parse_args()
    answer, families, counts = [], [], Counter()
    types = {}
    for s in (1,2):
        pairs = [t for t in combinations_with_replacement(range(1,42),2)
                 if sum(i*i for i in t)%83==s]
        assert len(pairs)==9+s
        expansions = [[expansion(t,h,200) for t in pairs] for h in range(1,42)]
        for size in range(1,len(pairs)+1):
            for support in combinations(range(len(pairs)),size):
                counts["supports"]+=1
                search(pairs,support,expansions,a.cutoff,answer,families,counts)
        types[s] = len(pairs)
    answer.sort()
    families.sort()
    data = dict(cutoff=a.cutoff,bucket_sizes=types,counts=counts,
                degrees=sorted(Counter(d for d,_,_ in answer).items()),
                points=answer, families=families)
    a.output.write_text(json.dumps(data,separators=(",",":"))+"\n")
    print("supports", counts, "sections", len(answer), "families",len(families))
    print("section residual degrees",data["degrees"])
    print("family residual degrees", sorted(Counter(d for d,_,_ in families).items()))


if __name__=="__main__":
    main()
