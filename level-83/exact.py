"""Compute a primitive elliptic and x-line norm in characteristic zero for each section."""
import argparse
from collections import Counter
import json
from pathlib import Path
from sage.all import PolynomialRing, ZZ
from linear import expansion


R = PolynomialRing(ZZ, "x")
x = R.gen()
cubic = x**3+x**2+x
branch = (-19-39*x-31*x**2+3*x**3+16*x**4+8*x**5+x**6,
          -11-13*x+4*x**2+10*x**3+4*x**4)


def multiply(f,g):
    a,b = f
    c,e = g
    return a*c+cubic*b*e, a*e+b*c-(x+1)*b*e


def canonical(polys):
    c = ZZ(0)
    for p in polys:
        c = c.gcd(p.content())
    if next(p for p in polys if p).leading_coefficient()<0:
        c = -c
    return tuple(R([a//c for a in p]) for p in polys)


def make_basis(series,M):
    lines = series.read_text().splitlines()
    assert int(lines[0])>=M+12
    xs,ys,z = [R([int(v) for v in s.split()[:M+12]]) for s in lines[1:]]
    basis = [R.zero()]*(M+1)
    basis[0],basis[2],basis[3] = R.one(),xs,ys
    for j in range(4,M+1):
        basis[j]=(basis[j-2]*xs).truncate(M+11)
    return basis, [(b*z).truncate(M+11) for b in basis]


def trace_parts(f,basis,anti):
    d=len(f)-1
    def cancel(v,rows,anti_offset):
        a,b = R.zero(),R.zero()
        for j in range(d,anti_offset-1,-1):
            k=j-(6 if anti_offset else 0)
            if k<0 or k==1:
                assert v[d-j]==0
                continue
            row=rows[k]
            assert abs(row[0])==1
            coefficient=row[0]*v[d-j]
            if k%2:
                b+=coefficient*x**((k-3)//2)
            else:
                a+=coefficient*x**(k//2)
            if coefficient:
                for i in range(j+1):
                    v[d-j+i]-=coefficient*row[i]
        return a,b
    rem=f[:]
    W=cancel(rem,anti,1)
    f[-1]+=rem[-1]
    T=cancel(f,basis,0)
    t=multiply(T,T)
    w=multiply(branch,multiply(W,W))
    return canonical((t[0]-w[0],t[1]-w[1]))


def section_norm(pairs,weights,d,basis,anti):
    expansions=[[(val,R(jet)) for val,jet in (expansion(p,h,d+140) for p in pairs)]
                for h in range(1,42)]
    product=R.one()
    effective=0
    for cusp in expansions:
        low=min(v for v,f in cusp)
        local=sum(w*f.shift(v-low) for w,(v,f) in zip(weights,cusp))
        g=int(local.valuation())
        assert 0<=g<140
        effective+=cusp[0][0]-low-g
        product=(product*local.shift(-g)).truncate(d+1)
    assert effective==d
    # Invert prod_{83 does not divide n}(1-t^n)^r over Z.
    r=len(pairs[0])
    assert all(len(m)==r for m in pairs)
    for n in range(1,d+1):
        if n%83:
            inv=R([int(i%n==0) for i in range(d+1)])
            product=(product*(inv**r).truncate(d+1)).truncate(d+1)
    a,b=trace_parts([product[i] for i in range(d+1)],basis,anti)
    (norm,)=canonical((a*(a-(x+1)*b)-cubic*b*b,))
    assert norm.degree()==d
    return a,b,norm


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("sections",type=Path)
    cli.add_argument("series",type=Path)
    cli.add_argument("output",type=Path)
    a=cli.parse_args()
    entries=json.loads(a.sections.read_text())["points"]
    basis,anti=make_basis(a.series,max(d for d,_,_ in entries))
    cache={}
    records=[]
    for d,pairs,w in entries:
        if not d:
            continue
        e0,e1,norm=section_norm(pairs,w,d,basis,anti)
        if norm not in cache:
            cache[norm]=sorted((int(f.degree()),int(m)) for f,m in norm.factor())
        ell=[[int(v) for v in f] for f in (e0,e1)]
        records.append(dict(degree=d,pairs=pairs,weights=w,elliptic=ell,
                            x_norm=[int(c) for c in norm],factors=cache[norm]))
    a.output.write_text(json.dumps(records,separators=(",",":"))+"\n")
    print("sections",len(records),"distinct x-norms",len(cache))
    print("factors of degrees at most 56",
          [(r["degree"],r["factors"]) for r in records if any(t<=56 for t,_ in r["factors"])])
    print("distinct elliptic norms at degree 57",
          len({str(r["elliptic"]) for r in records if r["degree"]==57}))


if __name__=="__main__":
    main()
