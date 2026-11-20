"""Rule out smaller closed points in individual Hecke correspondences on the two places."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from math import isqrt
from pathlib import Path
import subprocess


def is_prime(a):
    return a>=2 and all(a%j for j in range(2,isqrt(a)+1))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("engine",type=Path)
    p.add_argument("series",type=Path)
    p.add_argument("minima",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("--ells",default="2,3,5,7,11,13,17,19,23,29,31")
    p.add_argument("--jobs",type=int,default=4)
    a=p.parse_args()
    def run(ell,prime):
        return subprocess.check_output([str(a.engine.resolve()),str(ell),str(prime),
                                        str(a.series)],text=True).splitlines()
    polys=json.loads(a.minima.read_text())["x_norms"]
    prime=173
    for line,f in zip(run(1,prime),polys,strict=True):
        _,d,*monic=map(int,line.split())
        assert d==57 and monic==[c*pow(f[-1],-1,prime)%prime for c in f]
    def check(ell):
        possible=sorted({57*m//n for n in range(1,ell+2) for m in range(1,n)
                         if 57*m%n==0})
        assert possible and 1<=min(possible)<=max(possible)<=56
        pending=[set(possible),set(possible)]
        log=[]
        for prime in (n for n in range(5,20000) if n!=83 and (n-1)%ell==0 and is_prime(n)):
            patterns=[]
            for i,line in enumerate(run(ell,prime)):
                lab,d,*tail=line.split()
                d=int(d)
                assert int(lab)==i+1 and d==57*(ell+1)
                if not tail:
                    patterns.append(None)
                    continue
                modulus,fac=tail[0].split(":")
                assert len(tail)==1 and int(modulus)==prime
                factors=[]
                for token in fac.rstrip(",").split(","):
                    t,m=map(int,token.split("^"))
                    factors.extend([t]*m)
                assert sum(factors)==d
                acc=1
                for t in factors:
                    acc|=acc<<t
                pending[i]={k for k in pending[i] if (acc>>k)&1}
                patterns.append(sorted(factors))
            log.append(dict(prime=prime,factors=patterns))
            if not any(pending) or len(log)==32:
                break
        print(ell,"possible",possible,"tests",len(log),"remaining",list(map(sorted,pending)),flush=True)
        return dict(ell=ell,degrees=possible,reductions=log,remaining=list(map(sorted,pending)))
    with ThreadPoolExecutor(a.jobs) as pool:
        out=list(pool.map(check,map(int,a.ells.split(","))))
    a.output.write_text(json.dumps(out,indent=1)+"\n")


if __name__=="__main__":
    main()
