"""Certify the isogeny-class degree using integral reductions on the elliptic quotient."""
import argparse
from math import gcd
import json
from pathlib import Path
from sage.all import GF, PolynomialRing, ZZ


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("norms",type=Path)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    R=PolynomialRing(ZZ,"x")
    rows=[]
    for r in json.loads(a.norms.read_text()):
        d=r["degree"]
        f=R(r["x_norm"])
        assert d==f.degree() and f.content()==1 and f.is_irreducible()
        bound=d
        evidence=[]
        for prime,_ in f.leading_coefficient().factor():
            if prime<=3 or prime==83:
                continue
            drop=int(d-f.change_ring(GF(prime)).degree())
            assert 1<=drop<=d
            evidence.append([int(prime),drop])
            bound=gcd(bound,drop)
        assert bound==1
        rows.append(dict(degree=d,products=r["pairs"],weights=r["weights"],cusp_counts=evidence))
    a.output.write_text(json.dumps(rows,indent=1)+"\n")
    print("sections",len(rows),"with a singleton cusp reduction",
          sum(any(c==1 for p,c in r["cusp_counts"]) for r in rows))


if __name__=="__main__":
    main()
