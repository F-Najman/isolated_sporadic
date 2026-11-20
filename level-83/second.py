"""An isolation test and a residue-field comparison for the two degree-57 points."""
import argparse
import json
from pathlib import Path
import numpy as np
from sage.all import GF, PolynomialRing, ZZ, matrix, prime_range
from linear import expansion


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("jets",type=Path)
    cli.add_argument("minima",type=Path)
    cli.add_argument("output",type=Path)
    a=cli.parse_args()
    jets=np.load(a.jets)
    assert int(jets["prime"])==167
    forms=matrix(GF(167),jets["rows"])
    num,den=(2,9),(18,33)
    rows,contacts,poles=[],{},{}
    for c,(orbit,prec,start) in enumerate(jets["profile"]):
        orbit,prec,start=map(int,(orbit,prec,start))
        if c:
            v,A=expansion(num,c,90)
            w,B=expansion(den,c,90)
            poles[c]=max(0,w-v)
            if v==w:
                g=next(j for j in range(90) if A[j]+B[j])
                if g:
                    contacts[c]=g
        else:
            assert orbit==41  # the two-product leading values are independent
        m=max(0,poles.get(c,0)+1-contacts.get(c,0))
        assert m<=prec
        for t in range(orbit):
            rows.extend(range(start+t*prec,start+t*prec+m))
    degree=sum(poles.values())
    assert degree-sum(contacts.values())==57
    rank=int(forms.matrix_from_rows(rows).rank())
    upper=1+len(rows)-rank
    print("unit degree",degree,"cuspidal zeros",contacts,"deg B",len(rows),"rank",rank,"l(B) <=",upper)
    assert upper==2 and max(contacts.values())>1
    R=PolynomialRing(ZZ,"x")
    polys=[R(p) for p in json.loads(a.minima.read_text())["x_norms"]]
    assert len(polys)==2 and all(f.is_irreducible() for f in polys)
    irr=[None,None]
    split=None
    for prime in prime_range(2,5000):
        patterns=[]
        for i,poly in enumerate(polys):
            f=poly.change_ring(GF(prime))
            if f.degree()!=57 or f.gcd(f.derivative()).degree()!=0:
                patterns.append(None)
                continue
            parts=sorted(int(g.degree()) for g,e in f.factor())
            assert sum(parts)==57
            patterns.append(parts)
            if parts==[57] and irr[i] is None:
                irr[i]=int(prime)
        if split is None and all(patterns) and patterns[0]!=patterns[1]:
            split=(int(prime),patterns)
        if split and all(irr):
            break
    assert split
    rec=dict(unit_degree=degree,contacts=contacts,effective_degree=len(rows),
             jet_rank=rank,upper_bound=upper,irreducible_mod=irr,field_splitting=split)
    a.output.write_text(json.dumps(rec,indent=2)+"\n")
    print("irreducible reductions",irr,"different field splitting",split)


if __name__=="__main__":
    main()
