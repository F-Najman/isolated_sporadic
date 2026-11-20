"""Compare the elliptic x-coordinates of all degree-57 sections in the new enumerations."""
import argparse
from collections import Counter
import json
from pathlib import Path
from exact import make_basis,section_norm


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for field in ("sections","five","three","series","output"):
        p.add_argument(field,type=Path)
    p.add_argument("--fourterms",type=Path)
    a=p.parse_args()
    inputs=[]
    # Start with the example in the earlier level-83 note.
    inputs.append(("original",[(3,34),(35,40)],[1,-1]))
    for d,s,w in json.loads(a.sections.read_text())["points"]:
        if d==57:
            inputs.append(("two-product",list(map(tuple,s)),w))
    with a.five.open() as f:
        for line in f:
            vals=list(map(int,line.split()))
            if vals[-2]==57:
                inputs.append(("five-factor",[tuple(vals[:5]),tuple(vals[5:10])],[1,-vals[10]]))
    with a.three.open() as f:
        for line in f:
            vals=list(map(int,line.split()))
            if vals[-1]==57:
                inputs.append(("three-term",[tuple(vals[i:i+3]) for i in (0,3,6)],vals[9:12]))
    if a.fourterms:
        with a.fourterms.open() as f:
            for line in f:
                vals=list(map(int,line.split()))
                if vals[-1]==57:
                    inputs.append(("four-term",[tuple(vals[i:i+3]) for i in (0,3,6,9)],vals[12:16]))
    basis,anti=make_basis(a.series,57)
    norms={}
    records=[]
    for family,s,w in inputs:
        a0,b0,f=section_norm(s,w,57,basis,anti)
        key=tuple(map(int,f.list()))
        if key not in norms:
            assert f.is_irreducible()
            norms[key]=len(norms)+1
        records.append(dict(family=family,products=s,weights=w,x_class=norms[key]))
    a.output.write_text(json.dumps(dict(x_norms=list(norms),sections=records),
                                   separators=(",",":"))+"\n")
    print("distinct irreducible x-norms",len(norms))
    print("sections by class and family",sorted(Counter((r["x_class"],r["family"]) for r in records).items()))
    for k in range(1,len(norms)+1):
        print(k,next(r for r in records if r["x_class"]==k))


if __name__=="__main__":
    main()
