"""List all primitive sections in the recorded low-degree strata up to a coefficient height."""
import argparse
from collections import Counter, defaultdict
from itertools import product
import json
from math import gcd
from pathlib import Path
from linear import primitive


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("strata",type=Path)
    cli.add_argument("height",type=int)
    cli.add_argument("prefix",type=Path)
    a=cli.parse_args()
    data=json.loads(a.strata.read_text())
    fixed=defaultdict(dict)
    families=defaultdict(list)
    for d,s,w in data["points"]:
        key=tuple(map(tuple,s))
        fixed[key][tuple(w)]=d
    for d,s,rows in data["families"]:
        key=tuple(map(tuple,s))
        assert len(s)-len(rows)==2 and len(rows)<=1
        families[key].append((d,rows))
    H=a.height
    parameters=[(u,v) for u in range(1,H+1) for v in range(-H,H+1) if v and gcd(u,v)==1]
    streams={i:Path(f"{a.prefix}-{i}.dat").open("w") for i in (2,3,4)}
    counts=Counter()
    for support in sorted(fixed.keys() | families.keys()):
        points={w:d for w,d in fixed[support].items() if max(map(abs,w))<=H}
        n=len(support)
        for d,space in families[support]:
            for u,v in parameters:
                if not space:
                    w=(u,v)
                else:
                    p,row=space[0]
                    j,k=(t for t in range(n) if t!=p)
                    vals=[0]*n
                    vals[j],vals[k]=row[p]*u,row[p]*v
                    vals[p]=-row[j]*u-row[k]*v
                    if not vals[p]:
                        continue
                    w=primitive(vals)
                    if max(map(abs,w))>H:
                        continue
                # Any jump from the generic degree of a pencil is a
                # one-dimensional coefficient space in the cusp enumeration.
                actual=fixed[support].get(w,d)
                if w in points:
                    assert points[w]==actual
                points[w]=actual
        for w,d in sorted(points.items()):
            if not d:
                continue
            assert 56<d<=data["cutoff"]
            print(*(i for pair in support for i in pair),*w,d,d,file=streams[n])
            counts[n]+=1
    for f in streams.values():
        f.close()
    print("height",H,"sections by number of terms",sorted(counts.items()))


if __name__=="__main__":
    main()
