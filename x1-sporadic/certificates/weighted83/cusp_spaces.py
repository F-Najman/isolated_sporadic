"""Check cusp-supported Riemann--Roch classes near the degree-57 point."""
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from sage.all import GF, matrix


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('jets',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    n=83
    exponents={35:1,40:1,3:-1,34:-1}
    contacts={10:8,18:27,26:3}
    poles={}
    for c in range(1,42):
        order=sum(e*(i*c%n)*((i*c%n)-n) for i,e in exponents.items())
        assert order%(2*n)==0
        poles[c]=max(0,-order//(2*n))
    assert sum(poles.values())==95 and all(not poles[c] for c in contacts)
    jets=np.load(args.jets)
    assert int(jets['prime'])==167
    forms=matrix(GF(167),jets['rows'])
    base=[]
    for c,m in poles.items():
        w,prec,start=map(int,jets['profile'][c])
        assert w==1 and m<=prec
        base.extend(range(start,start+m))
    base_matrix=forms.matrix_from_rows(base)
    rank=base_matrix.rank()
    assert rank==94
    remaining=forms*base_matrix.right_kernel_matrix().transpose()
    results={}
    # For 2 points, the maximal 32-fold additions fit at every cusp. Also
    # test 3 points with a 12-fold addition for higher-support divisors.
    for size,mult in [(2,32),(3,12)]:
        blocks={}
        for c,pole in poles.items():
            w,prec,start=map(int,jets['profile'][c])
            extra=max(0,mult-contacts.get(c,0))
            assert pole+extra<=prec
            blocks[c]=range(start+pole,start+pole+extra)
        hist=Counter()
        undetermined=[]
        for cc in combinations(poles,size):
            rows=[i for c in cc for i in blocks[c]]
            mat=remaining.matrix_from_rows(rows)
            upper=1+len(base)+len(rows)-rank-mat.rank()
            hist[int(upper)]+=1
            # f-1 is a section; constants fail as at least one contact has
            # multiplicity greater than its allowed addition.
            assert any(m>(mult if c in cc else 0) for c,m in contacts.items())
            if upper!=2:
                undetermined.append((cc,int(upper)))
        results[str(size)]={'multiplicity':mult,'upper_bound_on_positive_part':dict(hist),
                            'unsettled':undetermined}
        print(size,mult,'Riemann--Roch upper bounds',dict(hist),flush=True)
    args.output.write_text(json.dumps(results,indent=1)+'\n')


if __name__=='__main__':
    main()
