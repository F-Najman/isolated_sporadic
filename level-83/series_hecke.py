"""Higher-precision expansion of the same elliptic and Fricke functions as in series.py."""
import argparse
from pathlib import Path
from sage.all import EllipticCurve, LaurentSeriesRing, O, QQ, ZZ


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("precision",type=int)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    n=a.precision
    E=EllipticCurve("83a1")
    ring=LaurentSeriesRing(QQ,"q",default_prec=n+64)
    q=ring.gen()
    x,y=map(ring,E.modular_parametrization().power_series(n+52))
    # Euler's pentagonal-number theorem.
    def euler(s):
        value=1+O(q**(n+35))
        k=1
        while s*k*(3*k-1)//2<n+35:
            value+=(-1)**k*(q**(s*k*(3*k-1)//2)+q**(s*k*(3*k+1)//2))
            k+=1
        return value
    z=ring(E.q_eigenform(n+40))/(q**7*euler(1)**2*euler(83)**2)
    branch=(-19-39*x-31*x**2+3*x**3+16*x**4+8*x**5+x**6
            +(-11-13*x+4*x**2+10*x**3+4*x**4)*y)
    assert (branch-z*z).valuation()>n-20
    with a.output.open("w") as f:
        print(n,file=f)
        for func,v in [(x,-2),(y,-3),(z,-6)]:
            assert func.valuation()==v
            print(*(ZZ(func[v+i]) for i in range(n)),file=f)


if __name__=="__main__":
    main()
