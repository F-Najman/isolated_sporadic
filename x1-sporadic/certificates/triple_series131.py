"""Write integral q-expansions of x, y and z on the Fricke model of X_0(131)."""
import sys
from pathlib import Path
from sage.all import EllipticCurve, LaurentSeriesRing, O, QQ, ZZ


def write_series(path, length=560):
    curve = EllipticCurve("131a1")
    assert curve.a_invariants() == (0,-1,1,1,0)
    ring = LaurentSeriesRing(QQ, "q", default_prec=length+30)
    q = ring.gen()
    x,y = (ring(v) for v in curve.modular_parametrization().power_series(length+30))
    eta = 1 + O(q**(length+20))
    for k in range(1,length+20):
        eta *= (1-q**k)**2
        if 131*k < length+20:
            eta *= (1-q**(131*k))**2
    z = ring(curve.q_eigenform(length+30))/(q**11*eta)
    assert (x.valuation(),y.valuation(),z.valuation()) == (-2,-3,-10)
    assert (y**2+y-x**3+x**2-x).valuation() > length
    with Path(path).open('w') as out:
        print(length,file=out)
        for series in (x,y,z):
            v = series.valuation()
            assert abs(series[v]) == 1 and series.prec() >= v+length
            print(*(ZZ(series[i]) for i in range(v,v+length)),file=out)


if __name__ == '__main__':
    write_series(sys.argv[1])
