"""sage -python series.py [--units WORKDIR]

Check the supersingular images on X_0^+(83) in characteristic two.
With --units, also verify the Fricke equation and write the input series
for the modular-unit norm computation.
"""
import argparse
from pathlib import Path

from sage.all import EllipticCurve, GF, LaurentSeriesRing, O, QQ, ZZ
from sage.modular.modform.j_invariant import j_invariant_qexp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--units', type=Path)
    args = parser.parse_args()
    prec = 390 if args.units else 135
    E = EllipticCurve('83a1')
    assert list(E.a_invariants()) == [1,1,1,1,0]
    ring = LaurentSeriesRing(QQ, 'q', default_prec=prec+50)
    q = ring.gen()
    x, y = [ring(w) for w in E.modular_parametrization().power_series(prec+50)]
    j = ring(j_invariant_qexp(prec+20))
    # At either cusp of X_0(83), j times its Fricke conjugate has pole 84.
    jq83 = q**-83 + 744 + 196884*q**83 + O(q**166)
    delta = j*jq83 - (x*(x+1)*y)**12
    assert delta.precision_absolute() > 0
    assert all(ZZ(delta[n]) % 2 == 0 for n in range(-84,1))
    K = GF(4, 'a')
    a = K.gen()
    E4 = E.change_ring(K)
    assert E.change_ring(GF(2)).cardinality() == 4
    assert E4.cardinality() == 8 and max(P.order() for P in E4.points()) == 8
    expected = {(K(0),K(0)):2, (K(0),K(1)):1, (K(1),K(1)):2,
                (a,K(0)):1, (a**2,K(0)):1}
    # For x and x+1, the divisor counts P and -P. For y, use
    # x(x^2+x+1)=0, whose roots are simple. Their combined degree is seven.
    actual = {}
    for P in E4.points():
        if P.is_zero():
            continue
        u,v = P.xy()
        order = (1+int(P == -P) if u in (0,1) else 0) + int(v == 0)
        if order:
            actual[u,v] = order
    assert actual == expected and sum(actual.values()) == 7
    split = sum(m for (u,v),m in actual.items() if u**2 == u and v**2 == v)
    assert split == 5 and (41*split+3)//4 == 52
    print('Seven supersingular images with multiplicity:', actual)
    print('Five lie above rational values on E; split lifts on X_1(83):', 41*split)
    if not args.units:
        return

    eta = q**7 * (1+O(q**(prec+30)))
    for n in range(1,prec+30):
        eta *= (1-q**n)**2
        if 83*n < prec+30:
            eta *= (1-q**(83*n))**2
    z = ring(E.q_eigenform(prec+40)) / eta
    assert [f.valuation() for f in (x,y,z)] == [-2,-3,-6]
    A = [-19,-39,-31,3,16,8,1]
    B = [-11,-13,4,10,4]
    branch = sum(c*x**i for i,c in enumerate(A)) + y*sum(c*x**i for i,c in enumerate(B))
    assert (z*z-branch).valuation() > 100
    args.units.mkdir(exist_ok=True, parents=True)
    with (args.units/'series.dat').open('w') as stream:
        print(prec, file=stream)
        for f in (x,y,z):
            print(*(ZZ(f[f.valuation()+i]) for i in range(prec)), file=stream)
    print('Verified the Fricke equation; wrote', args.units/'series.dat')


if __name__ == '__main__':
    main()
