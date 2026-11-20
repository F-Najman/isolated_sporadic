"""Run with sage -python degree52_83.py.

Verify the supersingular fibres in characteristics 2 and 3, and test
invariant forms for the degree-eight Lang cover in characteristic 2.
"""
from itertools import combinations_with_replacement

from sage.all import EllipticCurve, GF, LaurentSeriesRing, O, QQ, ZZ, matrix
from sage.modular.modform.j_invariant import j_invariant_qexp


def root_of_norm(p, x, y):
    if p == 2:
        return x * (x + 1) * y
    return y + x**2 + x*y - x**3 + x**2*y


def check_norms():
    curve = EllipticCurve("83a1")
    assert list(curve.a_invariants()) == [1, 1, 1, 1, 0]
    series = LaurentSeriesRing(QQ, "q", default_prec=170)
    q = series.gen()
    x, y = map(series, curve.modular_parametrization().power_series(180))
    j = series(j_invariant_qexp(160))
    norm = j * (q**-83 + 744 + 196884*q**83 + O(q**166))
    for p in [2, 3]:
        difference = norm - root_of_norm(p, x, y)**12
        assert difference.precision_absolute() > 0
        assert all(ZZ(difference[i]) % p == 0 for i in range(-84, 1))
    return curve


def check_fibres(curve, p):
    field = GF(p**2, "a")
    E = curve.change_ring(field)
    want = (0, 1, 2, 0, 2, 1, 1, 0) if p == 2 else (
        0, 0, 0, 2, 2, 0, 2, 0, 2, 0, 0, 0, 2, 2, 2)
    assert E.cardinality() == len(want)
    zeroes = {}
    for P in E.points():
        if P.is_zero():
            continue
        x, y = P.xy()
        if root_of_norm(p, x, y) != 0:
            continue
        if p == 3:
            # Seven distinct zeroes saturate a pole of order seven.
            zeroes[P] = 2
        else:
            # Div(x(x+1)y) at finite points. The three roots of y are simple.
            zeroes[P] = int(x in (0, 1))*(1 + int(P == -P)) + int(y == 0)
    assert sum(zeroes.values()) == (7 if p == 2 else 14)
    for generator in E.points():
        if generator.order() != len(want):
            continue
        counts = tuple(zeroes.get(i*generator, 0) for i in range(len(want)))
        if counts == want:
            print(f"p={p}, cyclic generator {generator}, multiplicities {counts}")
            return counts
    raise AssertionError("Unexpected cyclic distribution")


def restriction_nullity(shifts, y_degree):
    # Invariant polynomials of bidegree (1,r) in x_0,...,x_7 and
    # y_0,...,y_7 have an orbit-sum basis indexed by multisets of r
    # subscripts of y relative to the one distinguished subscript of x.
    columns = list(combinations_with_replacement(range(8), y_degree))
    rows = {}
    for col, offsets in enumerate(columns):
        for m in shifts:
            for i in range(8):
                monomial = tuple(sorted([(i+m) % 8] + [(i+j) % 8 for j in offsets]))
                row = rows.setdefault((m, monomial), [0]*len(columns))
                row[col] ^= 1
    return len(columns) - matrix(GF(2), list(rows.values())).rank()


def main():
    curve = check_norms()
    counts = check_fibres(curve, 2)
    all_shifts = [i for i, count in enumerate(counts) if count]
    even_shifts = [i for i in all_shifts if i % 2 == 0]
    assert sum(counts[i] for i in even_shifts)*41 == 205
    results = [restriction_nullity(even_shifts, 1),
               restriction_nullity(all_shifts, 1),
               restriction_nullity(all_shifts, 2)]
    assert results == [1, 0, 0]
    print("Lang-cover nullities (bilinear even, bilinear all, cubic all):", results)
    check_fibres(curve, 3)
    k = 7
    excess = (10*52*k + 5*82) - 574*k
    assert excess == 32 and excess < 41
    print("Degree of the residual intersection at 3 for k=7:", excess)


if __name__ == "__main__":
    main()
