"""Write the Fricke model of X_0(131) and its two j-maps as Magma functions.

Run with Sage's Python. The output path is the sole command-line argument.
"""
import sys
from pathlib import Path

from sage.all import EllipticCurve, LaurentSeriesRing, O, QQ, ZZ
from sage.modular.modform.j_invariant import j_invariant_qexp


def in_elliptic_coordinates(series, pole_order, x, y):
    # A basis of L(m O) uses each pole order other than 1 exactly once.
    terms = []
    remainder = series
    for order in list(range(pole_order, 1, -1)) + [0]:
        j = order % 2
        i = (order - 3 * j) // 2
        monomial = x**i * y**j
        coefficient = ZZ(remainder[-order] / monomial[-order])
        if coefficient:
            terms.append(f"({coefficient})*x^{i}*y^{j}")
            remainder -= coefficient * monomial
    assert remainder.valuation() > 15
    return "+".join(terms)


def main(output):
    precision = 305
    ring = LaurentSeriesRing(QQ, "q", default_prec=precision)
    q = ring.gen()
    curve = EllipticCurve("131a1")
    assert curve.a_invariants() == (0, -1, 1, 1, 0)
    x, y = [ring(t) for t in curve.modular_parametrization().power_series(precision + 10)]
    eta_product = q**11 * (1 + O(q**(precision + 1)))
    for k in range(1, precision):
        eta_product *= (1 - q**k)**2
        if 131 * k < precision:
            eta_product *= (1 - q**(131 * k))**2
    z = ring(curve.q_eigenform(precision + 10)) / eta_product
    j = ring(j_invariant_qexp(precision))
    j_dual = q**-131 + 744 + 196884 * q**131 + O(q**262)
    forms = [("R", z**2, 20), ("Jp", j + j_dual, 131), ("Jm", (j - j_dual) / z, 121)]
    lines = []
    for name, series, order in forms:
        polynomial = in_elliptic_coordinates(series, order, x, y)
        lines.append(f"function {name}(x,y) return {polynomial}; end function;")
    Path(output).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
