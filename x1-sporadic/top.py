"""Single-prime density bound from the full self-intersection of the incidence divisor."""
from fractions import Fraction
from functools import cache

import sympy as sp

from joint import X, Y, _raw, _trace_polynomials


@cache
def _correction(label, dimension, p, elliptic_trace):
    """Return the contribution beyond p^2+1 in the density denominator."""
    if dimension == 1:
        return Fraction(max(0, 2*p - elliptic_trace**2))

    data = _raw().get(label)
    if data is None or data['cyclotomic'] or 'a'+str(p) not in data:
        # Every conjugate of a_p^2/epsilon(p) - 2p is at least -2p.
        return Fraction(2*p, dimension)

    field, traces = _trace_polynomials(label)
    assert field.degree() == dimension
    t = traces[p]
    if t.degree() <= 0:
        return max(Fraction(0), -Fraction(t.nth(0))/dimension)

    poly = sp.Poly(sp.resultant(field.as_expr(), Y-t.as_expr(), X), Y).monic()
    assert poly.degree() == dimension
    # All roots are real. Strictly positive coefficients are equivalent to
    # strict negativity of all roots (multiplicities included).
    if not all(c > 0 for c in poly.all_coeffs()):
        return Fraction(0)
    # P'(0)/P(0) = sum_sigma ( -1 / sigma(t) ).
    return Fraction(poly.nth(0))/Fraction(poly.nth(1))


def top_bound(form, p, S):
    """Return the exact rational single-prime lower bound on the moving degree."""
    if form['level'] % p == 0:
        raise ValueError('The prime must be prime to the level of the newform')
    dimension = form['dim']
    ap = form['initial_traces'][p-1] if dimension == 1 else None
    denominator = p*p + 1 + _correction(form['label'], dimension, p, ap)
    return Fraction(S)/denominator
