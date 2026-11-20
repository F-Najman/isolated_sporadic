"""Exact residue fields of two fibres on X_1(138). Run with sage -python."""
from collections import Counter
from fractions import Fraction
from math import gcd
from sys import argv

from cypari2 import Pari
from siegel_contacts import _leading, cusp_contacts

N = 138
h = {14: 1, 37: 1, 32: -1, 55: -1}
pari = Pari()
pari.allocatemem(250_000_000, silent=True)
cyclotomic = pari.polcyclo(2*N, 'x')

pari('''unitseries(e, N, M) = {
    my(z = 1+O(t^(M+1)));
    for(n=1, M, if(e[n%N+1], z *= (1-t^n+O(t^(M+1)))^e[n%N+1]));
    z
}''')
pari('''poleresidues(ef, eh, vf, pole, sf, sh, contact, D, N) = {
    my(M=pole*D, f=sf*unitseries(ef,N,M+contact+2), H=sh*unitseries(eh,N,M));
    my(z=1+O(t^(M+1)), r=vector(D), w, logder);
    if(vf<0, f=f-t^(-vf), f=t^vf*f-1);
    w=valuation(f,t);
    f=f/t^w;
    logder=min(vf,0)+w+t*deriv(f,t)/f;
    for(k=1, D, z*=H; r[k]=polcoef(z*logder,pole*k,t));
    r
}''')


def primitive_expansion(unit, c):
    # The bottom row (c,276) is primitive. Its positive-degree factors
    # have integral coefficients in the uniformizer t=q^(1/138).
    exponents = [0]*N
    valuation = parity = 0
    for i, e in unit.items():
        shift, s = divmod(c*i, N)
        exponents[s] += e
        exponents[-s % N] += e
        valuation += e*s*(s-N)
        parity += e*shift
    assert valuation % (2*N) == 0
    return valuation//(2*N), (-1)**(parity % 2), exponents


def leading_value(c, d):
    numerator, denominator = _leading(N, c, d, h, sum(e*i*i for i, e in h.items()))

    def evaluate(terms):
        coefficients = [terms.get(i, 0) for i in range(2*N)]
        return pari.Mod(pari.Polrev(coefficients, 'x'), cyclotomic)

    return evaluate(numerator) / evaluate(denominator)


def residual_polynomial(f, expected_degree, expected_residual, prime):
    degree, hits, valuations = cusp_contacts(N, f, (1,))
    contacts = {(c, d): order for c, d, order in hits[1]}
    D = degree - sum(contacts.values())
    assert (degree, D) == (expected_degree, expected_residual)
    _, _, h_orders = cusp_contacts(N, h, ())
    h_orders = {(c, d): v for c, d, v in h_orders}
    sums = [Fraction(0)]*D
    regular_values = Counter()
    for c, d, vf in valuations:
        contact = contacts.get((c, d), 0)
        order = contact or min(vf, 0)
        vh = h_orders[c, d]
        if gcd(c, N) > 1:
            assert vh == 0
            if order:
                regular_values[str(leading_value(c, d).minpoly())] += order
            continue
        v, sf, ef = primitive_expansion(f, c)
        a, sh, eh = primitive_expansion(h, c)
        assert (v, a) == (vf, vh)
        if a > 0:
            continue
        if a == 0:
            for k in range(1, D+1):
                sums[k-1] -= order*sh**k
            continue
        residues = pari('poleresidues')(ef, eh, vf, -a, sf, sh, contact, D, N)
        for i, value in enumerate(residues):
            sums[i] -= Fraction(str(value))
    # The sum of cusp residues is rational. Taking the normalized field
    # trace replaces each regular cusp value by the average of its conjugates.
    for polynomial, weight in regular_values.items():
        polynomial = pari(polynomial)
        powers = polynomial.polsym(D)
        for k in range(1, D+1):
            sums[k-1] -= Fraction(weight, int(polynomial.poldegree()))*Fraction(str(powers[k]))
    coefficients = [Fraction(1)]
    for j in range(1, D+1):
        coefficients.append(-sum(coefficients[k]*sums[j-k-1] for k in range(j))/j)
    assert all(v.denominator == 1 for v in coefficients)
    polynomial = pari.Pol([v.numerator for v in coefficients], 'x')
    factors = polynomial.factormod(prime, 1)
    assert list(map(int, factors[0])) == [D] and list(map(int, factors[1])) == [1]
    print('degree, cusp contribution, residual, irreducibility prime, first four power sums',
          degree, degree-D, D, prime, [str(x) for x in sums[:4]], flush=True)
    if '--polynomials' in argv:
        print(polynomial, flush=True)
    return polynomial


if __name__ == '__main__':
    residual_polynomial({1: 1, 3: 1, 36: 1, 2: -1, 9: -1, 51: -1}, 185, 82, 43)
    residual_polynomial({1: 1, 2: 1, 56: 1, 4: -1, 8: -1, 41: -1}, 203, 86, 71)
