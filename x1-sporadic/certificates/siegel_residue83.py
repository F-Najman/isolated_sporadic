"""sage -python siegel_residue83.py (or Python with cypari2 and SymPy)."""
from fractions import Fraction
from cypari2 import Pari
from siegel_contacts import cusp_contacts

pari = Pari()
pari.allocatemem(100_000_000, silent=True)
pari('''unitseries(e, N, M) = {
    my(z=1+O(t^(M+1)));
    for(n=1, M, if(e[n%N+1], z *= (1-t^n+O(t^(M+1)))^e[n%N+1]));
    z
}''')

pari('''logcoefficients(u, v, M) = {
    my(w, g);
    if(v < 0, u = u-t^(-v), u = t^v*u-1);
    w = valuation(u,t);
    g = min(v,0)+w + t*deriv(u/t^w,t)/(u/t^w);
    vector(M+1, i, polcoef(g,i-1,t))
}''')

pari('''residues_at_one_cusp(u, g, pole, sign, D) = {
    my(z=1+O(t^(pole*D+1)), r=vector(D));
    for(k=1, D, z*=u*sign; r[k]=polcoef(z*g,pole*k,t));
    r
}''')


def expansion_data(e, c, N):
    # Representatives with c odd and bottom-right entry 2N have no
    # nontrivial cyclotomic coefficient in the positive-degree factors.
    if c % 2 == 0:
        c += N
    period = [0]*N
    val, parity = 0, 0
    for i, exponent in e.items():
        shift, s = divmod(c*i, N)
        period[s] += exponent
        period[-s % N] += exponent
        val += exponent*s*(s-N)
        parity += exponent*shift
    assert val % (2*N) == 0
    return val//(2*N), (-1)**(parity % 2), period


N = 83
f = {35:1, 40:1, 3:-1, 34:-1}
h = {13:1, 3:1, 6:-1, 15:-1}
degree, hits, _ = cusp_contacts(N, f, (1,))
contact = {c:k for c,d,k in hits[1]}
assert degree == 95 and contact == {10:8, 18:27, 26:3}
D = degree-sum(contact.values())
sums = [Fraction(0)]*D
for c in range(1, (N+1)//2):
    v, fsign, ep = expansion_data(f, c, N)
    a, hsign, eh = expansion_data(h, c, N)
    w = contact.get(c, min(v, 0))
    if a > 0:
        continue
    if a == 0:
        for k in range(1, D+1):
            sums[k-1] -= w*hsign**k
        continue
    precision = -a*D
    f_unit = fsign*pari('unitseries')(ep, N, precision+contact.get(c,0)+2)
    log_vector = pari('logcoefficients')(f_unit, v, precision)
    log_series = log_vector.Polrev('t').Ser('t', precision+1)
    h_unit = pari('unitseries')(eh, N, precision)
    residues = pari('residues_at_one_cusp')(h_unit, log_series, -a, hsign, D)
    for i, value in enumerate(residues):
        sums[i] -= Fraction(str(value))

assert sums[:5] == [6,680,-9324,293520,-5571254]
coefficients = [Fraction(1)]
for j in range(1,D+1):
    coefficients.append(-sum(coefficients[k]*sums[j-k-1] for k in range(j))/j)
assert all(v.denominator == 1 for v in coefficients)
polynomial = pari.Pol([v.numerator for v in coefficients])
factors = polynomial.factormod(11,1)
assert list(map(int,factors[0])) == [57] and list(map(int,factors[1])) == [1]
print('X1(83): the noncuspidal divisor is a single degree-57 closed point.')
