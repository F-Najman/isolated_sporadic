"""Exact level-123 certificate. Run with sage -python siegel_quotient123.py.
By default verify both the quotient and the exact degree-48 point upstairs.
Pass --full to also replay the pullback of the quotient coordinate.
"""
from collections import Counter
from fractions import Fraction
from math import gcd
import sys

from cypari2 import Pari
from sympy import totient
from siegel_contacts import cusp_contacts, _leading

N = 123
GROUP = (1, 14, 32, 38, 40, 44, 50, 55)
F = {4:1, 5:1, 26:1, 29:1, 37:1, 46:1, 53:1, 56:1,
     2:-1, 13:-1, 23:-1, 28:-1, 35:-1, 43:-1, 47:-1, 59:-1}
H = {6:2, 8:1, 10:1, 11:1, 17:1, 18:2, 31:1, 39:2, 49:1, 52:1, 54:2, 58:1,
     2:-1, 13:-1, 23:-1, 24:-2, 28:-1, 30:-2, 33:-2, 35:-1, 43:-1, 47:-1,
     51:-2, 59:-1}
pari = Pari()
pari.allocatemem(400_000_000, silent=True)
pari('''unitseries(e,N,M) = {
  my(z=1+O(t^(M+1)));
  for(n=1,M, if(e[n%N+1], z *= (1-t^n+O(t^(M+1)))^e[n%N+1])); z
}''')
pari('''localres(u,v,h,pole,sg,D) = {
  my(w,g,M=pole*D,z=1+O(t^(M+1)),r=vector(D));
  if(v<0, u=u-t^(-v), u=t^v*u-1);
  w=valuation(u,t); g=min(v,0)+w+t*deriv(u/t^w,t)/(u/t^w);
  for(k=1,D,z*=sg*h;r[k]=polcoef(z*g,pole*k,t)); r
}''')


def folded(i):
    return min(i % N, (-i) % N)


def check_invariance(e):
    assert sum(e.values()) == 0 and sum(i*i*a for i,a in e.items()) % N == 0
    for a in GROUP:
        if a % 2 == 0:
            a += N
        assert {folded(a*i):v for i,v in e.items()} == e
        assert sum(v*(a*i//N) for i,v in e.items()) % 2 == 0


def data_at_rational_cusp(e, c):
    if c % 2 == 0:
        c += N
    period = [0]*N
    order = parity = 0
    for i,v in e.items():
        q,r = divmod(i*c,N)
        period[r] += v
        period[-r % N] += v
        order += v*r*(r-N)
        parity += v*q
    assert order % (2*N) == 0
    return order//(2*N), (-1)**(parity % 2), period


def closed_cusp_sums(H, c, weight, D):
    # All geometric values are conjugate, and each conjugate occurs equally often.
    assert sum(v*(c*i%N)*(c*i%N-N) for i,v in H.items()) == 0
    num,den = _leading(N,c,1,H,sum(i*i*v for i,v in H.items()))
    mod = pari.polcyclo(2*N,'x')
    def elt(row):
        return pari.Mod(pari.Polrev([row.get(i,0) for i in range(2*N)]),mod)
    alpha = elt(num)/elt(den)
    pol = alpha.minpoly()
    g = int(pol.poldegree())
    assert weight % g == 0
    coeffs = [Fraction(str(v)) for v in reversed(pol.Vecrev())]
    tr = [Fraction(g)]
    for k in range(1,D+1):
        t = -sum(coeffs[j]*tr[k-j] for j in range(1,min(k-1,g)+1))
        if k <= g:
            t -= k*coeffs[k]
        tr.append(t)
    return [v*(weight//g) for v in tr[1:]]


assert {folded(a*b) for a in GROUP for b in GROUP} == set(GROUP)
for unit in (F,H):
    check_invariance(unit)
degree,all_hits,all_orders = cusp_contacts(N,F,(1,))
hits = {(a,b):r for a,b,r in all_hits[1]}
assert degree == 256 and Counter(hits.values()) == {1:80,2:16}
residual = degree-sum(hits.values())
assert residual == 144
def characteristic_polynomial(H, group_size):
    D = residual//group_size
    power_sums = [Fraction(0)]*D
    for c in range(N//2+1):
        cc = c or N
        orders = {v for a,b,v in all_orders if a == cc}
        assert len(orders) == 1
        v = orders.pop()
        tails = {hits.get((a,b),min(v,0)) for a,b,w in all_orders if a == cc}
        assert len(tails) == 1
        tail = tails.pop()
        if gcd(c,N)>1:
            if tail:
                weight = int(totient(gcd(c,N)))//(2 if c==0 else 1)
                local = [tail*q for q in closed_cusp_sums(H,cc,weight,D)]
            else:
                local = [0]*D
        else:
            vf,sf,ef = data_at_rational_cusp(F,c)
            vh,sh,eh = data_at_rational_cusp(H,c)
            assert v == vf
            if vh > 0:
                local = [0]*D
            elif vh == 0:
                local = [tail*sh**k for k in range(1,D+1)]
            else:
                precision = -vh*D
                u = sf*pari('unitseries')(ef,N,precision+max(0,tail)+2)
                h = pari('unitseries')(eh,N,precision)
                local = [Fraction(str(x)) for x in pari('localres')(u,vf,h,-vh,sh,D)]
        for j,x in enumerate(local):
            power_sums[j] -= x
    sums = [v/group_size for v in power_sums]
    co = [Fraction(1)]
    for j in range(1,D+1):
        co.append(-sum(co[k]*sums[j-k-1] for k in range(j))/j)
    pol = pari.Pol([pari(str(c)) for c in co])
    return pol

target = pari('x^6-399*x^5+426*x^4-911*x^3+300*x^2-66*x+1')
assert target.polisirreducible()
def factors(pol):
    fac = pol.factor()
    return fac, [(int(x.poldegree()),int(m)) for x,m in zip(fac[0],fac[1])]

quotient, pattern = factors(characteristic_polynomial(H, len(GROUP)))
assert pattern == [(6,1),(12,1)] and quotient[0][0] == target
print("degree-18 quotient divisor:",pattern)

# The ratio g_16/g_25 is a rational modular unit. Its primitive-cusp
# poles have order at most 13 and at every other cusp it has order zero.
coordinate = {16:1,25:-1}
assert (16**2-25**2) % N == 0
assert all(sum(v*(c*i%N)*(c*i%N-N) for i,v in coordinate.items()) == 0
           for c in range(N//2+1) if gcd(c,N)>1)
upstairs, pattern = factors(characteristic_polynomial(coordinate, 1))
assert pattern == [(48,1),(96,1)]
for factor in upstairs[0]:
    assert factor.polisirreducible()
print("degree-144 noncuspidal divisor:",pattern)
print("defining polynomial of the degree-48 point:",upstairs[0][0])

if '--full' in sys.argv:
    _, pattern = factors(characteristic_polynomial(H, 1))
    assert pattern == [(6,8),(12,8)]
    print("pullback of the quotient coordinate:",pattern)
