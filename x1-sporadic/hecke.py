"""Exact extrema of a_p^2 / epsilon(p) for weight-two newforms."""
from collections import deque
from fractions import Fraction as Q
import json, math
import sympy as sp
x=sp.symbols('x')
primes=(2,3,5,7,11,13,17,19)
def values(record, extra):
    n=record['level']
    if record['dim']==1:
        assert record['char_order']==1
        return {p:[Q(record['initial_traces'][p-1]**2)]*3 for p in primes if n%p}
    form=extra[record['label']]
    assert form['dim']==record['dim']
    assert not form['cyclotomic']
    f=sp.Poly.from_list(form['field_poly'][::-1],x,domain='QQ')
    basis=[sp.Poly(x**i,x,domain='QQ') for i in range(f.degree())]
    if not form['power_basis']:
        basis=[sp.Poly.from_list(row[::-1],x,domain='QQ')/d for row,d in zip(form['nums'],form['dens'])]
        basis=[sp.Poly(b,x,domain='QQ') for b in basis]
    def poly(coeffs):
        return sum((c*b for c,b in zip(coeffs,basis)),sp.Poly(0,x,domain='QQ')) % f
    one=sp.Poly(1,x,domain='QQ')
    data=form['char_values']
    if data is None:
        assert record['char_order']==1
        chi={a:one for a in range(n) if math.gcd(a,n)==1}
    else:
        generators=[(a,poly(c)) for a,c in data]
        q=deque([1]);chi={1:one}
        while q:
            a=q.popleft()
            for b,u in generators:
                y=a*b%n
                w=(chi[a]*u)%f
                if y not in chi:
                    chi[y]=w;q.append(y)
                else:
                    assert chi[y]==w
        assert len(chi)==sp.totient(n)
    out={}
    for p in primes:
        if n%p==0:continue
        a=poly(form['a'+str(p)])
        # Check the field encoding against the independent absolute Hecke trace.
        ap_char=sp.Poly(sp.resultant(f.as_expr(),sp.Symbol('Y')-a.as_expr(),x),sp.Symbol('Y'))
        assert -ap_char.nth(ap_char.degree()-1)==record['initial_traces'][p-1]
        v=(a*a*chi[p%n].invert(f))%f
        h=sp.resultant(f.as_expr(),sp.Symbol('Y')-v.as_expr(),x)
        P=sp.Poly(h,sp.Symbol('Y'),domain='QQ').monic()
        intervals=P.intervals(eps=sp.Rational(1,10**8))
        assert sum(k for _,k in intervals)==record['dim']
        lower=Q(intervals[0][0][0]);upper_at_least=Q(intervals[-1][0][0])
        mean=Q(-P.nth(P.degree()-1))/P.degree()
        assert 0<=lower<=mean<=4*p and 0<=upper_at_least<=4*p
        out[p]=[lower,mean,upper_at_least]
    return out


def eigenvalue_lower_bounds(records, data):
    extra={f['label']:f for f in json.loads((data/'hecke_spectral.json').read_text())}
    return {r['label']:{p:v[2] for p,v in values(r,extra).items()} for r in records
            if not (r['analytic_rank_proved'] and r['analytic_rank']==0)}
