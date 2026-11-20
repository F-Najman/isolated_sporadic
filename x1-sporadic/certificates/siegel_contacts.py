"""Exact cusp contacts of products of g_(0,i/N) on X_1(N), for N >= 5."""
from functools import cache
from math import gcd
from collections import defaultdict
import sympy as sp


@cache
def _modulus(M):
    x = sp.Symbol('x')
    return tuple(reversed([int(v) for v in sp.Poly(sp.cyclotomic_poly(M, x), x).all_coeffs()]))


def _zero_at_zeta(M, coeff):
    """Check a rational integral relation among powers of a primitive M-th root."""
    A = [0]*M
    for k,v in coeff.items():
        A[k % M] += v
    F = _modulus(M)
    s = len(F) - 1
    for i in range(M-1, s-1, -1):
        t = A[i]
        if t:
            for j,a in enumerate(F):
                A[i-s+j] -= t*a
    return not any(A[:s])


def _factor(A, k, M):
    B = A.copy()
    for j,t in A.items():
        y = (j+k)%M
        B[y] = B.get(y,0)-t
    return {j:t for j,t in B.items() if t}


@cache
def cusp_representatives(N):
    """Bottom rows of scaling matrices: c modulo N and d modulo gcd(c,N)."""
    reps = []
    for c0 in range(N//2+1):
        c = c0 or N
        h = gcd(c,N)
        for residue in range(h):
            if gcd(residue,h) != 1:
                continue
            if 2*c0 % N == 0 and residue > (-residue)%h:
                continue
            d = residue
            while gcd(c,d) != 1:
                d += h
            reps.append((c,d,h))
    return tuple(reps)


def _leading(N,c,d,e,s2):
    """Value of the leading coefficient in Q(zeta_(2N))."""
    M = 2*N
    phase = c*d*(s2//N)
    num,den = {0:1},{0:1}
    for i,t in e.items():
        k,s = divmod(c*i,N)
        phase += t*((k+1)*N-d*i*(1+2*k))
        if s == 0:
            exponent = 2*d*i % M
            assert exponent
            for _ in range(abs(t)):
                if t > 0:
                    num = _factor(num,exponent,M)
                else:
                    den = _factor(den,exponent,M)
    return {(j+phase)%M:t for j,t in num.items()},den


def _first_log_term(N,c,d,h,e,degree):
    """Order in the actual uniformizer q^(h/N) after dividing out the leading term."""
    terms = defaultdict(list)
    for i,t in e.items():
        s, b = c*i % N, d*i % N
        terms[s].append((b,t))
        terms[(-s)%N].append((-b,t))
    for contact in range(1,degree+1):
        n = h*contact
        coeff = defaultdict(int)
        for r in sp.divisors(n):
            for b,t in terms.get(r%N,()):
                coeff[b*(n//r)%N] += r*t
        if not _zero_at_zeta(N,coeff):
            return contact
    raise ValueError('The full degree bound contains no nonzero logarithmic coefficient.')


def cusp_contacts(N, exponents, targets=(1,-1)):
    """Return degree, lists of (c,d,contact) over each target, and cusp orders."""
    e = {i:t for i,t in exponents.items() if t}
    if N < 5 or not e or min(e) <= 0 or 2*max(e)>N:
        raise ValueError('Use a nonempty exponent map on 1 <= i <= N/2, N >= 5.')
    s2 = sum(i*i*t for i,t in e.items())
    assert sum(e.values())%12==0 and s2%(N*gcd(2,N))==0
    assert set(targets) <= {1,-1}
    cusp_data=[]
    for c,d,h in cusp_representatives(N):
        # ord_{q^(h/N)} g_(0,i/N) = (6s²-6Ns+N²)/(12Nh).
        b=sum(t*(6*(c*i%N)**2-6*N*(c*i%N)+N*N) for i,t in e.items())
        v,rem=divmod(b,12*N*h)
        assert not rem
        cusp_data.append((c,d,h,v))
    degree = sum(max(-v,0) for _,_,_,v in cusp_data)
    assert degree==sum(max(v,0) for _,_,_,v in cusp_data)
    if not degree:
        raise ValueError('The proposed Siegel unit is constant.')
    hits={target:[] for target in targets}
    for c,d,h,v in cusp_data:
        if v:
            continue
        A,B=_leading(N,c,d,e,s2)
        for target in targets:
            delta=A.copy()
            for j,t in B.items():
                delta[j]=delta.get(j,0)-target*t
            if _zero_at_zeta(2*N,delta):
                k=_first_log_term(N,c,d,h,e,degree)
                hits[target].append((c,d,k))
    assert all(sum(k for _,_,k in v)<=degree for v in hits.values())
    return degree,hits,[(c,d,v) for c,d,_,v in cusp_data]


def cusp_fibres(N, exponents, targets=(1,-1)):
    degree,hits,_=cusp_contacts(N,exponents,targets)
    return degree,{t:sum(k for _,_,k in v) for t,v in hits.items()}


if __name__ == '__main__':
    U = {12:1, 27:1, 6:-1, 21:-1}
    degree,hits,_ = cusp_contacts(99, U, (1,))
    assert degree == 114 and len(hits[1]) == 14
    assert sorted(k for _,_,k in hits[1]) == [3]*10 + [6]*2 + [9]*2
    residual = degree - sum(k for _,_,k in hits[1])
    assert residual == 54
    print('X1(99): degree', degree, 'cuspidal degree', degree-residual,
          'noncuspidal degree', residual)

    degree, cusps = cusp_fibres(138, {1:1, 3:1, 36:1, 2:-1, 9:-1, 51:-1}, (1,))
    assert degree == 185 and cusps[1] == 103
    print('X1(138): noncuspidal degree at most', degree - cusps[1])
