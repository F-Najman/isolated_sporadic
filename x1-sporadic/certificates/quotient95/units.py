"""Split-cusp unit lattice of X_1(131)/H_5 and ternary norms to X_0(131)."""
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
import numpy as np
from sympy import Matrix, eye, ones

HERE = Path(__file__).resolve().parent
LEVEL = 131
DEGREE = 95
PRECISION = 160


def fold(a):
    return min(a % LEVEL, -a % LEVEL)


group = sorted({fold(2**(13*i)) for i in range(5)})
orbits = [sorted(fold(2**j*h) for h in group) for j in range(13)]
labels = [min(orbit) for orbit in orbits]
assert len(set(sum(orbits, []))) == 65
orders = np.array([[sum((a*i % LEVEL)*((a*i % LEVEL)-LEVEL) for i in orbit)
                    for orbit in orbits] for a in labels])
diff = orders - orders[:, -1, None]
assert not np.any(diff % (2*LEVEL))
valuation_matrix = diff // (2*LEVEL)


def cyclic_minimum(e):
    return min(e[j:]+e[:j] for j in range(13))


def all_units():
    """Enumerate the entire lattice ball, up to a diamond and inversion."""
    square = Matrix(valuation_matrix[:-1, :-1].tolist())
    assert abs(square.det()) == 5548598204321224637
    inverse = square.inv()
    inverse = inverse.col_join(-ones(1, 12)*inverse)
    # A degree-d unit divisor is in the convex hull of d*(e_i-e_j).
    identity = eye(13)
    vertices = [DEGREE*inverse*(identity[:12,i]-identity[:12,j])
                for i in range(13) for j in range(i)]
    assert max(abs(a) for v in vertices for a in v) < 3
    assert max(sum(abs(a) for a in v) for v in vertices) < 14
    assert max(sum(a*a for a in v) for v in vertices) < 17

    def halves(length):
        groups = defaultdict(list)
        for v in product(range(-2, 3), repeat=length):
            l1 = sum(map(abs, v))
            square_norm = sum(t*t for t in v)
            if l1 <= 12 and square_norm <= 16:
                groups[sum(v), l1, square_norm].append(v)
        return {key:np.array(value, dtype=np.int16) for key,value in groups.items()}

    first, second = halves(6), halves(7)
    right_vals = {key:(v @ valuation_matrix[:,6:].T).astype(np.int16)
                  for key,v in second.items()}
    found = set()
    tested = 0
    for (s,l,q),left in first.items():
        left_vals = (left @ valuation_matrix[:,:6].T).astype(np.int16)
        for key,right in second.items():
            if s+key[0] or l+key[1] > 12 or q+key[2] > 16:
                continue
            tested += len(left)*len(right)
            for e,v in zip(left,left_vals):
                for i in np.where(np.abs(v+right_vals[key]).sum(axis=1) <= 2*DEGREE)[0]:
                    exp = tuple(map(int,e)) + tuple(map(int,right[i]))
                    if any(exp):
                        found.add(min(cyclic_minimum(exp),cyclic_minimum(tuple(-t for t in exp))))
    assert tested == 8004595 and len(found) == 18
    return sorted(found)


def one(n):
    result = np.zeros(n, dtype=np.int64)
    result[0] = 1
    return result


def multiply(a,b,n,modulus=3):
    return np.convolve(a[:n],b[:n])[:n] % modulus


def reciprocal(a,modulus=3):
    n = len(a)
    b = one(n)
    b[0] = a[0]
    for i in range(1,n):
        b[i] = (-b[0]*(a[1:i+1] @ b[i-1::-1])) % modulus
    return b


class Norms:
    def __init__(self,modulus=3):
        self.modulus=modulus
        rows = (HERE/"fricke-series.txt").read_text().splitlines()
        x,y,z = (np.array([int(a)%self.modulus for a in row.split()][:PRECISION+21])
                 for row in rows[1:])
        self.basis = [None]*(PRECISION+1)
        self.basis[0] = one(PRECISION+1)
        self.basis[2], self.basis[3] = x,y
        for i in range(4,PRECISION+1):
            self.basis[i] = multiply(self.basis[i-2],x,PRECISION+1,self.modulus)
        self.anti = [None if v is None else multiply(z,v,PRECISION+1,self.modulus) for v in self.basis]

        self.signs = np.empty((13,13),dtype=int)
        self.local = {}
        self.powers = {}
        for i,label in enumerate(labels):
            a = label if label%2 else label+LEVEL
            for j,orbit in enumerate(orbits):
                self.signs[i,j] = sum(a*h//LEVEL for h in orbit)%2
                series = one(PRECISION+1)
                for h in orbit:
                    r = a*h%LEVEL
                    for start in (r,LEVEL-r):
                        for n in range(start,PRECISION+1,LEVEL):
                            series[n:] = (series[n:]-series[:-n])%self.modulus
                self.local[i,j] = series

    def local_expansion(self,exp,i):
        unit = one(PRECISION+1)
        for j,power in enumerate(exp):
            if not power:
                continue
            key = (i,j,power)
            if key not in self.powers:
                factor = self.local[i,j]
                if power < 0:
                    factor = reciprocal(factor,self.modulus)
                value = one(PRECISION+1)
                for _ in range(abs(power)):
                    value = multiply(value,factor,PRECISION+1,self.modulus)
                self.powers[key] = value
            unit = multiply(unit,self.powers[key],PRECISION+1,self.modulus)
        sign = (-1)**int((self.signs[i]@exp)%2)
        return int(valuation_matrix[i]@exp), sign*unit % self.modulus

    def signed_fibre(self,exp,epsilon):
        norm = one(PRECISION+1)
        cusp_orders, contacts = [], []
        for i in range(13):
            v,a = self.local_expansion(exp,i)
            first = min(v,0)
            shift = v-first
            expansion = np.zeros(PRECISION+1,dtype=np.int64)
            expansion[shift:] = a[:PRECISION+1-shift]
            expansion[-first] -= epsilon
            expansion %= self.modulus
            contact = int(np.flatnonzero(expansion)[0])
            cusp_orders.append(first+contact)
            contacts.append(contact)
            normalized = np.pad(expansion[contact:],(0,contact))
            norm = multiply(norm,normalized,PRECISION+1,self.modulus)
        d = -sum(cusp_orders)
        assert 0 < d <= DEGREE and max(contacts)+d <= PRECISION
        return self.recover(norm[:d+1]),cusp_orders

    def recover(self,norm):
        """Write twice the norm as T(x,y) plus or minus z*W(x,y)."""
        coeff = np.asarray(norm,dtype=np.int64).copy()
        d = len(coeff)-1
        rem = coeff.copy()
        w,t = np.zeros(d+1,dtype=int),np.zeros(d+1,dtype=int)
        for j in range(d,0,-1):
            if j<10 or j==11:
                assert rem[d-j]==0
                continue
            basis = self.anti[j-10][:j+1]
            lead = int(rem[d-j]*basis[0])%self.modulus
            w[j-10] = lead
            rem[d-j:] = (rem[d-j:]-lead*basis)%self.modulus
        # The Fricke conjugate is regular at this cusp.
        coeff[-1] = (coeff[-1]+rem[-1])%self.modulus
        for j in range(d,-1,-1):
            if j==1:
                assert coeff[d-j]==0
                continue
            basis = self.basis[j][:j+1]
            lead = int(coeff[d-j]*basis[0])%self.modulus
            t[j] = lead
            coeff[d-j:] = (coeff[d-j:]-lead*basis)%self.modulus
        assert not np.any(coeff)
        return (d,*map(int,t),*map(int,w))


def write_rows(path,rows):
    path.write_text("".join(" ".join(map(str,row))+"\n" for row in rows))


def prepare(directory,engine):
    exponents = all_units()
    translates = []
    fibres = []
    for i,e in enumerate(exponents):
        v = valuation_matrix @ e
        for signed in (v,-v):
            translates.extend(tuple(np.roll(signed,j)) for j in range(13))
        for j,epsilon in enumerate((1,-1)):
            poly,_ = engine.signed_fibre(e,epsilon)
            fibres.append((2*i+j,*poly))
    assert len(set(translates)) == 468
    write_rows(directory/"unit-divisors.dat",translates)
    write_rows(directory/"unit-norms.dat",fibres)
    return dict(sorted(Counter(sum(abs(x) for x in valuation_matrix@e)//2 for e in exponents).items()))
