"""sage -python canonical117.py. Count supersingular zeros of the 117a differential."""
import numpy as np
from sage.all import CuspForms, CyclotomicField, EllipticCurve, GF, ModularSymbols, QQ, ZZ, divisors, matrix, vector

N = 117
precision = 75
S = CuspForms(N, 2)
forms = list(S.q_integral_basis(precision))
f = EllipticCurve('117a1').q_eigenform(precision)

# Good reduction and the q-expansion principle identify the reduced integral
# basis with the regular differentials of X0(117) at 2.
F2 = GF(2)
B = matrix(F2, [[a[i] for i in range(precision)] for a in forms])
assert B.rank() == len(forms) == 11
f_coordinates = B.solve_left(vector(F2, [f[i] for i in range(precision)]))
pivot = next(j for j,v in enumerate(f_coordinates) if v)
forms = [f] + [a for i,a in enumerate(forms) if i != pivot]
pairs = [(i,j) for i in range(11) for j in range(i,11)]
products = [forms[i]*forms[j] for i,j in pairs]
products = matrix(F2, [[a[i] for i in range(precision)] for a in products])
assert products.rank() == 30
# The kernel gives quadratic equations. A nonzero quadratic differential has
# at most 4g-4=40 zeros, so the series precision suffices. Product rank 30
# excludes a hyperelliptic canonical image, which would have rank at most 21.
quadrics = [[pair for pair,value in zip(pairs,row) if value]
            for row in products.left_kernel_matrix()]

# Encode 0,1,z,z+1 in F4 by 0,1,2,3, respectively. Addition is bitwise XOR.
times = np.array([[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]], dtype=np.uint8)
points = zeros = 0
for pivot in range(11):
    size = 4**(10-pivot)
    for start in range(0,size,65536):
        ids = np.arange(start,min(start+65536,size),dtype=np.int32)
        coords = [np.zeros_like(ids)]*pivot + [np.ones_like(ids)]
        coords += [(ids >> (2*i)) & 3 for i in range(10-pivot)]
        coords = np.vstack(coords).astype(np.uint8)
        for quadric in quadrics:
            value = np.zeros(coords.shape[1],dtype=np.uint8)
            for i,j in quadric:
                value ^= times[coords[i],coords[j]]
            coords = coords[:,value == 0]
            if not coords.shape[1]:
                break
        points += coords.shape[1]
        zeros += int(np.count_nonzero(coords[0] == 0))
assert (points,zeros) == (22,6)

# Four rational cusps form the Atkin--Lehner orbit of infinity, where a1=1.
# The four others form the orbit of 1/3 under w13 and Gal(Q(zeta_3)/Q).
# The translation U: tau -> tau+1/3 normalizes Gamma0(117). Since U w117
# sends infinity to 1/3, test the first coefficient of (f|U)|w117.
K = CyclotomicField(3)
z = K.gen()
basis = []
fricke_first = []
for level in divisors(N):
    M = ModularSymbols(level,2,sign=1).cuspidal_submodule().new_submodule()
    for A in M.decomposition():
        if not A.dimension():
            continue
        W = A.atkin_lehner_operator(level).matrix()
        w = W[0,0]
        assert W == w and abs(w) == 1
        for a in A.q_expansion_basis(precision):
            for d in divisors(N//level):
                other = N//(level*d)
                basis.append([a[i//d] if i % d == 0 else 0 for i in range(1,precision)])
                # a(d tau)|w_N = w (other/d) a(other tau), for weight two.
                fricke_first.append(w*QQ(other)/d*a[1] if other == 1 else 0)
shifted = vector(K, [f[i]*z**i for i in range(1,precision)])
coeffs = matrix(K,basis).solve_left(shifted)
initial = coeffs.dot_product(vector(K,fricke_first))
assert initial == -(z+2)/3 and initial.norm() == QQ(1)/3

# All 8 geometric cusps are F4-rational; the 14 supersingular points are as
# well, since Frob4 acts as -2 on a scalar model. There are 168 cyclic
# 117-subgroups of the unique supersingular curve, and its 12 automorphisms
# modulo sign act freely (no roots of X^2+1 mod 3 or X^2+X+1 mod 9).
assert all((i*i+1) % 3 for i in range(3))
assert all((i*i+i+1) % 9 for i in range(9))
assert (117*(3+1)*(13+1))//(3*13) == 168 and 168//12+8 == points
print('X0(117)(F4) has 8 cusps and 14 supersingular points.')
print('The 117a differential vanishes at exactly six of them, all supersingular.')
