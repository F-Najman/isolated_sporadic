"""Run python3 quotient121.py (NumPy and a C++17 compiler are required)."""
from functools import cache
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from binary121 import field

# Base cusp order O,(0,0),(0,-1),(1,0),(1,-1). The first four columns give
# valuations on the elliptic base of x,y,x-1,x-y-1. Pullback multiplies by 11.
BASE_ORDERS = ((-2,-3,-2,-3), (1,2,0,0), (1,0,0,1), (0,1,1,2), (0,0,1,0))
T_ORDERS = (-19,5,4,9,1)
POLE = 50
DIRECTORY = Path(__file__).resolve().parent


def valuation(exponents, cusp):
    return (11 * sum(a*b for a,b in zip(exponents[:4], BASE_ORDERS[cusp]))
            + exponents[4] * T_ORDERS[cusp])


def compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for value in range(total + 1):
        for tail in compositions(total - value, length - 1):
            yield (value,) + tail


@cache
def monomials():
    # We only need 225 independent monomials in L(50R), where deg R = 5.
    # Distinct valuations at the first cusp prove linear independence.
    by_order = {}
    for s in range(11):
        bounds = [(POLE + s*v)//11 for v in T_ORDERS]
        for partition in compositions(sum(bounds), 5):
            _, a, b, c, d = [u-v for u,v in zip(partition,bounds)]
            j, remainder = divmod(c-d-2*b+2*a, 5)
            if remainder:
                continue
            i = a - 2*j
            powers = (i, j, d, b-i, s)
            order = valuation(powers, 0)
            if order not in by_order or sum(map(abs,powers)) < sum(map(abs,by_order[order])):
                by_order[order] = powers
    basis = tuple(by_order[v] for v in sorted(by_order))
    assert len(basis) == 5*POLE + 1-26 == 225
    assert all(valuation(z,c) >= -POLE for z in basis for c in range(5))
    return basis


@cache
def exponent_matrix():
    return np.array(monomials(), dtype=np.int32)


def affine_values(F, point, length):
    if length == 1:
        x,y,t = point
        logs = F.log[[x, y, x ^ 1, x ^ y ^ 1, t]]
        return F.exp[(exponent_matrix() @ logs) % F.period][None,:]
    series = F.affine_jets(point, length)

    @cache
    def factor(i, power):
        return F.series_power(series[i], power)

    @cache
    def product(powers):
        if len(powers) == 1:
            return factor(0, powers[0])
        return F.series_product(product(powers[:-1]), factor(len(powers)-1, powers[-1]))

    return np.array([product(z) for z in monomials()], dtype=np.int32).T


def cusp_values(cusp, length):
    # Let pi be a parameter on X_1(11) and z^11=pi. Trivialize L(50R)
    # by multiplying each function by z^50. For length <= 22 only the
    # first two coefficients of the unit part in pi can contribute.
    assert length <= 22
    F = field(1)
    precision = 12
    if cusp == 0:
        # pi=x/y, x=U/pi^2, y=U/pi^3, U^2+(1+pi^2)U+pi^3=0, U(0)=1.
        U = np.zeros(precision, dtype=np.int32)
        U[0] = 1
        for i in range(1, precision):
            U[i] = (U[i//2] if i%2 == 0 else 0) ^ (U[i-2] if i >= 2 else 0) ^ (i == 3)
        Xm = U.copy()
        Xm[2] ^= 1
        W = U.copy()
        W[1:] ^= U[:-1]
        W[3] ^= 1
        units = [U, U, Xm, W]
    else:
        x,y = [(0,0),(0,1),(1,0),(1,1)][cusp-1]
        X = np.zeros(precision, dtype=np.int32)
        X[:2] = [x,1]
        Xm = X.copy()
        Xm[0] ^= 1
        rhs = F.series_product(F.series_power(X,2), Xm)
        Y = np.zeros(precision, dtype=np.int32)
        Y[0] = y
        for i in range(1, precision):
            Y[i] = rhs[i] ^ (Y[i//2] if i%2 == 0 else 0)
        units = []
        for series, order in zip([X,Y,Xm,Xm^Y], BASE_ORDERS[cusp]):
            assert np.nonzero(series)[0][0] == order
            units.append(series[order:])
    units = [a[:2] for a in units]
    rhs = F.series_product(F.series_product(units[0], F.series_power(units[1],2)),
                           F.series_product(units[2], F.series_power(units[3],3)))
    units.append(F.eleventh_root_series(rhs, 1))
    result = np.zeros((length,len(monomials())), dtype=np.uint8)
    for j,z in enumerate(monomials()):
        order = POLE + valuation(z,cusp)
        if order < length:
            result[order,j] = 1
        if order + 11 < length:
            result[order+11,j] = sum(e*a[1] for e,a in zip(z,units)) % 2
    return result


def right_kernel(matrix,F):
    """Return a column basis for the null space over F."""
    reduced = matrix.copy()
    rows,cols = matrix.shape
    pivots = []
    for col in range(cols):
        i = len(pivots)
        candidates = np.flatnonzero(reduced[i:,col])
        if len(candidates) == 0:
            continue
        pivot = i + int(candidates[0])
        reduced[[i,pivot]] = reduced[[pivot,i]]
        reduced[i] = F.table[F.power(reduced[i,col],-1),reduced[i]]
        multipliers = reduced[:,col].copy()
        multipliers[i] = 0
        reduced ^= F.table[multipliers[:,None], reduced[i]]
        pivots.append(col)
        if len(pivots) == rows:
            break
    free = [j for j in range(cols) if j not in pivots]
    basis = np.zeros((cols,len(free)), dtype=np.int32)
    basis[free,:] = np.eye(len(free), dtype=np.int32)
    basis[pivots,:] = reduced[:len(pivots),free]
    for v in basis.T:
        assert not np.any(np.bitwise_xor.reduce(F.table[matrix, v[None,:]], axis=1))
    return basis


def twisted_spaces():
    # T=H(Phi)(sum of five translates of P - R); Phi is 2-Frobenius.
    # H(X) modulo X^7-1 is 16 + X + 12 X^2 + 8 X^4 + 3 X^6.
    F = field(7)
    assert F.modulus == 131
    coefficients = (16,1,12,0,8,0,3)
    points = set(F.affine_points())
    bases = []
    roots = []
    while points:
        root = min(points)
        orbit = {}
        current = root
        for j in range(7):
            rotated = current
            for _ in range(5):
                orbit[rotated] = j
                rotated = F.rotate(rotated)
            assert rotated == current
            current = F.frobenius(current)
        assert current == root and set(orbit) <= points and len(orbit) == 35
        points.difference_update(orbit)
        conditions = [affine_values(F,p,coefficients[j]) for p,j in orbit.items() if coefficients[j]]
        matrix = np.vstack(conditions)
        assert matrix.shape == (200,225)
        basis = right_kernel(matrix,F)
        assert basis.shape == (225,25)
        bases.append(basis)
        roots.append(root)
    assert roots == [(2,100,10),(2,101,91),(3,102,37),(3,103,16)]
    return np.array(bases, dtype=np.uint8)


def evaluation_data(path, bound):
    bases = twisted_spaces()
    twists,rows,columns = bases.shape
    matrix = bases.transpose(1,0,2).reshape(rows,-1)
    # An F_2-linear combination of the raw evaluations becomes an F_128 row.
    # Cache combinations in blocks of eight to avoid repeating matrix products.
    tables = []
    for start in range(0,rows,8):
        block = matrix[start:start+8]
        table = np.zeros((256,twists*columns), dtype=np.uint8)
        for mask in range(1,1 << len(block)):
            j = (mask & -mask).bit_length() - 1
            table[mask] = table[mask ^ (1<<j)] ^ block[j]
        tables.append(table)
    counts = [0,5]
    with path.open('wb') as stream:
        np.array([twists,columns,bound],dtype='<i4').tofile(stream)

        def write_place(degree, values):
            packed = np.packbits(values,axis=1,bitorder='little')
            result = np.zeros((len(values),twists*columns), dtype=np.uint8)
            for j,table in enumerate(tables):
                result ^= table[packed[:,j]]
            np.array([degree,len(values)//degree],dtype='<i4').tofile(stream)
            result.tofile(stream)

        for cusp in range(5):
            write_place(1, cusp_values(cusp,bound))
        assert not list(field(1).affine_points())
        for degree in range(2,bound+1):
            F = field(degree)
            count = 0
            for point in F.places_of_exact_degree():
                jets = affine_values(F,point,bound//degree)
                bits = (jets[:,None,:] >> np.arange(degree)[None,:,None]) & 1
                write_place(degree,bits.reshape(-1,rows))
                count += 1
            counts.append(count)
    print('Closed places in degrees 1 through',bound,counts[1:],flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, default=19, choices=[16,17,18,19])
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        data, rank_program = tmp/'evaluation.bin', tmp/'rank121'
        evaluation_data(data,args.degree)
        subprocess.run(['c++','-std=c++17','-O3',str(DIRECTORY/'rank121.cpp'),'-o',str(rank_program)],
                       check=True)
        subprocess.run([str(rank_program),str(data)],check=True)
