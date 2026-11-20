"""Finite fields and local coordinates for Fisher's level-11 curve in characteristic two."""
from functools import cache
import numpy as np


class BinaryField:
    """Use the power basis of the first primitive binary polynomial of the given degree."""

    def __init__(self, degree):
        self.degree = degree
        self.order = 1 << degree
        self.period = self.order - 1
        for modulus in range(self.order + 1, 2 * self.order, 2):
            x = 1
            powers = [1]
            for _ in range(self.period):
                x <<= 1
                if x & self.order:
                    x ^= modulus
                if x == 1:
                    break
                powers.append(x)
            if len(powers) == self.period:
                break
        assert len(set(powers)) == self.period and x == 1
        self.modulus = modulus
        self.exp = np.array(powers * 2 + [1], dtype=np.int32)
        self.log = np.zeros(self.order, dtype=np.int32)
        self.log[powers] = np.arange(self.period)
        # Higher-order jets are needed only at points of small degree.
        if degree <= 9:
            self.table = self.exp[self.log[:, None] + self.log[None, :]]
            self.table[0, :] = 0
            self.table[:, 0] = 0

    def power(self, value, exponent):
        if exponent == 0:
            return 1
        if value == 0:
            assert exponent > 0
            return 0
        return int(self.exp[(int(self.log[value]) * exponent) % self.period])

    def multiply(self, a, b):
        return int(self.exp[self.log[a] + self.log[b]]) if a and b else 0

    def affine_points(self):
        """All points away from the five ramified cusps; x, y, t are field elements."""
        nonzero = np.arange(1, self.order, dtype=np.int32)
        squares = self.exp[(2 * self.log[nonzero]) % self.period]
        cubes = self.exp[(3 * self.log[nonzero]) % self.period]
        rhs = np.zeros(self.order, dtype=np.int32)
        rhs[1:] = squares ^ cubes
        artin_schreier = np.zeros(self.order, dtype=np.int32)
        artin_schreier[1:] = squares ^ nonzero
        y_root = np.full(self.order, -1, dtype=np.int32)
        y_root[artin_schreier] = np.arange(self.order)
        t_root = np.full(self.order, -1, dtype=np.int32)
        t_root[self.exp[(11 * self.log[nonzero]) % self.period]] = nonzero
        for x in range(2, self.order):
            first_y = int(y_root[rhs[x]])
            if first_y == -1:
                continue
            for y in (first_y, first_y ^ 1):
                w = x ^ y ^ 1
                assert y != 0 and w != 0
                log_rhs = self.log[x] + 2 * self.log[y] + self.log[x ^ 1] + 3 * self.log[w]
                t = int(t_root[self.exp[log_rhs % self.period]])
                if t == -1:
                    continue
                yield x, y, t
                if self.period % 11 == 0:
                    for i in range(1, 11):
                        yield x, y, int(self.exp[(self.log[t] + i * (self.period // 11)) % self.period])

    def frobenius(self, point):
        return tuple(self.power(a, 2) for a in point)

    def rotate(self, point):
        """Cyclically permute Fisher's five Klein coordinates, then recover (x,y,t)."""
        x, y, t = point
        w = x ^ y ^ 1

        def product(*factors):
            return int(self.exp[sum(int(self.log[a]) for a in factors) % self.period])

        def ratio(a, b):
            return self.multiply(a, self.power(b, -1))

        klein = (product(x, y, x ^ 1), ratio(product(y, self.power(t, 4)), w),
                 ratio(self.power(t, 7), self.power(w, 2)), product(x, y, t), product(y, t, t))
        a, b, c, d, e = klein[-1:] + klein[:-1]
        return (ratio(product(a, b, d), product(c, c, e)),
                ratio(product(a, b, b, b), product(c, c, c, e)),
                ratio(product(a, b), product(c, c)))

    def places_of_exact_degree(self):
        remaining = set(self.affine_points())
        for point in sorted(remaining):
            if point not in remaining:
                continue
            length = 0
            current = point
            while current in remaining:
                remaining.remove(current)
                length += 1
                current = self.frobenius(current)
            assert current == point and self.degree % length == 0
            if length == self.degree:
                yield point

    def series_product(self, a, b):
        result = np.zeros(len(a), dtype=np.int32)
        for i, coefficient in enumerate(a):
            if coefficient:
                result[i:] ^= self.table[coefficient, b[:len(a)-i]]
        return result

    def series_inverse(self, a):
        result = np.zeros(len(a), dtype=np.int32)
        result[0] = self.power(a[0], -1)
        for i in range(1, len(a)):
            cross = np.bitwise_xor.reduce(self.table[a[1:i+1], result[:i][::-1]])
            result[i] = self.multiply(result[0], cross)
        return result

    def series_power(self, a, exponent):
        if exponent < 0:
            return self.series_power(self.series_inverse(a), -exponent)
        result = np.zeros(len(a), dtype=np.int32)
        result[0] = 1
        while exponent:
            if exponent % 2:
                result = self.series_product(result, a)
            exponent //= 2
            if exponent:
                a = self.series_product(a, a)
        return result

    def eleventh_root_series(self, series, constant):
        root = np.zeros(len(series), dtype=np.int32)
        root[0] = constant
        inverse_derivative = self.power(constant, -10)
        for i in range(1, len(series)):
            error = series[i] ^ self.series_power(root, 11)[i]
            root[i] = self.multiply(error, inverse_derivative)
        assert np.array_equal(self.series_power(root, 11), series)
        return root

    def affine_jets(self, point, length):
        """Unit series for x,y,x-1,x-y-1,t in the local parameter x-x(P)."""
        x, y, t = point
        X = np.zeros(length, dtype=np.int32)
        X[0] = x
        if length > 1:
            X[1] = 1
        Xm = X.copy()
        Xm[0] ^= 1
        rhs = self.series_product(self.series_product(X, X), Xm)
        Y = np.zeros(length, dtype=np.int32)
        Y[0] = y
        for i in range(1, length):
            Y[i] = rhs[i] ^ (self.power(Y[i//2], 2) if i % 2 == 0 else 0)
        W = Xm ^ Y
        radicand = self.series_product(self.series_product(X, self.series_power(Y, 2)),
                                       self.series_product(Xm, self.series_power(W, 3)))
        return X, Y, Xm, W, self.eleventh_root_series(radicand, t)


@cache
def field(degree):
    return BinaryField(degree)
