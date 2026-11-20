"""Compute Fourier coefficients of Gamma1(N) orbit sums of the normalized Weierstrass function from its Tate expansion. Indices are pairs in (Z/NZ)^2 and exponents use the geometric cusp width."""

import math
import numpy as np
from functools import lru_cache

def extra_eis(N, p, roots):
    uv = []
    ind = []
    nor = 0
    for u in range(N):
        h = math.gcd(N, u)
        for b in range(h):
            if u == 0 and b == 0:
                continue
            for v in range(b, N, h):
                uv.append((u, v))
                ind.append(nor)
            nor += 1
    uv = np.array(uv)
    ind = np.array(ind)
    root = np.array(roots)
    beta = np.arange(N)

    @lru_cache(None)
    def tor(den, k):
        m = den * k
        A = np.zeros((N, N), np.int64)
        if k == 0:
            A += pow(12, -1, p)
            for b in range(1, N):
                A[0, b] += roots[b] * pow((1 - roots[b]) ** 2, -1, p) % p
        else:
            for x in range(1, m + 1):
                if m % x:
                    continue
                nr = m // x
                A[x % N, :] += nr * root[beta * nr % N]
                A[-x % N, :] += nr * root[-beta * nr % N]
                if x % N == 0:
                    A -= 2 * nr
        return A % p

    def get(c, b, h, k):
        a, t = next(((a, t) for a in range(N) for t in range(N) if (a * b - t * c) % N == 1))
        A = tor(h, k)
        F = A[(uv[:, 0] * a + uv[:, 1] * c) % N, (uv[:, 0] * t + uv[:, 1] * b) % N]
        return (np.bincount(ind, weights=F, minlength=nor).astype(np.int64) % p).tolist()
    return get
