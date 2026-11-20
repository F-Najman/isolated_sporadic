"""Python 3. Verify a differential at the reduction of the degree-57 place."""
from collections import Counter


N = 83
F = (35, 40, 3, 34)
# Each row gives U = g_a g_b/(g_c g_d), followed by the coefficients
# of 83 dlog(g_i/g_1), 2 <= i <= 41, in F_3.
ROWS = [
    (1, 1, 1, 1, '0002212210000122101002020110102002112212'),
    (3, 34, 35, 40, '1200020200022101220000000212210000120111'),
    (5, 12, 35, 40, '1121002200220110210102201101100021100111'),
    (3, 26, 6, 20, '0102100000000000012100100002210221221202'),
    (18, 22, 9, 35, '1002001200022110011111121121200022112221'),
    (3, 38, 1, 37, '0201022211122010010121001002111111002020'),
    (3, 35, 7, 40, '0000210100100022000120202201000200200020'),
    (3, 33, 9, 41, '0111010010220212122002201002010100112100'),
    (3, 26, 35, 37, '1100210022122021022010122000012120200000'),
    (3, 28, 40, 40, '1111112101022110200011021102220011200000'),
]


def unit_at_cusp(indices, c):
    powers = Counter()
    order = parity = 0
    for i, multiplicity in zip(indices, (1, 1, -1, -1)):
        k, s = divmod(c*i, N)
        order += multiplicity*s*(s-N)
        parity += multiplicity*k
        powers[s] += multiplicity
        powers[N-s] += multiplicity
    assert order % (2*N) == 0
    return order//(2*N), (-1)**(parity % 2), powers


def unit_coefficients(data, length):
    """Coefficients through T^(length-1) after removing the initial power."""
    if length <= 0:
        return []
    _, sign, powers = data
    series = [sign] + [0]*(length-1)
    for residue, exponent in powers.items():
        for j in range(residue, length, N):
            for _ in range(abs(exponent)):
                if exponent > 0:
                    for k in range(length-1, j-1, -1):
                        series[k] = (series[k]-series[k-j]) % 3
                else:
                    for k in range(j, length):
                        series[k] = (series[k]+series[k-j]) % 3
    return series


def log_coefficients(i, c, length):
    """83 T g_i'/g_i, with the constant common to all i omitted."""
    s = c*i % N
    series = [s*(s-N)//2] + [0]*(length-1)
    for residue in (s, N-s):
        for j in range(residue, length, N):
            for k in range(j, length, j):
                series[k] -= N*j
    return series


trace = [0]*6
for c0 in range(1, 42):
    c = c0 if c0 % 2 else c0 + N
    f = unit_at_cusp(F, c)
    # The characteristic-zero order of f-1 is needed even if it jumps mod 3.
    v = {10:8, 18:27, 26:3}.get(c0, min(f[0], 0))
    last = 6-v
    alpha = Counter()  # coefficients relative to dT/T
    for *indices, weights in ROWS:
        assert len(weights) == 40
        u = unit_at_cusp(indices, c)
        length = last-u[0]+1
        if length <= 0:
            continue
        uc = unit_coefficients(u, length)
        first = log_coefficients(1, c, length)
        lc = [0]*length
        for i, a in enumerate(map(int, weights), 2):
            if a:
                li = log_coefficients(i, c, length)
                for k in range(length):
                    lc[k] += a*(li[k]-first[k])
        for k in range(length):
            alpha[u[0]+k] += sum(uc[j]*lc[k-j] for j in range(k+1))

    assert all(value % 3 == 0 for k, value in alpha.items() if k <= -v)
    low = min(alpha, default=last+1)
    fc = unit_coefficients(f, 7-low-f[0])
    for k in range(1, 7):
        trace[k-1] += sum(a*fc[k-t-f[0]] for t, a in alpha.items()
                          if 0 <= k-t-f[0] < len(fc)) - alpha[k]

trace = [a % 3 for a in trace]
assert sum(a*b for a, b in zip(trace, (2,-1,-2,-2,-1,-1))) % 3 == 1
print('First six coefficients of the differential traced to X0(83) at 0:', trace)
print('Its elliptic projection is nonzero modulo 3.')
