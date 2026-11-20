"""Exact orders of the modular unit and its special fibre at the rational cusps.

The functions E_g and their transformation laws are in Propositions 1–3
of Yifan Yang, arXiv:0712.0629.  We use the model with a mu_N marking.
Its infinity cusps are rational, and dual isogeny identifies it with X_1(N).
"""
from collections import Counter


def cusp_fibres(p, numerator, denominator):
    e = Counter(numerator)
    e.subtract(denominator)
    e = {i: n for i, n in e.items() if n}
    assert sum(e.values()) == 0 and sum(i*i*n for i,n in e.items()) % p == 0
    degree = 0
    fibres = {1: 0, -1: 0}
    for cusp in range(1, (p+1)//2):
        a = cusp if cusp % 2 else cusp+p
        residues = {i: a*i % p for i in e}
        v, remainder = divmod(sum(n*residues[i]*(residues[i]-p)
                                  for i,n in e.items()), 2*p)
        assert not remainder
        degree += max(-v, 0)
        if v:
            continue
        # For odd a, the exponential term of Yang's transformation is 1.
        sign = (-1)**(sum(n*((a*i)//p) for i,n in e.items()) % 2)
        contact = min(min(r, p-r) for r in residues.values())
        assert sum(min(r, p-r)==contact for r in residues.values()) == 1
        fibres[sign] += contact
    return degree, fibres


d, fibres = cusp_fibres(107, (10,52), (3,21))
assert d == 153 and fibres == {1: 0, -1: 60}
assert d-fibres[-1] == 93 and 93 % 53 != 0
print("N=107, noncuspidal degree at most", d-fibres[-1])

d, fibres = cusp_fibres(83, (35,40), (3,34))
assert d == 95 and fibres == {1: 38, -1: 0}
assert d-fibres[1] == 57 and 57 % 41 != 0
print("N=83, noncuspidal degree at most", d-fibres[1])
