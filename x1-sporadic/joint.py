"""Exact simultaneous-prime incidence inequalities in a Hecke coefficient field."""
from collections import deque
from fractions import Fraction
from functools import cache
from pathlib import Path
import json
import sympy as sp

X, Y = sp.symbols('X Y')
DATA = Path(__file__).resolve().parent / 'data'


@cache
def _raw():
    return {r['label']: r for r in json.loads((DATA / 'hecke_spectral.json').read_text())}


@cache
def _trace_polynomials(label):
    """Return the field polynomial and the elements a_p² / epsilon(p) - 2p."""
    data = _raw()[label]
    assert not data['cyclotomic']
    f = sp.Poly.from_list(data['field_poly'][::-1], X, domain='QQ')
    basis = [sp.Poly(X**i, X, domain='QQ') for i in range(f.degree())]
    if not data['power_basis']:
        basis = [sp.Poly.from_list(row[::-1], X, domain='QQ').mul_ground(sp.Rational(1, d))
                 for row, d in zip(data['nums'], data['dens'])]
    zero, one = sp.Poly(0, X, domain='QQ'), sp.Poly(1, X, domain='QQ')

    def element(coords):
        return sum((b*c for b, c in zip(basis, coords)), zero) % f

    level = int(label.split('.')[0])
    character = {1:one}
    if data['char_values'] is not None:
        generators = [(a, element(coords)) for a, coords in data['char_values']]
        queue = deque([1])
        while queue:
            a = queue.popleft()
            for b, v in generators:
                ab = a*b % level
                value = (character[a]*v) % f
                if ab not in character:
                    character[ab] = value
                    queue.append(ab)
                else:
                    assert character[ab] == value
        assert len(character) == sp.totient(level)
    traces = {}
    for key, coords in data.items():
        if not (key.startswith('a') and key[1:].isdigit()):
            continue
        p = int(key[1:])
        if level % p == 0:
            continue
        chi = character[p % level] if data['char_values'] is not None else one
        a = element(coords)
        traces[p] = a*a*chi.invert(f) % f - 2*p
    return f, traces


@cache
def _positive_conjugate(label, i, j, wi, wj):
    f, ts = _trace_polynomials(label)
    v = ts[i]*wi + ts[j]*wj
    if v.degree() <= 0:
        return v.nth(0) > 0
    P = sp.Poly(sp.resultant(f.as_expr(), Y-v.as_expr(), X), Y).monic()
    # Every root is real. All roots are nonpositive iff all coefficients
    # of this monic polynomial are nonnegative.
    return any(a < 0 for a in P.all_coeffs())


def joint_excluded(form, d, counts):
    """Exclude a positive-rank factor by a pair of Frobenius inequalities."""
    label = form['label']
    gaps = {p:S - (p*p+1)*d for p,S in counts.items()}
    for i, positive in gaps.items():
        if positive <= 0:
            continue
        for j, nonpositive in gaps.items():
            if nonpositive > 0:
                continue
            # Scale the weights to (wi, wj) with wj=1.
            wi = Fraction(-nonpositive, positive)
            if form['dim'] == 1:
                tr = form['initial_traces']
                if wi*(tr[i-1]**2 - 2*i) + (tr[j-1]**2 - 2*j) > 0:
                    return True
            elif _positive_conjugate(label, i, j, wi, 1):
                return True
    return False
