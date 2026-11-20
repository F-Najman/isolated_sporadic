"""Run with sage -python degree52_cusps83.py [EIS83_NPZ].

Use the supplied cusp valuations and jets, or read the valuations from the
earlier level-83 archive. Regenerate the jets with eisenjets.py 83 124.
"""
import argparse
import ast
from collections import Counter
from itertools import combinations_with_replacement
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from sage.all import (CuspForms, DirichletGroup, EllipticCurve, GF,
                      ModularSymbols, QQ, ZZ, gcd, identity_matrix, matrix,
                      vector)


def check_quotient_orders(cusp_index):
    chi = (DirichletGroup(83).gen()**2).minimize_base_ring()
    space = ModularSymbols(chi, 2, sign=1).cuspidal_submodule()
    trivial = ModularSymbols(83, 2, sign=1).cuspidal_submodule()
    assert space.dimension() == 6 and trivial.dimension() == 7
    curve = EllipticCurve("83a1")
    running = ZZ(0)
    for p in (3, 5, 7, 29):
        prym = ZZ(((1 + p*chi(p))*identity_matrix(chi.base_ring(), 6)
                   - space.hecke_matrix(p)).det().norm())
        base = ZZ(((1+p)*identity_matrix(QQ, 7) - trivial.hecke_matrix(p)).det())
        elliptic = p + 1 - curve.ap(p)
        order, remainder = (base*prym).quo_rem(elliptic)
        assert remainder == 0 and order % cusp_index == 0
        if p == 3:
            assert gcd(order, 7*(9**41 + 1)) == 1
        running = gcd(running, order)
        print(f"Through p={p}, quotient-order gcd / cusp index = "
              f"{running//cusp_index}", flush=True)
    assert running == cusp_index


def diamond(divisor, multiplier):
    image = [divisor[0]] + [0]*41
    for j in range(1, 42):
        t = multiplier*j % 83
        image[min(t, 83-t)] = divisor[j]
    return vector(ZZ, image)


def cusp(j):
    v = vector(ZZ, 42)
    v[j] = 1
    return v


class CuspLattice:
    def __init__(self, valuations):
        units = matrix(ZZ, valuations)
        assert units.dimensions() == (41, 42) and units[0, 0] == -1
        assert all(41*r[0] + sum(r[1:]) == 0 for r in units)
        self.units = units
        minor = units.delete_columns([1])
        self.index = abs(minor.det())
        self.minor_inverse = minor.inverse()
        assert gcd(self.index, 7*(9**41 + 1)) == 1
        rational = units.matrix_from_columns([0]).left_kernel_matrix()*units
        self.basis = rational.delete_columns([0]).LLL()
        self.rational_inverse = self.basis.delete_columns([0]).inverse()
        self.integer_basis = np.array(self.basis, dtype=np.int64)
        orthogonal, self.triangular = np.linalg.qr(self.integer_basis.T)
        self.orthogonal = orthogonal

    def principal(self, v):
        assert 41*v[0] + sum(v[1:]) == 0
        return all(x in ZZ for x in vector([v[0], *v[2:]])*self.minor_inverse)

    def reduce_large(self, v):
        v = v + v[0]*self.units.row(0)
        assert v[0] == 0
        v = vector(QQ, v[1:])
        centre = vector(QQ, [sum(v)/41]*41)
        coords = vector((v-centre)[1:])*self.rational_inverse
        small = vector(QQ, [c - c.round() for c in coords])*self.basis + centre
        return np.array(vector(ZZ, small), dtype=np.int64)

    def shorten(self, v):
        # The approximation only selects integer combinations of certified
        # modular units. The divisor and the row operations stay integral.
        coords = self.orthogonal.T @ (v - v.sum()/41)
        for j in range(39, -1, -1):
            coeff = int(np.rint(coords[j]/self.triangular[j, j]))
            if coeff:
                v = v - coeff*self.integer_basis[j]
                coords[:j+1] -= coeff*self.triangular[:j+1, j]
        return v

    def preimage(self, target, degree, mult):
        # Invert 7(9+tau) on the finite degree-zero cusp-class group.
        base = degree*cusp(1)
        error = target - 7*(9*base + diamond(base, mult))
        assert 41*error[0] + sum(error[1:]) == 0
        total = vector(ZZ, 42)
        term = error
        for j in range(41):
            total += (-1)**j * ZZ(9)**(40-j) * term
            term = diamond(term, mult)
        factor = ZZ(7*(9**41+1)).inverse_mod(self.index)
        v = self.reduce_large(base + factor*total)
        full = vector(ZZ, [0, *v])
        assert self.principal(7*(9*full + diamond(full, mult)) - target)
        return v


class JetTest:
    def __init__(self, path):
        with np.load(path) as data:
            assert int(data["prime"]) == 167
            self.forms = matrix(GF(167), data["rows"])
            self.profile = [tuple(map(int, r)) for r in data["profile"]]
        self.cusps = self.profile[1:]
        assert len(self.cusps) == 41 and all(w == 1 for w, _, _ in self.cusps)

    def excludes(self, v):
        if min(v) >= 0:
            return False
        rows = []
        for m, (_, precision, start) in zip(v, self.cusps):
            if m > precision:
                return False
            rows.extend(range(start, start + max(0, int(m))))
        return self.forms.matrix_from_rows(rows).rank() == len(rows)

    def check_all_incidence_degrees(self):
        # F*O contains every geometric cusp with multiplicity one.
        base_rows = [start+i*prec for w, prec, start in self.profile for i in range(w)]
        base = self.forms.matrix_from_rows(base_rows)
        base_rank = base.rank()
        remaining = self.forms*base.right_kernel_matrix().transpose()

        def upper(extra):
            rows = []
            for c, mult in extra.items():
                w, prec, start = self.profile[c]
                assert 1+mult <= prec
                for i in range(w):
                    first = start+i*prec+1
                    rows.extend(range(first, first+mult))
            return 1+len(base_rows)+len(rows)-base_rank-remaining.matrix_from_rows(rows).rank()

        assert upper({j: 4 for j in range(1, 42)}) == 2
        # In dimension at most two the rational diamond representation is
        # trivial. Compute the exact dimension on X_0(83) using W83=-U83.
        expansions = CuspForms(83, 2).q_expansion_basis(prec=84)
        base_jets = matrix(QQ, [[f[i] for i in (1, 2, 3, 4, 5, 83)]
                                for f in expansions])
        assert base_jets.rank() == 6
        hist = Counter([upper({1: 52}), upper({0: 1, 1: 11})])
        for j in range(2, 42):
            inv = int(ZZ(j).inverse_mod(83))
            stop = 27 if j <= min(inv, 83-inv) else 26
            for m in range(1, stop):
                hist[upper({1: 52-m, j: m})] += 1
        assert hist == {2: 1022}
        print("All k: upper bound 2 for each of 1022 two-cusp orbits; "
              "dimension 1 for F*O plus four times the rational cusp sum.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    folder = Path(__file__).resolve().parent
    parser.add_argument("jets", type=Path, nargs="?", default=folder/"eis83.npz")
    parser.add_argument("--archive", type=Path, default=Path(__file__).with_name(
        "level-83-certificate.zip"))
    args = parser.parse_args()
    valuations = folder/"certificates"/"valuations83.txt"
    if valuations.exists():
        table = valuations.read_text()
    else:
        with ZipFile(args.archive) as archive:
            table = archive.read("certificates/valuations83.txt").decode()
    lattice = CuspLattice(ast.literal_eval(table))
    check_quotient_orders(lattice.index)
    jets = JetTest(args.jets)
    jets.check_all_incidence_degrees()
    fibre = vector(ZZ, [1]*42)
    for mult in (3, int(ZZ(3).inverse_mod(83))):
        first = lattice.preimage(32*cusp(1) + 44*fibre, 52, mult)
        steps = [np.zeros(41, dtype=np.int64)]
        steps.extend(lattice.preimage(cusp(j) - cusp(1), 0, mult)
                     for j in range(2, 42))
        tested = 0
        unresolved = []
        for j in range(1, 42):
            step = steps[j-1]
            # Up to diamonds one cusp is 1 and has weight at least 16.
            # In the equal-weight case, identify a ratio with its inverse.
            inv = int(ZZ(j).inverse_mod(83))
            stop = 17 if j <= min(inv, 83-inv) else 16
            for m in ([0] if j == 1 else range(1, stop)):
                v = lattice.shorten(first + m*step)
                assert v.sum() == 52
                tested += 1
                if not jets.excludes(v):
                    unresolved.append((j, m))
        assert tested == 621
        print(f"Diamond multiplier {mult}: {tested} classes, "
              f"{len(unresolved)} inconclusive", unresolved, flush=True)
        assert not unresolved
        tested = 0
        unresolved = []
        for i, j, h in combinations_with_replacement(range(41), 3):
            v = lattice.shorten(first + steps[i] + steps[j] + steps[h])
            assert v.sum() == 52
            tested += 1
            if not jets.excludes(v):
                unresolved.append((i+1, j+1, h+1))
        assert tested == 12341
        print(f"Diamond multiplier {mult}: {tested} degree-three residual tails, "
              f"{len(unresolved)} inconclusive", unresolved, flush=True)
        assert not unresolved


if __name__ == "__main__":
    main()
