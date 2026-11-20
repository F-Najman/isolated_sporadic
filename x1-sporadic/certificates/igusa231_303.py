"""Numerical checks for the elliptic incidence and CM proof at 231 and 303."""

import argparse
from fractions import Fraction
from math import ceil


def density_bound(level, nodes, supersingular_at_two):
    """First integer surviving the necessary inequalities in the proof."""
    if level not in (231, 303) or supersingular_at_two != 4 * nodes:
        raise ValueError("Expected the level and counts of one of the two Igusa models")
    # At 77 the trace term is strictly negative since the sum map moves.
    small_branch = nodes // 10 + 1 if level == 231 else ceil(Fraction(nodes, 12))
    noninvariant = ceil(Fraction(nodes, 2))
    first_with_noninvariant = noninvariant + min(noninvariant, 2 * small_branch)

    for d in range(1, first_with_noninvariant + 1):
        # Larger simple factors and the two new elliptic classes at 303.
        if supersingular_at_two <= 7 * d:
            return d
        if level == 303 and (4 * supersingular_at_two <= 27 * d
                             or supersingular_at_two <= 5 * d):
            return d
        if d == first_with_noninvariant:
            return d
        if d >= 8 * small_branch:  # one Frobenius comparison is an identity
            return d
        # For a fixed sum b1+b2 the largest possible minimum of the two
        # coincidence bounds occurs when the component degrees differ by at most one.
        for total in range(2 * small_branch, d // 2 + 1):
            matched_max = min(nodes, 2 * total + 4 * (total // 2) - small_branch)
            missing_min = nodes - matched_max
            if missing_min == 1:
                missing_min = 2  # neither tame curve has a rational supersingular mark
            if 2 * total + missing_min <= d:
                return d


def verify_curves():
    """Enumerate conductor divisors using Cremona and recompute ranks and traces."""
    from sage.all import EllipticCurve, QuadraticField, divisors
    from sage.databases.cremona import CremonaDatabase

    database = CremonaDatabase()
    assert QuadraticField(-6, "s").class_number(proof=True) == 2
    conductors = sorted(set(divisors(231)) | set(divisors(303)))
    positive = {}
    for conductor in conductors:
        for label, (coefficients, listed_rank, _) in database.allcurves(conductor).items():
            if label != label.rstrip("0123456789") + "1":
                continue
            curve = EllipticCurve(coefficients)
            assert curve.conductor() == conductor
            rank = curve.rank(use_database=False, proof=True)
            assert rank == listed_rank
            if rank:
                name = str(conductor) + label
                positive[name] = (int(rank), *(int(curve.ap(p)) for p in (2, 3, 5)))
    assert positive == {
        "77a1": (1, 0, -3, -1),
        "101a1": (1, 0, -2, -1),
        "303a1": (1, 0, 1, -3),
        "303b1": (1, -2, 1, -1),
    }
    for name, (rank, a2, a3, a5) in positive.items():
        print(f"{name}: rank {rank}, traces at 2, 3, 5 are {a2}, {a3}, {a5}")


if __name__ == "__main__":
    from levels210_420 import cm_orbit, supersingular_nodes

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-curves", action="store_true",
                        help="also enumerate isogeny classes and prove the ranks with Sage")
    args = parser.parse_args()
    if args.verify_curves:
        verify_curves()

    for level, expected in [(231, (480, 1920, 266, 60, 240)),
                            (303, (850, 3400, 462, 100, 400))]:
        nodes = supersingular_nodes(level // 3, 12)
        at_two = supersingular_nodes(level, 24)
        density = density_bound(level, nodes, at_two)
        assert (216**2 + 6) % level == 0
        hilbert_degree, degree = cm_orbit(level, root=216)
        assert (nodes, at_two, density, hilbert_degree, degree) == expected
        assert degree < density
        print(f"X_1({level}): {nodes} nodes at 3, {at_two} supersingular points at 2, "
              f"density at least {density}, CM orbit degree {degree}")
