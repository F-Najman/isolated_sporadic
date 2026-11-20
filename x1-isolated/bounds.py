"""Proven upper multiples for the four rank-zero Jacobian orders."""
from sage.all import J1


for level, expected in [(13, 19), (16, 20), (18, 21), (20, 60)]:
    torsion = J1(level).rational_torsion_subgroup()
    bound = torsion.multiple_of_order_using_frobp(maxp=73)
    assert bound == expected
    print(level, bound)
