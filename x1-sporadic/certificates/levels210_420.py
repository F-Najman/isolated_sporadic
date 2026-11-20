"""Numerical checks for the proof of sporadic CM points at 210 and 420."""

from math import gcd


def small_prime_gonality(prime, nodes):
    """Lower bound from the two-component model in the proof (p = 2 or 3)."""
    # A function of degree b on the quotient must satisfy (p^2+1)b >= nodes.
    base_degree = (nodes + prime**2) // (prime**2 + 1)
    if prime == 2:
        distinct_comparisons = (2*nodes + 2) // 3
        identical_comparison = 3*base_degree
        return min(distinct_comparisons, identical_comparison)
    if prime == 3:
        noninvariant_degree = (nodes + 1) // 2
        other_degree = min(noninvariant_degree, 2*base_degree)
        return min(nodes, 8*base_degree, noninvariant_degree + other_degree)
    raise ValueError("The two formulas apply to primes 2 and 3.")


def modulo_sign(x, y, level):
    return min((x % level, y % level), (-x % level, -y % level))


def supersingular_nodes(tame_level, automorphisms):
    primitive = sum(gcd(gcd(a, b), tame_level) == 1
                    for a in range(tame_level) for b in range(tame_level))
    nodes, remainder = divmod(primitive, automorphisms)
    assert remainder == 0
    return nodes


def cm_orbit(level, root=48):
    # Use the basis 1,s with s^2=-6. The marked point is (s-root)/level.
    orbit = set()
    for a in range(level):
        for b in range(level):
            if gcd(a*a + 6*b*b, level) == 1:
                orbit.add(modulo_sign(-root*a - 6*b, a - root*b, level))
    conjugate = {modulo_sign(x, -y, level) for x, y in orbit}
    assert orbit.isdisjoint(conjugate)
    # The two reduced forms of discriminant -24 give two conjugate j-values.
    return len(orbit), 2*len(orbit | conjugate)


if __name__ == "__main__":
    for n, prime, aut_order, expected in (
        (210, 2, 24, (384, 231, 116, 24, 96)),
        (420, 3, 12, (1152, 808, 404, 96, 384)),
    ):
        nodes = supersingular_nodes(n // prime, aut_order)
        gonality = small_prime_gonality(prime, nodes)
        density = (gonality + 1) // 2
        hilbert_orbit, point_degree = cm_orbit(n)
        assert (nodes, gonality, density, hilbert_orbit, point_degree) == expected
        assert point_degree < density
        print(f"X_1({n}): {nodes} nodes at {prime}, geometric gonality >= {gonality}, "
              f"density >= {density}, CM point degree {point_degree}.")

    nodes = supersingular_nodes(70, 12)
    assert nodes == 288 and small_prime_gonality(3, nodes) == 202
    print("Independent reduction at 3 for X_1(210): geometric gonality >= 202, density >= 101.")
