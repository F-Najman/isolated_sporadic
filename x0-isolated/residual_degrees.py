"""Factor the residual zero divisors of the canonical differentials.

Run with sage -python residual_degrees.py.
"""
from sage.all import CuspForms, Gamma0, Matrix, PolynomialRing, QQ


def residual_degrees(p):
    forms = CuspForms(Gamma0(p), 2).echelon_basis()
    g = len(forms)
    length = 4*(2*g-2) + 1
    ring = PolynomialRing(QQ, "q")
    expansions = [f.q_expansion(length+1) for f in forms]
    sections = [ring([f[i] for i in range(1, length+1)]) for f in expansions]

    def product(space, multipliers):
        rows = [(ring(list(row))*s).truncate(length).padded_list(length)
                for row in space.rows() for s in multipliers]
        return Matrix(QQ, rows)

    def reduce_modulo(space, subspace):
        return space - space.matrix_from_columns(subspace.pivots())*subspace

    v1 = Matrix(QQ, [f.padded_list(length) for f in sections])
    v2 = product(v1, sections).row_space().basis_matrix()
    v3 = product(v2, sections).row_space().basis_matrix()
    v4 = product(v3, sections).row_space().basis_matrix()
    assert [v.nrows() for v in (v2,v3,v4)] == [(2*k-1)*(g-1) for k in (2,3,4)]

    # Restriction to the zero scheme of the last echelon form.
    multiples3 = product(v2, [sections[-1]]).row_space().basis_matrix()
    multiples4 = product(v3, [sections[-1]]).row_space().basis_matrix()
    source = reduce_modulo(v3, multiples3).row_space().basis_matrix()
    assert source.nrows() == v4.nrows()-multiples4.nrows() == 2*g-2

    def multiplication(s):
        return reduce_modulo(product(source, [s]), multiples4)

    # The first differential trivializes the canonical bundle on the zero
    # scheme. Multiplication by the ratio of the second to the first is an
    # endomorphism, with an eigenvalue 0 of multiplicity at least g-1 at infinity.
    denominator = multiplication(sections[0])
    assert denominator.rank() == 2*g-2
    pivots = denominator.pivots()
    numerator = multiplication(sections[1]).matrix_from_columns(pivots)
    operator = numerator * denominator.matrix_from_columns(pivots).inverse()
    characteristic = operator.charpoly("T")
    T = characteristic.parent().gen()
    residual, remainder = characteristic.quo_rem(T**(g-1))
    assert remainder == 0 and residual.degree() == g-1
    assert residual.is_irreducible()
    print(p, g, residual.degree(), residual, flush=True)


print("N genus degree_of_closed_point residual_polynomial", flush=True)
for level in [53, 61, 79, 89, 101]:
    residual_degrees(level)
