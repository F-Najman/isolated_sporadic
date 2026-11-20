// Characteristic polynomial of a modular unit on a residual cusp fibre.
// vanhoeij131.py supplies cusp expansions and an auxiliary prime p.
#include <flint/nmod_poly.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using Integer = long long;

// Product of (1-t^m)^e[m mod level] through degree precision-1.
// Compute t times the logarithmic derivative, integrate, and exponentiate.
void unit_series(nmod_poly_t result, const std::vector<int>& exponents,
                 int precision, int p, const std::vector<int>& inverses) {
    std::vector<Integer> logarithmic(precision);
    int level = exponents.size();
    for (int m = 1; m < precision; ++m) {
        int e = exponents[m % level];
        if (e == 0) continue;
        for (int k = m; k < precision; k += m)
            logarithmic[k] -= Integer(m)*e;
    }
    nmod_poly_t log_unit;
    nmod_poly_init(log_unit, p);
    for (int k = 1; k < precision; ++k) {
        Integer numerator = (logarithmic[k] % p + p) % p;
        nmod_poly_set_coeff_ui(log_unit, k, numerator*inverses[k] % p);
    }
    nmod_poly_exp_series(result, log_unit, precision);
    nmod_poly_clear(log_unit);
}

// If h = sign * t^(-order) * H, subtract the constant term of
// h^k * G from moment k. The baby/giant powers share series products.
void subtract_residues(std::vector<Integer>& moments, nmod_poly_t H,
                       nmod_poly_t G, int order, int sign, int p) {
    int degree = moments.size();
    int precision = degree*order + 1;
    int step = int(std::sqrt(degree)) + 1;
    std::vector<std::vector<int>> babies(step, std::vector<int>(precision));
    nmod_poly_t power, block;
    nmod_poly_init(power, p);
    nmod_poly_init(block, p);
    nmod_poly_one(power);
    for (int i = 0; i < step; ++i) {
        for (int a = 0; a < nmod_poly_length(power); ++a)
            babies[i][a] = nmod_poly_get_coeff_ui(power, a);
        nmod_poly_mullow(power, power, H, precision);
    }
    nmod_poly_set(block, G);
    for (int j = 0; j <= degree/step; ++j) {
        for (int i = 0; i < step; ++i) {
            int k = step*j+i;
            if (k == 0 || k > degree) continue;
            int target = order*k;
            int first = std::max(0, target-int(nmod_poly_length(block))+1);
            Integer coefficient = 0;
            for (int a = first; a <= target; ++a)
                coefficient += Integer(babies[i][a])*nmod_poly_get_coeff_ui(block, target-a);
            Integer signed_value = (sign < 0 && k % 2) ? -coefficient : coefficient;
            moments[k-1] -= signed_value % p;
        }
        nmod_poly_mullow(block, block, power, precision);
    }
    nmod_poly_clear(power);
    nmod_poly_clear(block);
}

int main() {
    int level, degree, p, cusps;
    std::cin >> level >> degree >> p >> cusps;
    assert(degree > 0 && degree < p && p < 1000000);
    std::vector<int> inverses(p);
    inverses[1] = 1;
    for (int k = 2; k < p; ++k)
        inverses[k] = p-Integer(p/k)*inverses[p%k] % p;
    std::vector<Integer> moments(degree);
    for (auto& m : moments) std::cin >> m;
    for (int cusp = 0; cusp < cusps; ++cusp) {
        int v, sign_f, w, pole_h, sign_h;
        std::cin >> v >> sign_f >> w >> pole_h >> sign_h;
        std::vector<int> f_exponents(level), h_exponents(level);
        for (auto& e : f_exponents) std::cin >> e;
        for (auto& e : h_exponents) std::cin >> e;
        int cancellation = w-std::min(v, 0);
        int precision = degree*pole_h+1;
        assert(cancellation >= 0 && pole_h > 0 && precision+cancellation+2 < p);
        nmod_poly_t f, h, derivative, G;
        nmod_poly_init(f, p);
        nmod_poly_init(h, p);
        nmod_poly_init(derivative, p);
        nmod_poly_init(G, p);

        // f becomes the unit part of sign*U-1, with valuation w.
        unit_series(f, f_exponents, precision+cancellation+2, p, inverses);
        if (sign_f < 0) nmod_poly_neg(f, f);
        if (v > 0) nmod_poly_shift_left(f, f, v);
        int one_exponent = std::max(-v, 0);
        auto old = nmod_poly_get_coeff_ui(f, one_exponent);
        nmod_poly_set_coeff_ui(f, one_exponent, (old+p-1) % p);
        nmod_poly_shift_right(f, f, cancellation);
        assert(nmod_poly_get_coeff_ui(f, 0) != 0);
        nmod_poly_derivative(derivative, f);
        nmod_poly_shift_left(derivative, derivative, 1);
        nmod_poly_div_series(G, derivative, f, precision);
        nmod_poly_set_coeff_ui(G, 0, (w % p+p) % p);
        unit_series(h, h_exponents, precision, p, inverses);
        subtract_residues(moments, h, G, pole_h, sign_h, p);
        nmod_poly_clear(f);
        nmod_poly_clear(h);
        nmod_poly_clear(derivative);
        nmod_poly_clear(G);
    }
    assert(std::cin);
    for (auto& m : moments) m = (m % p+p) % p;
    std::vector<Integer> coeffs(degree+1);
    coeffs[0] = 1;
    for (int k = 1; k <= degree; ++k) {
        Integer total = 0;
        for (int i = 0; i < k; ++i)
            total += coeffs[i]*moments[k-1-i];
        coeffs[k] = ((p-total % p) % p)*inverses[k] % p;
    }
    for (auto c : coeffs) std::cout << c << ' ';
    std::cout << '\n';
}
