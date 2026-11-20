// Exact, conservative enumeration of lattice vectors in rational balls.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>

using Integer = std::int64_t;
using Wide = __int128;
constexpr int capacity = 128;
using Vector = std::array<Integer, capacity>;
int rank, dimension, degree;
Integer scale, budget_scale, norm_scale;
Vector basis[capacity], gram_coeff[capacity], offset, coeff, norm_lower;
std::uint64_t nodes = 0, leaves = 0;

void check_vector() {
    ++leaves;
    Integer norm = 0;
    for (int j = 0; j < dimension; ++j) {
        Integer x = 0;
        for (int h = 0; h < rank; ++h) x += coeff[h] * basis[h][j];
        norm += std::abs(x);
        if (norm > 2 * degree) return;
    }
    if (norm) {
        std::cerr << "Found a nonzero vector of L1 norm " << norm << '\n';
        std::exit(2);
    }
}

void search(int i, Integer budget, Integer tail_norm) {
    if (i < 0) {
        check_vector();
        return;
    }
    Integer center = offset[i];
    Integer nearest = (std::abs(center) + scale/2) / scale * (center > 0 ? -1 : 1);
    int direction = center + nearest * scale <= 0 ? 1 : -1;
    for (Integer step = 0; ; ++step) {
        // Enumerate in order of increasing absolute projected coordinate.
        Integer z = nearest + (step % 2 ? 1 : -1) * ((step+1)/2) * direction;
        assert(std::abs(z) < 1000000 && tail_norm < 1000000);
        Integer approx = std::abs(center + z * scale);
        // Each rounded Gram coefficient, and the rounded target, errs by at most 1.
        Integer lower = std::max<Integer>(0, approx - 1 - tail_norm);
        assert(lower < 1000000000);
        // scale^2 * norm_scale / budget_scale = 2^27.
        Wide cost = (Wide(lower) * lower * norm_lower[i]) >> 27;
        if (cost > budget) break;
        ++nodes;
        coeff[i] = z;
        if (z) for (int h = 0; h < i; ++h) offset[h] += z * gram_coeff[i][h];
        search(i-1, budget-Integer(cost), tail_norm + std::abs(z));
        if (z) for (int h = 0; h < i; ++h) offset[h] -= z * gram_coeff[i][h];
    }
}

int main(int argc, char** argv) {
    int pairs;
    std::cin >> rank >> dimension >> degree >> scale >> budget_scale >> norm_scale >> pairs;
    assert((rank == 47 && dimension == 48 && degree == 60 && pairs == 50) ||
           (rank == 90 && dimension == 91 && degree == 117 && pairs == 326));
    assert(scale == 32768 && budget_scale == scale && norm_scale == 4096);
    for (int i = 0; i < rank; ++i) for (int j = 0; j < dimension; ++j) {
        std::cin >> basis[i][j];
        assert(std::abs(basis[i][j]) < 1000000);
    }
    for (int i = 0; i < rank; ++i) for (int j = 0; j < rank; ++j) {
        std::cin >> gram_coeff[i][j];
        assert(std::abs(gram_coeff[i][j]) < 1000000);
    }
    for (int i = 0; i < rank; ++i) {
        std::cin >> norm_lower[i];
        assert(norm_lower[i] > 0 && norm_lower[i] < 1000000000000LL);
    }
    for (int k = 1; k <= pairs; ++k) {
        int a, b;
        std::cin >> a >> b;
        for (int h = 0; h < rank; ++h) {
            std::cin >> offset[h];
            assert(std::abs(offset[h]) < 1000000);
            offset[h] = -offset[h];
        }
        assert(std::cin);
        if (argc > 1 && k < std::atoi(argv[1])) continue;
        if (argc > 2 && k > std::atoi(argv[2])) continue;
        search(rank-1, Integer(degree)*degree*budget_scale/2, 0);
        std::cout << k << ' ' << a << ' ' << b << " nodes=" << nodes
                  << " leaves=" << leaves << std::endl;
    }
}
