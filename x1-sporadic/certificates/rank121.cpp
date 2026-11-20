// Enumerate rational effective divisors and test their evaluations over F_128.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

using Byte = std::uint8_t;
constexpr int columns = 25;
constexpr int twists = 4;
using Row = std::array<Byte, columns>;
using Evaluations = std::array<Row, twists>;

int read_integer(std::istream& in) {
    std::array<Byte,4> bytes;
    in.read(reinterpret_cast<char*>(bytes.data()), bytes.size());
    if (!in || bytes[2] || bytes[3]) throw std::runtime_error("invalid integer");
    return bytes[0] + 256 * bytes[1];
}

struct Place {
    int degree;
    std::vector<Evaluations> jets;
};

struct Field {
    Byte product[128][128]{};
    Byte inverse[128]{};

    Field() {
        for (int i = 0; i < 128; ++i) {
            for (int j = 0; j < 128; ++j) {
                int a = i, b = j, c = 0;
                while (b) {
                    if (b & 1) c ^= a;
                    b >>= 1;
                    a <<= 1;
                    if (a & 128) a ^= 131;  // z^7 + z + 1
                }
                product[i][j] = c;
                if (c == 1) inverse[i] = j;
            }
        }
    }
};

class Divisors {
    Field field;
    std::vector<Place> places;
    std::vector<std::pair<int,int>> support;
    std::array<int, twists> survivors{};
    int degree;
    int total = 0;

    bool full_rank(int twist) const {
        std::array<Row,columns> pivots{};
        std::array<bool,columns> present{};
        for (auto [index, multiplicity] : support) {
            const Place& place = places[index];
            for (int j = 0; j < multiplicity * place.degree; ++j) {
                Row row = place.jets[j][twist];
                int col = 0;
                for (; col < columns; ++col) {
                    Byte c = row[col];
                    if (!c) continue;
                    if (!present[col]) {
                        for (int k = col; k < columns; ++k)
                            pivots[col][k] = field.product[field.inverse[c]][row[k]];
                        present[col] = true;
                        break;
                    }
                    for (int k = col+1; k < columns; ++k)
                        row[k] ^= field.product[c][pivots[col][k]];
                }
                if (col == columns) return false;  // A dependent row.
            }
        }
        return true;
    }

    void enumerate(int first, int left) {
        if (left == 0) {
            ++total;
            for (int i = 0; i < twists; ++i) {
                if (full_rank(i)) return;
                ++survivors[i];
            }
            return;
        }
        for (int i = first; i < static_cast<int>(places.size()); ++i) {
            int d = places[i].degree;
            if (d > left) break;
            for (int m = 1; m*d <= left; ++m) {
                support.emplace_back(i,m);
                enumerate(i+1, left-m*d);
                support.pop_back();
            }
        }
    }

public:
    explicit Divisors(std::istream& in) {
        if (read_integer(in) != twists || read_integer(in) != columns)
            throw std::runtime_error("invalid matrix header");
        degree = read_integer(in);
        if (degree < 16 || degree > 19) throw std::runtime_error("invalid degree");
        while (in.peek() != std::char_traits<char>::eof()) {
            int d = read_integer(in), m = read_integer(in);
            if (d < 1 || d > degree || m != degree/d || (!places.empty() && places.back().degree > d))
                throw std::runtime_error("invalid place");
            Place p{d, std::vector<Evaluations>(m*d)};
            for (auto& jet : p.jets) {
                for (Row& row : jet) {
                    in.read(reinterpret_cast<char*>(row.data()), columns);
                    if (!in) throw std::runtime_error("truncated row");
                    for (Byte c : row) if (c >= 128) throw std::runtime_error("invalid field element");
                }
            }
            places.push_back(std::move(p));
        }
    }

    void verify() {
        enumerate(0,degree);
        std::cout << "Degree " << degree << ", effective divisors " << total << ", remaining";
        for (int count : survivors) std::cout << " " << count;
        std::cout << std::endl;
        // These four totals are independently obtained by factors121.py.
        constexpr int expected[] = {319800, 639480, 1278840, 2557560};
        if (total != expected[degree-16] || survivors.back() != 0)
            throw std::runtime_error("certificate failed");
    }
};

int main(int argc, char** argv) {
    if (argc != 2) throw std::runtime_error("pass the binary evaluation file");
    std::ifstream data(argv[1], std::ios::binary);
    if (!data) throw std::runtime_error("missing input");
    Divisors(data).verify();
}
