#include <complex>
#include <vector>
#include <cassert>
#include <algorithm>

struct polynom_n {
typedef std::complex<double> com;

    void check() {
        for (auto &x : a)
            assert(x.size() == a.size());
    }

    polynom_n(const std::vector<std::vector<com>> &a) : a(a) {
        check();
    }
    polynom_n(const polynom_n &other) : a(other.a) {}
    polynom_n& operator=(const polynom_n &other) {
        a = other.a;
        return *this;
    }

    size_t size() const {
        return a.size();
    }

    polynom_n& operator*=(const polynom_n &other) {
        assert(size() == other.size());
        size_t n = a.size();
        std::vector<std::vector<com>> temp(n << 1);
        std::generate(temp.begin(), temp.end(), [n]() {return std::vector<com>(n << 1);});

        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                for (size_t ii = 0; ii < n; ++ii)
                    for (size_t jj = 0; jj < n; ++jj)
                        temp[i + ii][j + jj] += a[i][j] * other.a[i][j];
        a = temp;
        return *this;
    }
    friend polynom_n operator* (polynom_n a, const polynom_n &b);

    com get(int i, int j) const {
        return a[i][j];
    }

    std::vector<std::vector<com>> a;
};

polynom_n operator*(polynom_n a, const polynom_n &b) {
    return a *= b;
}
