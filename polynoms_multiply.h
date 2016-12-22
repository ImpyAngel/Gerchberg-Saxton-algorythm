#ifndef POLYNOMS_MULTIPLY_H
#define POLYNOMS_MULTIPLY_H

struct polynom_n {
typedef std::complex<double> com;

    void check();

    polynom_n(const std::vector<std::vector<com>> &a);
    polynom_n(const polynom_n &other);
    polynom_n& operator=(const polynom_n &other);

    size_t size() const;

    polynom_n& operator*=(const polynom_n &other);

    friend polynom_n operator* (polynom_n a, const polynom_n &b);

    com get(int i, int j) const;

private:
    std::vector<std::vector<com>> a;
};

polynom_n operator*(polynom_n a, const polynom_n &b);

#endif // POLYNOMS_MULTIPLY_H
