#include <iostream>
#include <cstdio>
#include <vector>
#include <complex>
#include <cassert>

typedef std::complex<double> com;

const double M_PI = 3.14159265358979323846;
const com I(0,1);
const int div = 100;
const double lambda = 1;
const double k = 2 * M_PI/lambda;
const double z = 1;
const int magic = 100;

struct discrete_fourier_transform {
    enum invert { eFFT = 1, eIFT = -1 };
private:
    discrete_fourier_transform() {}

    static void one_dimensional_fft(std::vector<com> &a, invert flag) {
        size_t n = a.size(), n2 = n >> 1;
        if (n == 1)
            return;
        assert((n & 1) == 0);
        std::vector<com> a0, a1;
        for (size_t i = 0; i < n; ++i) {
            std::vector<com> &cur = (i == 0 ? a0 : a1);
            cur.emplace_back(a[i]);
        }

        one_dimensional_fft(a0, flag);
        one_dimensional_fft(a1, flag);
        double ang = M_PI * 2 / n * flag;
        com w(1), wn(cos(ang), sin(ang));
        for (size_t i = 0; i < n2; ++i) {
            a[i] = a0[i] + w * a1[i];
            a[i + n2] = a0[i] - w * a1[i];
            w *= wn;
            if (flag == eIFT) {
                a[i] /= 2;
                a[i + n2] /= 2;
            }
        }
    }

public:
    /// a.size() should be power of 2: 2^x for unsigned x.
    /// flag should eFFT for FFT or eIFT for IFT.
    /// a.size() == a[i].size for all i: 0 <= i < a.size().
    static std::vector<std::vector<com>> ft(std::vector<std::vector<com>> a, invert flag) {
        assert(__builtin_popcount(a.size()) == 1);
        for (size_t t = 0; t < 2; ++t) {
            for (size_t i = 0; i < a.size(); ++i) {
                assert(a.size() == a[i].size());
                one_dimensional_fft(a[i], flag);
            }
            for (size_t i = 0; i < a.size(); ++i) {
                for (size_t j = i + 1; j < a.size(); ++j) {
                    swap(a[i][j], a[j][i]);
                }
            }
        }
        return a;
    }
};


std::vector<std::vector<com> > retrieved(div), A0(div);
std::vector<com> temp(div);

com sqr(com a) {
    return a * a;
}

void init(){
	std::vector<com> temp1(div, com(1, 0));
	std::vector<com> temp2(div);

	for (int i = 0; i < div/ 3; ++i)    {
		temp2[i] = com(0, 0);
		temp2[i + div/3] = com(1, 0);
		temp2[i + 2*(div/3)] = com(0, 0);
	}

	retrieved.assign(div, temp2);

	for (int i = 0; i < div/3; ++i) {
		retrieved[i + div/3] = temp1;
	}

	A0.assign(div,  temp1);
}
com H(int i, int j) {
	return exp(i * k / 2 * z *(sqr(i) + sqr(j)));
}


void Gerchberg_Saxton() {
	std::vector<std::vector<com> >fi(div,temp), W(div,temp), F(div,temp), F1(div,temp),
	W_else(div, temp);
	for (int i = 0; i < div; ++i) {
		for (int j = 0; j < div; ++j) {
			W[i][j] = A0[i][j] * exp(I* fi[i][j]);
		}
	}

	for (int k = 0; k < magic; ++k) {
		// F = FFT(W);
		for (int i = 0; i < div; ++i) {
			for (int j = 0; j < div; ++j) {
				F1[i][j] = retrieved[i][j] * F[i][j] / abs(F[i][j]);
			}
		}
		for (int i = 0; i < div; ++i) {
			for (int j = 0; j < div; ++j) {
				W_else[i][j] = F1[i][j] * conj(exp(I* fi[i][j]));
			}
		}
		// W = IFT(W_else);
		for (int i = 0; i < div; ++i) {
			for (int j = 0; j < div; ++j) {
				W[i][j] = A0[i][j] * W[i][j] / abs(W[i][j]);
			}
		}
	}
}

int main() {
	init();
	Gerchberg_Saxton();
}
