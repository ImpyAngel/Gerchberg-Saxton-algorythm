#include <iostream>
#include <cstdio>
#include <vector>

double sqr(double a) {
    return a * a;
}

const int div = 4;
const double i_max = 100000; 
int main() {
    freopen("input.in", "w", stdout);
    for (int i = 0; i < div; ++i) {
    	for (int j = 0; j < div; ++j) {
    		double k = 1 + sqr(sqr(i - div / 2) + sqr(j - div / 2));
    		std::cout << i_max / k << ' ';
    	}
    	std::cout << '\n';
    }
}