#ifndef UTILLITY_H
#define UTILLITY_H

#include <vector>
#include <string_view>
#include <numbers>
#include <cmath>
#include <iostream>
#include <complex>

template<typename T>
void printVector(std::string_view s, std::vector<T> v) {
	std::cout << s;
	for (auto& ele : v) {
		if constexpr (std::is_same_v<T, std::complex<double>>) {
            std::cout << ele << " ";
        } else {
            std::cout << ele << " ";
        }
	}
	std::cout << "\n";
}

// Return sample points of a sin wave
std::vector<double> sampleSin(const int N, const double freq);

// Validate the result of the FFT against a known good
int validateFFT(std::vector<std::complex<double>> &X, std::vector<std::complex<double>> &Y);



#endif // UTILLITY_H