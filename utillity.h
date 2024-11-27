#ifndef UTILLITY_H
#define UTILLITY_H

#include <vector>
#include <string_view>
#include <numbers>
#include <cmath>
#include <iostream>
#include <complex>
#include <random>
#include <format>

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

//Auto validate function
void autoValidate(std::vector<std::complex<double>> (*func)(const std::vector<double>&), const char* func_name ,std::vector<double> &X, std::vector<std::complex<double>> &knownGood);
void autoValidate(std::vector<std::complex<double>> (*func)(const std::vector<double>&, std::complex<double>), const char* func_name ,std::vector<double> &X, std::complex<double> omega,  std::vector<std::complex<double>> &knownGood);
void autoValidate(std::vector<std::complex<double>> (*func)(const std::vector<std::complex<double>>&, bool), const char* func_name ,std::vector<std::complex<double>> &X, bool inverse,  std::vector<std::complex<double>> &knownGood);





std::vector<std::complex<double>> butterflyAdd(std::complex<double> a, std::complex<double> b);

std::vector<double> generateRandomVector(size_t size, double min, double max, unsigned int seed);

#endif // UTILLITY_H