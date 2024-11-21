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

/*
// Return sample points of a sin wave
std::vector<double> sampleSin(const int N, const double freq) {
    double samplingRate = N;

    std::vector<double> samples(N);    
    for (int i = 0; i < N; ++i) {
        double t = i / samplingRate;
        samples[i] = std::sin(2 * std::numbers::pi * freq * t);
    }
    return samples;
}

*/
#endif // UTILLITY_H