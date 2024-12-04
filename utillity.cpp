#include "utillity.h"

#define VALIDATE_EPSILON "1e-12"



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

// Validate the result of the FFT against a known good
int validateFFT(std::vector<std::complex<double>> &X, std::vector<std::complex<double>> &Y){
    // Return 0 if the two vectors are equal, 1 otherwise

    // Must be the same size. 
    if (X.size() != Y.size()) {
        return 1;
    }

    // Compare each element
    for (int i = 0; i < X.size(); ++i) {
        if (std::abs(X[i].real() - Y[i].real()) > 1e-4 || 
            std::abs(X[i].imag() - Y[i].imag()) > 1e-4) {
            std::cout <<  "X:" << X[i] << " Y: " << Y[i] << "\n";
            std::cout << "Difference of " << std::abs(X[i] - Y[i]) << " at index " << i << "\n";
            return 1;
        }
    }

    return 0;
}



void autoValidate(std::vector<std::complex<double>> (*func)(const std::vector<double>&), const char* func_name ,std::vector<double> &X, std::vector<std::complex<double>> &knownGood){
    std::vector<std::complex<double>> Y = func(X);
    if (validateFFT(knownGood, Y) == 0) {
        std::cout<< func_name << " passed\n";
    } else {
        std::cout << func_name<< " failed\n";
    }
}
void autoValidate(std::vector<std::complex<double>> (*func)(const std::vector<double>&, std::complex<double>), const char* func_name ,std::vector<double> &X, std::complex<double> omega, std::vector<std::complex<double>> &knownGood){
    std::vector<std::complex<double>> Y = func(X, omega);
    if (validateFFT(knownGood, Y) == 0) {
        std::cout<< func_name << " passed\n";
    } else {
        std::cout << func_name<< " failed\n";
    }
}

void autoValidate(std::vector<std::complex<double>> (*func)(const std::vector<std::complex<double>>&, bool), const char* func_name ,std::vector<std::complex<double>> &X, bool inverse, std::vector<std::complex<double>> &knownGood){
    std::vector<std::complex<double>> Y = func(X, inverse);
    if (validateFFT(knownGood, Y) == 0) {
        std::cout<< func_name << " passed\n";
    } else {
        std::cout << func_name<< " failed\n";
    }
}
void autoValidate(std::vector<std::complex<double>> (*func)(const std::vector<std::complex<double>>&, std::complex<double>), const char* func_name ,std::vector<std::complex<double>> &X, std::complex<double> omega, std::vector<std::complex<double>> &knownGood){
    std::vector<std::complex<double>> Y = func(X, omega);
    if (validateFFT(knownGood, Y) == 0) {
        std::cout<< func_name << " passed\n";
    } else {
        std::cout << func_name<< " failed\n";
    }
}



std::vector<std::complex<double>> butterflyAdd(std::complex<double> a, std::complex<double> b){
    // Return {a + b, a - b}
    return {a + b, a - b};
}

std::vector<double> generateRandomVector(size_t size, double min, double max, unsigned int seed) {
    std::vector<double> result(size);
    std::mt19937 gen(seed); // Initialize generator with a seed
    std::uniform_real_distribution<double> dis(min, max);

    for (double& val : result) {
        val = dis(gen);
    }
    return result;
}






unsigned int reverseBits(unsigned int num, int len) {
    // Reverse bits of num
	unsigned int reversed = 0;
	for (int i = 0; i < len; ++i) {
		unsigned int bit = (num >> i) & 1;
		reversed = (reversed << 1) | bit;
	}
	return reversed;
};




unsigned int getFirstNBits(unsigned int num, int n) {
    // Extract the first N bits
	unsigned int mask = (1U << n) - 1;
	return num & mask;  
};