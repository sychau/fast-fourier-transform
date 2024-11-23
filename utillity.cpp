#include "utillity.h"

#define VALIDATE_EPSILON "1e-11"



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
        if (X[i] != Y[i] && std::abs(X[i] - Y[i]) > 1e-6) {
            std::cout << "Difference of " << std::abs(X[i] - Y[i]) << " at index " << i << "\n";
            return 1;
        }
    }

    return 0;
}
