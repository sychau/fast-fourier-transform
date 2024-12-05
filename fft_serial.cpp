#include "fft_serial.h"
#include "utillity.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <numbers>
#include <string>
#include <bitset>
#include <iomanip>
#include <format>
#include <fftw3.h>
#include <map>

// IPC Ch.13
// O(n^2)
std::vector<std::complex<double>> naiveDFT(const std::vector<double> &X) {
    const int n = X.size();
    std::vector<std::complex<double>> Y(X.size());

	const std::complex<double> img(0.0, 1.0);
	std::complex<double> omega = std::exp(-2.0 * img * std::numbers::pi / (double)n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Y[i] += X[j] * std::pow(omega, i * j);
		}
    }
	return Y;
}

// IPC Ch.13
// O(n*log(n)) recursive
std::vector<std::complex<double>> recursiveFFT(const std::vector<double> &X, std::complex<double> omega) {
	std::vector<std::complex<double>> Y(X.size());
	const int n = X.size();
	if (n == 1) {
		Y[0] = X[0];
		return Y;
	}

	// Populate new vectors containing only odd or even indexed elements
	std::vector<double> xEven;
    for (int i = 0; i < n; i += 2) {
        xEven.push_back(X[i]);
    }
	std::vector<double> xOdd;
    for (int i = 1; i < n; i += 2) {
        xOdd.push_back(X[i]);
    }
	
	std::vector<std::complex<double>> left = recursiveFFT(xEven, std::pow(omega, 2));
	std::vector<std::complex<double>> right = recursiveFFT(xOdd, std::pow(omega, 2));
	for (int i = 0; i < n; ++i) {
		Y[i] = left[i % (n / 2)] + std::pow(omega, i) * right[i % (n / 2)];
	}
	return Y;
}

// https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
// O(n*log(n)) iterative
// If isInverse is true, IFFT will be performed
// Input must be length of power of 2, complex number
std::vector<std::complex<double>> iterativeFFT(const std::vector<std::complex<double>> &X, bool isInverse) {
	const int n = X.size();
	const int r = std::log2(n); // number of bits to represent n

	double exponentSign = isInverse ? 1.0 : -1.0;
	const std::complex<double> img(0.0, 1.0);

	// Use bit reversal order
	std::vector<std::complex<double>> Y(n);
	for (int i = 0; i < n; ++i) {
		Y[reverseBits(i, r)] = X[i];
	}

	for (int s = 1; s <= r; ++s) {
		int m = 1 << s;
		std::complex<double> omegaM = std::exp(exponentSign * 2.0 * img * std::numbers::pi / (double)m);

		for (int k = 0; k < n; k += m) {
			std::complex<double> omega = 1.0;

			for (int j = 0; j < m / 2; ++j) {
				std::complex<double> t = omega * Y[k + j + m / 2];
				std::complex<double> u = Y[k + j];
				Y[k + j] = u + t;
				Y[k + j + m / 2] = u - t;
				omega *= omegaM;
			}
		}
	}
	// Normalize result if performing IFFT
    if (isInverse) {
        for (int i = 0; i < n; ++i) {
            Y[i] /= n;
        }
    }
	
	return Y;
}

// This iterative version is from ICP, it is easier to decompose tasks to different threads/ processes
std::vector<std::complex<double>> iterativeIcpFft(const std::vector<std::complex<double>> &X, bool isInverse) {
	const int n = X.size();
	const int r = std::log2(n);

	const std::complex<double> img(0.0, 1.0);
	double exponentSign = isInverse ? 1.0 : -1.0;
	std::complex<double> omegaBase = std::exp(exponentSign * 2.0 * img * std::numbers::pi / (double)n);

	std::vector<std::complex<double>> R(X); // Result array
    std::vector<std::complex<double>> S(R); // Auxillary array to hold previous value of R

	std::map<int, std::complex<double>> omegaCache; // Store previously calculated omega

	// Outer loop O(log n), m represent stage
    for (int m = 0; m < r; ++m) {
		// Inner loop O(n)
		for (int i = 0; i < n; ++i) {
			// Let (b_0 b_1 ... b_r-1) be the binary representation of i 
            unsigned int j = i & ~(1 << (r - 1 - m)); // Clear the bit, j:= (b_0 ... b_m−1 0 b_m+1 ... b_r−1)
            unsigned int k = i | (1 << (r - 1 - m));  // Set the bit, k:= (b_0 ... b_m−1 1 b_m+1 ... b_r−1)

			unsigned int omegaExp = getFirstNBits(reverseBits(i, r), m + 1) << (r - 1 - m); // (b_m b_m-1 ... b-0 ... 0 ... 0)
			std::complex<double> omega;
            if (omegaCache.find(omegaExp) != omegaCache.end()) {
                omega = omegaCache[omegaExp];
            } else {
                omega = std::pow(omegaBase, omegaExp);
                omegaCache[omegaExp] = omega;
            }
            R[i] = S[j] + S[k] * omega;
		}
		// Update S with new values
		std::copy(R.begin(), R.end(), S.begin());
    }

    // In-place bit-reversal reordering
    for (int i = 0; i < n; ++i) {
        int reversedIndex = reverseBits(i, r);
        if (i < reversedIndex) {
            std::swap(R[i], R[reversedIndex]);
        }
    }

	// Normalize result if performing IFFT
    if (isInverse) {
        for (int i = 0; i < n; ++i) {
            R[i] /= n;
        }
    }

    return R;
}

// Using fftw library
std::vector<std::complex<double>> fftwR2c(std::vector<double> &x) {
    int N = x.size();
    std::vector<std::complex<double>> y(N / 2 + 1);
    fftw_plan p = fftw_plan_dft_r2c_1d(N, const_cast<double*>(x.data()),
                                       reinterpret_cast<fftw_complex*>(y.data()), FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    return y;
}

std::vector<double> ifftwC2r(std::vector<std::complex<double>> &y) {
	// The input gets changed sometimes.
	int N = y.size();
    std::vector<double> x(N);
    fftw_plan p = fftw_plan_dft_c2r_1d(N, reinterpret_cast<fftw_complex*>(y.data()),
                                       x.data(), FFTW_ESTIMATE);

    fftw_execute(p);
    fftw_destroy_plan(p);

    // Normalize the result
    for (int i = 0; i < N; ++i) {
        x[i] /= N;
    }

    return x;
}

// fftw exploits symmetry to save spaces, it returns a vector that is half the size of the input
// This function recover the orignal size
std::vector<std::complex<double>> expandFftwResult(
	const std::vector<std::complex<double>>& result
) {
    int N = (result.size() - 1) * 2;  // Calculate the original size of the input
    std::vector<std::complex<double>> fullResult(N);

    // Copy the real-to-complex result to the first half
    for (int i = 0; i <= N / 2; ++i) {
        fullResult[i] = result[i];
    }
    // Fill the second half with complex conjugates in reverse order
    for (int i = N / 2 + 1; i < N; ++i) {
        fullResult[i] = std::conj(result[N - i]);
    }

    return fullResult;
}
