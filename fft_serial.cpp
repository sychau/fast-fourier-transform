#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <numbers>
#include <string>
#include <bitset>
#include <iomanip>

#include <fftw3.h>

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
	const int nBitsLen = std::log2(n); // number of bits to represent n

	double exponentSign = isInverse ? 1.0 : -1.0;
	const std::complex<double> img(0.0, 1.0);

	 
	// Reorder the input vector using bit reversal
	auto bitReverseCopy = [](const std::vector<std::complex<double>> A) {
		const int n = A.size();
		const int nBitsLen = std::log2(n); // number of bits to represent n
		std::vector<std::complex<double>> R(n);

		// Reverse the bits of k, for example 01001 will become 10010
		auto rev = [] (int k, int bits) {
			int reversed_index = 0;
			for (int i = 0; i < bits; ++i) {
				if (k & (1 << i)) {
					reversed_index |= (1 << (bits - 1 - i));
				}
			}
			return reversed_index;
		};

		for (int i = 0; i < n; ++i) {
			R[rev(i, nBitsLen)] = A[i];
		}
		return R;
	};

	std::vector<std::complex<double>> Y = bitReverseCopy(X);

	for (int s = 1; s <= nBitsLen; ++s) {
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
int main() {
    // std::vector<double> samples = sampleSin(8, 1);
	std::vector<double> samples = {1, 2, 3, 4, 5, 6, 7, 8};
	
    std::cout << std::setprecision(3);
	printVector("samples: ", samples);

	std::vector<std::complex<double>> naiveResult = naiveDFT(samples);
	printVector("naive dft: ", naiveResult);
	
	const std::complex<double> j(0.0, 1.0);
	std::complex<double> omega = std::exp(-j * 2.0 * std::numbers::pi / (double)samples.size());
	std::vector<std::complex<double>> fftResult = recursiveFFT(samples, omega);
	printVector("recursive fft: ", fftResult);

	std::vector<std::complex<double>> complexSamples(samples.begin(), samples.end());
	std::vector<std::complex<double>> iterativeFftResult = iterativeFFT(complexSamples, false);
	printVector("iterative fft: ", iterativeFftResult);

	std::vector<std::complex<double>> iterativeIfftResult = iterativeFFT(iterativeFftResult, true);
	printVector("iterative ifft: ", iterativeIfftResult);

	std::vector<std::complex<double>> fftwResult = expandFftwResult(fftwR2c(samples));
	printVector("fftw: ", fftwResult);

	std::vector<double> ifftwResult = ifftwC2r(fftwResult);
	printVector("ifftw: ", ifftwResult);
}