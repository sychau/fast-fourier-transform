#include "utillity.h"
#include "fft_multiThreaded.h"
#include "fft_serial.h"

#include "core/cxxopts.h"

#include <chrono>
#include <thread>
#include <iomanip>
#include <numbers>

#define DEFAULT_ARRAY_SIZE_EXP "4"
#define DEFAULT_SEED "0"

int main(int argc, char* argv[]){
	// Initialize command line arguments
    cxxopts::Options options("Fast Fourier Transfer Serial", "Perform Fast Fourier Transfer on a 1D array");
    options.add_options(
        "custom",
        {	
            {"sizeExp", "Array Size Exponent (Max 31)",         
            cxxopts::value<uint>()->default_value(DEFAULT_ARRAY_SIZE_EXP)},
            {"seed", "Random seed",
            cxxopts::value<uint>()->default_value(DEFAULT_SEED)},
        }
    );
    auto cl_options = options.parse(argc, argv);
    uint arraySizeExp = cl_options["sizeExp"].as<uint>();
    uint arraySize = 1 << arraySizeExp;
    uint seed = cl_options["seed"].as<uint>();

    // Print information
	std::cout << "Array size: " <<  arraySize << std::endl;
	std::cout << "Seed: " << seed << std::endl;
	std::cout << "\n";

	if (!isPowerOfTwo(arraySize)) {
		std::cout << "Incorrect input size or number of threads\n";
		return 1;
	}
    // Generate random samples
    std::vector<double> samples = generateRandomVector(arraySize, 0.0, 10.0, seed);
    std::vector<std::complex<double>> complexSamples(samples.begin(), samples.end());
	
	// Perform FFT and time it
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::vector<std::complex<double>> fftwResult = expandFftwResult(fftwR2c(samples));
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
	std::cout << "fftw: " << duration << " milliseconds\n";

	t1 = std::chrono::high_resolution_clock::now();
    const std::complex<double> img(0.0, 1.0);
	std::complex<double> omega = std::exp(-2.0 * img * std::numbers::pi / (double)arraySize);
	std::vector<std::complex<double>> recursiveResult = recursiveFFT(samples, omega);
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
	std::cout << "recursiveFFT: " << duration << " milliseconds\n";
	
	t1 = std::chrono::high_resolution_clock::now();
	std::vector<std::complex<double>> iterativeResult = iterativeFFT(complexSamples, false);
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
	std::cout << "iterativeFFT: " << duration << " milliseconds\n";

	t1 = std::chrono::high_resolution_clock::now();
	std::vector<std::complex<double>> iterativeIcpResult = iterativeIcpFft(complexSamples, false);
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
	std::cout << "iterativeIcpFft: " << duration << " milliseconds\n";
    std::cout << '\n';

	// Verify result
	if (validateFFT(fftwResult, recursiveResult) == 0) {
		std::cout << "recursiveFFT passed!\n";
	} else {
		std::cout << "recursiveFFT failed!\n";
	}
	if (validateFFT(fftwResult, iterativeResult) == 0) {
		std::cout << "iterativeFFT passed!\n";
	} else {
		std::cout << "iterativeFFT failed!\n";
	}
    if (validateFFT(fftwResult, iterativeIcpResult) == 0) {
		std::cout << "iterativeIcpFft passed!\n";
	} else {
		std::cout << "iterativeIcpFft failed!\n";
	}
    return 0;
}