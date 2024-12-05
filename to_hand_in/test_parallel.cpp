#include "utillity.h"
#include "fft_multiThreaded.h"
#include "fft_serial.h"

#include "core/cxxopts.h"

#include <chrono>
#include <thread>
#include <iomanip>
#include <numbers>

#define FUNC_NAME(x) x, #x
#define DEFAULT_THREADS "1"
#define DEFAULT_ARRAY_SIZE_EXP "4"
#define DEFAULT_SEED "0"

int main(int argc, char* argv[]){
	// Initialize command line arguments
    cxxopts::Options options("Fast Fourier Transfer Multi-threaded", "Perform Fast Fourier Transfer on a 1D array");
    options.add_options(
        "custom",
        {	
			{"nThreads", "Number of threads",
			cxxopts::value<uint>()->default_value(DEFAULT_THREADS)},
            {"sizeExp", "Array Size Exponent (Max 31)",         
            cxxopts::value<uint>()->default_value(DEFAULT_ARRAY_SIZE_EXP)},
            {"seed", "Random seed",
            cxxopts::value<uint>()->default_value(DEFAULT_SEED)},
        }
    );
    auto cl_options = options.parse(argc, argv);
	uint nThreads = cl_options["nThreads"].as<uint>();
    uint arraySizeExp = cl_options["sizeExp"].as<uint>();
    uint arraySize = 1 << arraySizeExp;
    uint seed = cl_options["seed"].as<uint>();

    // Print information
	std::cout << std::format("Number of threads: {}\n", nThreads);
	std::cout << std::format("Array size: {}\n", arraySize);
	std::cout << std::format("Seed: {}\n", seed);
	std::cout << "\n";

	if (!isPowerOfTwo(arraySize) || nThreads <= 0 || arraySize < nThreads) {
		std::cout << "Incorrect input size or number of threads\n";
		return 1;
	}
    // Generate random samples
    std::vector<double> samples = generateRandomVector(arraySize, 0.0, 10.0, seed);
    std::vector<std::complex<double>> complexSamples(samples.begin(), samples.end());
	

	// IterativeIcpFFT and time it as the baseline serial
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::vector<std::complex<double>> iterativeIcpFftResult = iterativeIcpFft(complexSamples, false);
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "iterativeIcpFft: " << duration << " microseconds\n";

	// FFTW for verification and time comparison
	std::vector<std::complex<double>> fftwResult = expandFftwResult(fftwR2c(samples));

	// Multithreaded iterativeICPFFT and time it
	t1 = std::chrono::high_resolution_clock::now();
	std::vector<std::complex<double>> multithreadedIterativeIcpFftResult = multithreaded_iterativeIcpFft(complexSamples, false, nThreads);
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "multithreaded_iterativeIcpFft ("<< nThreads <<" threads): " << duration << " microseconds\n";
	std::cout << "\n";

	// Verify result
	if (validateFFT(fftwResult, multithreadedIterativeIcpFftResult) == 0) {
		std::cout << "multithreaded_iterativeIcpFft passed!\n";
	} else {
		std::cout << "multithreaded_iterativeIcpFft failed!\n";
	}

	// Verify result
	if (validateFFT(fftwResult, iterativeIcpFftResult) == 0) {
		std::cout << "iterativeIcpFft passed!\n";
	} else {
		std::cout << "iterativeIcpFft failed!\n";
	}


    return 0;
}