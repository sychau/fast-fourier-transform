#include "utillity.h"
#include "fft_multiThreaded.h"
#include "fft_serial.h"

#include "core/cxxopts.h"

#include <chrono>
#include <thread>
#include <iomanip>

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
	std::cout << "Number of threads: "<< nThreads << std::endl;
	std::cout << "Array size: " << arraySize << std::endl;
	std::cout << "Seed: " << seed << std::endl;
	std::cout << "\n";

	if (!isPowerOfTwo(arraySize) || nThreads <= 0 || arraySize < nThreads) {
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
	std::vector<std::complex<double>> iterativeFftResult = multithreaded_iterativeIcpFft(complexSamples, false, nThreads);
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
	std::cout << "multithreaded_iterativeIcpFft ("<< nThreads <<" threads): " << duration << " milliseconds\n";
	std::cout << '\n';
	
	// Verify result
	if (validateFFT(fftwResult, iterativeFftResult) == 0) {
		std::cout << "multithreaded_iterativeIcpFft passed!\n";
	} else {
		std::cout << "multithreaded_iterativeIcpFft failed!\n";
	}

    return 0;
}