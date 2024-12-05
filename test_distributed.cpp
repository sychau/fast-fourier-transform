#include "fft_distributed.h"
#include "fft_serial.h"
#include "utillity.h"

#include "core/cxxopts.h"

#include <chrono>

#define DEFAULT_ARRAY_SIZE_EXP "4"
#define DEFAULT_SEED "0"

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);

    // Initialize command line arguments
    cxxopts::Options options("Fast Fourier Transfer Distributed", "Perform Fast Fourier Transfer on a 1D array");
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

    int processes, currProcId;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);  // Total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &currProcId);  // Current process ID

    // Print information
    if (currProcId == 0) {
        std::cout << std::format("Number of processes: {}\n", processes);
        std::cout << std::format("Array size: {}\n", arraySize);
        std::cout << std::format("Seed: {}\n", seed);
        std::cout << "\n";
        if (!isPowerOfTwo(processes) || arraySize < processes) {
            std::cout << "Incorrect input size or number of processes\n";
            return 1;
        }
    }

    // Generate random samples
    std::vector<double> samples = generateRandomVector(arraySize, 0.0, 10.0, seed);
    std::vector<std::complex<double>> complexSamples(samples.begin(), samples.end());

    // Perform FFT and time it

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::vector<std::complex<double>> fftwResult = expandFftwResult(fftwR2c(samples));
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    if (currProcId == 0) {
	    std::cout << "fftw: " << duration << " milliseconds\n";
    }
    t1 = std::chrono::high_resolution_clock::now();
	std::vector<std::complex<double>> iterFftMpiRes = icpFftDistributed(complexSamples, false);
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    if (currProcId == 0) {
        std::cout << "icpFftDistributed: " << duration << " milliseconds\n";
        std::cout << "\n";
    }

    // Verify result
    if (currProcId == 0) {
        if (validateFFT(fftwResult, iterFftMpiRes) == 0) {
            std::cout << "icpFftDistributed passed!\n";
        } else {
            std::cout << "icpFftDistributed failed!\n";
        }
    }
	MPI_Finalize();
}
