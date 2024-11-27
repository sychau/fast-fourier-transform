#include "fft_distributed.h"
#include "fft_serial.h"
#include "utillity.h"
#include "core/get_time.h"
#include "core/cxxopts.h"

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
        if (!isPowerOfTwo(arraySize) || !isPowerOfTwo(processes) || arraySize < processes) {
            std::cout << "Incorrect input size or number of processes\n";
            return 1;
        }
    }

    // Generate random samples
    std::vector<double> samples = generateRandomVector(arraySize, 0.0, 10.0, seed);
    std::vector<std::complex<double>> complexSamples(samples.begin(), samples.end());

    // Perform FFT and time it
    timer fftTimer;
    fftTimer.start();
	std::vector<std::complex<double>> iterFftMpiRes = icpFftDistributed(complexSamples, false);
    double fftTime = fftTimer.stop();

	std::vector<std::complex<double>> iterIfftMpiRes = icpFftDistributed(iterFftMpiRes, true);

    // Verify result
    if (currProcId == 0) {
        timer fftwTimer;
        fftwTimer.start();
        std::vector<std::complex<double>> fftwResult = expandFftwResult(fftwR2c(samples));
        double fftwTime = fftwTimer.stop();

        std::cout << std::format("Time taken(fftw) ms: {}\n", fftwTime * 1000);
        std::cout << std::format("Time taken(icp fft) ms: {}\n", fftTime * 1000);
        std::cout << "\n";

        if (validateFFT(fftwResult, iterFftMpiRes) == 0) {
            std::cout << "iterative fft distributed passed!\n";
        } else {
            std::cout << "iterative fft distributed failed!\n";
        }

        std::vector<double> ifftwResult = ifftwC2r(fftwResult);
        std::vector<std::complex<double>> complexIfftwResult(ifftwResult.begin(), ifftwResult.end());

        if (validateFFT(complexIfftwResult, iterIfftMpiRes) == 0) {
            std::cout << "iterative ifft distributed passed!\n";
        } else {
            std::cout << "iterative ifft distributed failed!\n";
        }
    }
	MPI_Finalize();
}
