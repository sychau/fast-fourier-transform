#include "fft_distributed.h"
#include "fft_serial.h"
#include "utillity.h"

#define FUNC_NAME(x) x, #x


int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);

    std::vector<double> samples = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    std::vector<std::complex<double>> complexSamples(samples.begin(), samples.end());

	std::vector<std::complex<double>> iterativeFftResultDistributed = iterativeFFTDistributed(complexSamples, false);

    int worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    if (worldRank == 0) {
        std::vector<std::complex<double>> fftwResult = expandFftwResult(fftwR2c(samples));
        printVector("fftw result: ", fftwResult);
        printVector("iterative fft distributed result: ", iterativeFftResultDistributed);
        if (validateFFT(fftwResult, iterativeFftResultDistributed) == 0) {
            std::cout << "iterative fft distributed passed!\n";
        } else {
            std::cout << "iterative fft distributed failed!\n";
        }
    }
	MPI_Finalize();
}
