#include "fft_distributed.h"

std::vector<std::complex<double>> iterativeFFTDistributed(const std::vector<std::complex<double>> &X, bool isInverse) {
    
    int worldSize;
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    int worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    std::cout << "Number of processes: " <<  worldSize << "\n";
    std::cout << "Current process: " << worldRank << "\n";

    return {};
}