#ifndef FFT_DISTRIBUTED
#define FFT_DISTRIBUTED

#include <vector>
#include <complex>
#include <mpi.h>

std::vector<std::complex<double>> iterativeFFTDistributed(const std::vector<std::complex<double>> &X, bool isInverse);

#endif // FFT_DISTRIBUTED