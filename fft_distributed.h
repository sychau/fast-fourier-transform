#ifndef FFT_DISTRIBUTED
#define FFT_DISTRIBUTED

#include <vector>
#include <complex>
#include <bitset>
#include <mpi.h>
#include <format>
#include "utillity.h"

std::vector<std::complex<double>> icpFftDistributed(const std::vector<std::complex<double>> &X, bool isInverse);
bool isPowerOfTwo(unsigned int i);
unsigned int reverseBits (unsigned int num, int len);
unsigned int getFirstNBits (unsigned int num, int n);

#endif // FFT_DISTRIBUTED