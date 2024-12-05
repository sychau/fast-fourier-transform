#ifndef FFT_DISTRIBUTED_H
#define FFT_DISTRIBUTED_H

#include "utillity.h"
#include "fft_serial.h"

#include <vector>
#include <complex>
#include <bitset>
#include <mpi.h>
#include <format>

std::vector<std::complex<double>> icpFftDistributed(const std::vector<std::complex<double>> &X, bool isInverse);

#endif // FFT_DISTRIBUTED_H