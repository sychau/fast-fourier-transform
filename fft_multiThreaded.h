#ifndef FFT_MULTITHREADED
#define FFT_MULTITHREADED

#include <vector>
#include <complex>

#include "core/utils.h"

void dft_thread_function(const std::vector<double> &X, std::vector<std::complex<double>> &Y, const int start_i, const int end_i, const int n, const std::complex<double> omega);
std::vector<std::complex<double>> naiveDFT_multiThreaded(const std::vector<double> &X);



void multithreaded_iterativeIcpFft_loop_function(const int r, const int n, std::vector<std::complex<double>> &R, std::vector<std::complex<double>> &S, const std::complex<double> omega, CustomBarrier &b, int thread_id, int block_size = 1 << 5);
std::vector<std::complex<double>> multithreaded_iterativeIcpFft(const std::vector<std::complex<double>> &X, bool isInverse, int nThreads, int block_size= 1 << 5);


void multithreaded_iterativeFft_loop_function(const int r, double exponentSign, const int n, std::vector<std::complex<double>> &Y, CustomBarrier &b, int thread_id);
std::vector<std::complex<double>> multithreaded_iterativeFFT(const std::vector<std::complex<double>> &X, bool isInverse, int nThreads);

#endif // FFT_MULTITHREADED