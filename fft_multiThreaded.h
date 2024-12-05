#ifndef FFT_MULTITHREADED
#define FFT_MULTITHREADED

#include <vector>
#include <complex>
#include <barrier>


void dft_thread_function(const std::vector<double> &X, std::vector<std::complex<double>> &Y, const int start_i, const int end_i, const int n, const std::complex<double> omega);
std::vector<std::complex<double>> naiveDFT_multiThreaded(const std::vector<double> &X);



void multithreaded_iterativeIcpFft_innner_loop_function(const int r, const int n, std::vector<std::complex<double>> &R, std::vector<std::complex<double>> &S, const std::complex<double> omega, std::barrier<> &b, int thread_id, int block_size =1);
std::vector<std::complex<double>> multithreaded_iterativeIcpFft(const std::vector<std::complex<double>> &X, bool isInverse, int nThreads, int block_size=1);



#endif // FFT_MULTITHREADED