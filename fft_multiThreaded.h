#ifndef FFT_MULTITHREADED
#define FFT_MULTITHREADED

#include <vector>
#include <complex>


void dft_thread_function(const std::vector<double> &X, std::vector<std::complex<double>> &Y, const int start_i, const int end_i, const int n, const std::complex<double> omega);
std::vector<std::complex<double>> naiveDFT_multiThreaded(const std::vector<double> &X);
//std::vector<std::complex<double>> recursiveFFT_multiThreaded(const std::vector<double> &X, std::complex<double> omega);
//std::vector<std::complex<double>> iterativeFFT_multiThreaded(const std::vector<std::complex<double>> &X, bool isInverse);



#endif // FFT_MULTITHREADED