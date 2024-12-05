#ifndef FFT_SERIAL_H
#define FFT_SERIAL_H

#include <vector>
#include <complex>

std::vector<std::complex<double>> naiveDFT(const std::vector<double> &X);
std::vector<std::complex<double>> recursiveFFT(const std::vector<double> &X, std::complex<double> omega);
std::vector<std::complex<double>> iterativeFFT(const std::vector<std::complex<double>> &X, bool isInverse);
std::vector<std::complex<double>> iterativeIcpFft(const std::vector<std::complex<double>> &X, bool isInverse);
std::vector<std::complex<double>> fftwR2c(std::vector<double> &x);
std::vector<double> ifftwC2r(std::vector<std::complex<double>> &y);
std::vector<std::complex<double>> expandFftwResult(const std::vector<std::complex<double>> &result);

#endif // FFT_SERIAL_H