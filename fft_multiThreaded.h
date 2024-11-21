#include <vector>

std::vector<std::complex<double>> naiveDFT_multiThreaded(const std::vector<double> &X);
std::vector<std::complex<double>> recursiveFFT_multiThreaded(const std::vector<double> &X, std::complex<double> omega);
std::vector<std::complex<double>> iterativeFFT_multiThreaded(const std::vector<std::complex<double>> &X, bool isInverse);