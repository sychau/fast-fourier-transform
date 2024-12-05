#ifndef UTILLITY_H
#define UTILLITY_H

#include <vector>
#include <string_view>
#include <numbers>
#include <cmath>
#include <iostream>
#include <complex>
#include <random>
#include <bitset>

// Return sample points of a sin wave
std::vector<double> sampleSin(const int N, const double freq);

// Validate the result of the FFT against a known good
int validateFFT(std::vector<std::complex<double>> &X, std::vector<std::complex<double>> &Y);

//Auto validate function
void autoValidate(std::vector<std::complex<double>> (*func)(const std::vector<double>&), const char* func_name ,std::vector<double> &X, std::vector<std::complex<double>> &knownGood);
void autoValidate(std::vector<std::complex<double>> (*func)(const std::vector<double>&, std::complex<double>), const char* func_name ,std::vector<double> &X, std::complex<double> omega,  std::vector<std::complex<double>> &knownGood);
void autoValidate(std::vector<std::complex<double>> (*func)(const std::vector<std::complex<double>>&, bool), const char* func_name ,std::vector<std::complex<double>> &X, bool inverse,  std::vector<std::complex<double>> &knownGood);
void autoValidate(std::vector<std::complex<double>> (*func)(const std::vector<std::complex<double>>&, std::complex<double>), const char* func_name ,std::vector<std::complex<double>> &X, std::complex<double> omega, std::vector<std::complex<double>> &knownGood);

std::vector<std::complex<double>> butterflyAdd(std::complex<double> a, std::complex<double> b);

std::vector<double> generateRandomVector(size_t size, double min, double max, unsigned int seed);

// Check if a number is power of 2
bool isPowerOfTwo(unsigned int i);

// Reverse bits of num
unsigned int reverseBits (unsigned int num, unsigned int len);

// Extract the first N bits
unsigned int getFirstNBits (unsigned int num, int n);

#endif // UTILLITY_H