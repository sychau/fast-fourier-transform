#include "fft_multiThreaded.h"
#include <thread>
#include <iostream>
#include <numbers>

void dft_thread_function(const std::vector<double> &X, std::vector<std::complex<double>> &Y, const int start_i, const int end_i, const int n, const std::complex<double> omega) {
    
    for (int i = start_i; i < end_i; ++i) {
        for (int j = 0; j < n; ++j) {
            Y[i] += X[j] * std::pow(omega, i * j);
        }
    }
}


std::vector<std::complex<double>> naiveDFT_multiThreaded(const std::vector<double> &X){
    const int n = X.size();
    std::vector<std::complex<double>> Y(X.size());

	const std::complex<double> img(0.0, 1.0);
	std::complex<double> omega = std::exp(-2.0 * img * std::numbers::pi / (double)n);

    std::vector<std::thread> threads;
    const int nThreads = std::thread::hardware_concurrency(); // This gets the amount of availble threads.
    std::cout << "Number of threads: " << nThreads << std::endl;

    const int chunkSize = n / nThreads;
    const int remainder = n % nThreads;
    std::vector<int> end_boundaries;

    for (int j =0; j < nThreads; ++j){
        if (j ==0){
            end_boundaries.push_back(chunkSize + (j < remainder ? 1 : 0));
        } else {
            end_boundaries.push_back(end_boundaries[j-1] + chunkSize + (j < remainder ? 1 : 0));
        }
    }

    for (int i = 0; i < nThreads; ++i) {
        threads.push_back(std::thread(dft_thread_function, std::ref(X), std::ref(Y), i==0 ? 0: end_boundaries[i-1], end_boundaries[i], n, omega));
    }

    for (auto &t : threads) {
        t.join();
    }

    return Y;
}


