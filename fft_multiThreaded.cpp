#include "fft_multiThreaded.h"
#include "utillity.h"
#include <atomic>
#include <bitset>
#include <iostream>
#include <thread>

#define _USE_MATH_DEFINES
#include <math.h>


void dft_thread_function(const std::vector<double> &X, std::vector<std::complex<double>> &Y, const int start_i, const int end_i, const int n,
                         const std::complex<double> omega) {

    for (int i = start_i; i < end_i; ++i) {
        for (int j = 0; j < n; ++j) {
            Y[i] += X[j] * std::pow(omega, i * j);
        }
    }
}
std::vector<std::complex<double>> naiveDFT_multiThreaded(const std::vector<double> &X) {
    const int n = X.size();
    std::vector<std::complex<double>> Y(X.size());

    const std::complex<double> img(0.0, 1.0);
    std::complex<double> omega = std::exp(-2.0 * img * M_PI / (double)n);

    std::vector<std::thread> threads;
    const int nThreads = std::thread::hardware_concurrency(); // This gets the amount of availble threads.
    std::cout << "Number of threads: " << nThreads << std::endl;

    const int chunkSize = n / nThreads;
    const int remainder = n % nThreads;
    std::vector<int> end_boundaries;

    for (int j = 0; j < nThreads; ++j) {
        if (j == 0) {
            end_boundaries.push_back(chunkSize + (j < remainder ? 1 : 0));
        } else {
            end_boundaries.push_back(end_boundaries[j - 1] + chunkSize + (j < remainder ? 1 : 0));
        }
    }

    for (int i = 0; i < nThreads; ++i) {
        threads.push_back(std::thread(dft_thread_function, std::ref(X), std::ref(Y), i == 0 ? 0 : end_boundaries[i - 1], end_boundaries[i], n, omega));
    }

    for (auto &t : threads) {
        t.join();
    }

    return Y;
}

std::atomic<int> iterativeIcpFft_covered_i(0);
void multithreaded_iterativeIcpFft_innner_loop_function(const int r, const int n, std::vector<std::complex<double>> &R,
                                                        std::vector<std::complex<double>> &S, const std::complex<double> omega, CustomBarrier &b, int thread_id, int block_size) {

    int start_i;
    int m;
    int i;
    unsigned int omegaExp;


    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto barrier1_time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    auto barrier2_time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    auto thread_0_stuff = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    auto inner_loop = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    inner_loop = thread_0_stuff= barrier1_time = barrier2_time = 0;


    // Outer loop O(log n)
    for (m = 0; m < r; ++m) {
        if (thread_id == 0) {
            t1= std::chrono::high_resolution_clock::now();
            std::copy(R.begin(), R.end(), S.begin());
            iterativeIcpFft_covered_i = 0;
            t2 = std::chrono::high_resolution_clock::now();
            thread_0_stuff += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        }


        t1 = std::chrono::high_resolution_clock::now();
        b.wait();
        t2 = std::chrono::high_resolution_clock::now();
        barrier1_time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();


        // Inner loop O(n)
        t1 = std::chrono::high_resolution_clock::now();

        start_i = iterativeIcpFft_covered_i.fetch_add(block_size);
        while (start_i < n) {
            for (i = start_i; i < std::min(start_i + block_size, n); ++i) {
                // Let (b_0 b_1 ... b_r-1) be the binary representation of i
                std::bitset<32> j(i);
                j.reset(r - 1 - m); // (b_0 ... b_m−1 0 b_m+1 ... b_r−1)
                std::bitset<32> k(i);
                k.set(r - 1 - m); // (b_0 ... b_m−1 1 b_m+1 ... b_r−1)

                omegaExp = getFirstNBits(reverseBits(i, r), m + 1) << (r - 1 - m); // (b_m b_m-1 ... b-0 ... 0 ... 0)
                R[i] = S[j.to_ulong()] + S[k.to_ulong()] * std::pow(omega, omegaExp);
            }
            start_i = iterativeIcpFft_covered_i.fetch_add(block_size);
        }

        t2 = std::chrono::high_resolution_clock::now();
        inner_loop += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();


        t1 = std::chrono::high_resolution_clock::now();
        b.wait();
        t2 = std::chrono::high_resolution_clock::now();
        barrier2_time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    }
    //std::cout << "Thread " << thread_id << " barrier1_time: " << barrier1_time << " microseconds\n";
    //std::cout << "Thread " << thread_id << " barrier2_time: " << barrier2_time << " microseconds\n";
    //std::cout << "Thread " << thread_id << " inner_loop: " << inner_loop << " microseconds\n";
    if (thread_id == 0) {
        std::cout << "Thread 0 stuff: " << thread_0_stuff << " microseconds\n";
    }
}

std::vector<std::complex<double>> multithreaded_iterativeIcpFft(const std::vector<std::complex<double>> &X, bool isInverse, int nThreads, int block_size) {
    const int n = X.size();
    const int r = std::log2(n);
    CustomBarrier b(nThreads);

    const std::complex<double> img(0.0, 1.0);
    double exponentSign = isInverse ? 1.0 : -1.0;
    std::complex<double> omega = std::exp(exponentSign * 2.0 * img * M_PI / (double)n);

    std::vector<std::complex<double>> R(X); // Result array
    std::vector<std::complex<double>> S(n); // Auxillary array to hold previous value of R


    std::vector<std::thread> threads;
    for (int i = 0; i < nThreads; ++i) {
        threads.push_back(std::thread(multithreaded_iterativeIcpFft_innner_loop_function, r, n, std::ref(R), std::ref(S), omega, std::ref(b), i, block_size));
    }

    for (auto &t : threads) {
        t.join();
    }

    // In-place bit-reversal reordering
    for (int i = 0; i < n; ++i) {
        int reversedIndex = reverseBits(i, r);
        if (i < reversedIndex) {
            std::swap(R[i], R[reversedIndex]);
        }
    }

    // Normalize result if performing IFFT
    if (isInverse) {
        for (int i = 0; i < n; ++i) {
            R[i] /= n;
        }
    }

    return R;
}