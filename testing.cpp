#include "utillity.h"
#include "fft_multiThreaded.h"
#include "fft_serial.h"

#include <chrono>
#include <thread>
#include <iomanip>
#include <numbers>


int main(int argc, char* argv[]){
    std::vector<double> samples = sampleSin(8, 1);
    // std::vector<double> samples = {1, 2, 3, 4, 5, 6, 7, 8};

    //std::cout << std::setprecision(3);
	printVector("samples: ", samples);

	std::vector<std::complex<double>> naiveResult = naiveDFT(samples);
	printVector("naive dft: ", naiveResult);
	
	const std::complex<double> j(0.0, 1.0);
	std::complex<double> omega = std::exp(-j * 2.0 * std::numbers::pi / (double)samples.size());
	std::vector<std::complex<double>> fftResult = recursiveFFT(samples, omega);
	printVector("recursive fft: ", fftResult);

	std::vector<std::complex<double>> complexSamples(samples.begin(), samples.end());
	std::vector<std::complex<double>> iterativeFftResult = iterativeFFT(complexSamples, false);
	printVector("iterative fft: ", iterativeFftResult);

	std::vector<std::complex<double>> iterativeIfftResult = iterativeFFT(iterativeFftResult, true);
	printVector("iterative ifft: ", iterativeIfftResult);

	std::vector<std::complex<double>> fftwResult = expandFftwResult(fftwR2c(samples));
	printVector("fftw: ", fftwResult);

	std::vector<double> ifftwResult = ifftwC2r(fftwResult);
	printVector("ifftw: ", ifftwResult);

    std::vector<std::complex<double>> naiveResult_multiThreaded = naiveDFT_multiThreaded(samples);
    printVector("naive dft multiThreaded: ", naiveResult_multiThreaded);

    if (validateFFT(fftwResult, naiveResult_multiThreaded) == 0) {
        std::cout << "naiveDFT_multiThreaded passed\n";
    } else {
        std::cout << "naiveDFT_multiThreaded failed\n";
    }

    return 0;
}