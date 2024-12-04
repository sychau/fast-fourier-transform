#include "utillity.h"
#include "fft_multiThreaded.h"
#include "fft_serial.h"
#include "core/get_time.h"
#include <chrono>
#include <thread>
#include <iomanip>
#include <numbers>

#define FUNC_NAME(x) x, #x

int main(int argc, char* argv[]){
    std::vector<double> samples = sampleSin(1 << 16, 1);
    // std::vector<double> samples = {1, 2, 3, 4, 5, 6, 7, 8};

    // //std::cout << std::setprecision(3);
	// printVector("samples: ", samples);

	// std::vector<std::complex<double>> naiveResult = naiveDFT(samples);
	// printVector("naive dft: ", naiveResult);
	
	const std::complex<double> j(0.0, 1.0);
	std::complex<double> omega = std::exp(-j * 2.0 * std::numbers::pi / (double)samples.size());
	//std::vector<std::complex<double>> fftResult = recursiveFFT(samples, omega);
	// printVector("recursive fft: ", fftResult);

	// std::vector<std::complex<double>> complexSamples(samples.begin(), samples.end());
	// std::vector<std::complex<double>> iterativeFftResult = iterativeFFT(complexSamples, false);
	// printVector("iterative fft: ", iterativeFftResult);

	// std::vector<std::complex<double>> iterativeIfftResult = iterativeFFT(iterativeFftResult, true);
	// printVector("iterative ifft: ", iterativeIfftResult);

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::vector<std::complex<double>> fftwResult = expandFftwResult(fftwR2c(samples));
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "fftw time: " << duration << " microseconds\n";
	//printVector("fftw: ", fftwResult);

	// std::vector<double> ifftwResult = ifftwC2r(fftwResult);
	// printVector("ifftw: ", ifftwResult);
	// fftwResult = expandFftwResult(fftwR2c(samples));// Get the result back.



	// //std::vector<std::complex<double>> iterativeIfftIcpResult = iterativeIcpFft(iterativeFftIcpResult, true);
	// // printVector("iterative icp ifft: ", iterativeIfftIcpResult);

	// t1 = std::chrono::high_resolution_clock::now();
    // std::vector<std::complex<double>> naiveResult_multiThreaded = naiveDFT_multiThreaded(samples);
	// t2 = std::chrono::high_resolution_clock::now();
	// duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	// std::cout << "dft multithreaded time: " << duration << " microseconds\n";
    // printVector("naive dft multiThreaded: ", naiveResult_multiThreaded);



	t1 = std::chrono::high_resolution_clock::now();
	std::vector<std::complex<double>> fftResult = recursiveFFT(samples, omega);
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "fft recursive fft time: " << duration << " microseconds\n";


	// t1 = std::chrono::high_resolution_clock::now();
	// std::vector<std::complex<double>> naiveResult = naiveDFT(samples);	
	// t2 = std::chrono::high_resolution_clock::now();
	// duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	// std::cout << "naive dft time: " << duration << " microseconds\n";





	autoValidate(FUNC_NAME(recursiveFFT), samples, omega, fftwResult);
    //autoValidate(FUNC_NAME(naiveDFT_multiThreaded), samples, fftwResult);
	//autoValidate(FUNC_NAME(naiveDFT), samples, fftwResult);



	std::vector<std::complex<double>> complexSamples(samples.begin(), samples.end());

	t1 = std::chrono::high_resolution_clock::now();
	std::vector<std::complex<double>> iterativeFftResult = iterativeIcpFft(complexSamples, false);
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	std::cout << "iterativefft Time: " << duration << " microseconds\n";

	for (int i = 1; i < 12; i+=1) {
		t1 = std::chrono::high_resolution_clock::now();
		std::vector<std::complex<double>> iterativeFftResult = multithreaded_iterativeIcpFft(complexSamples, false, i);
		t2 = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
		std::cout << "multithreaded_iterativeIcpFft for "<< i <<" threads Time: " << duration << " microseconds\n";
	}







    return 0;
}