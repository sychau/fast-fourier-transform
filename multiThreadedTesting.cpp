#include "utillity.h"
#include "fft_multiThreaded.h"

#include <chrono>
#include <thread>
#include <iomanip>


int main(int argc, char* argv[]){
    std::vector<double> samples = {1, 2, 3, 4, 5, 6, 7, 8};
    
    std::cout << std::setprecision(3);
    printVector("samples: ", samples);

    std::vector<std::complex<double>> naiveResult_multiThreaded = naiveDFT_multiThreaded(samples);
    printVector("naive dft multiThreaded: ", naiveResult_multiThreaded);

    return 0;
}