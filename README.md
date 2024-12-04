Build
```
mkdir build
cd build
cmake ..
cmake --build .
```

Run single process test
```
./fft_test
```

Run distributed test
Input number of processes, array size exponent(power of 2) and random seed
```
mpirun -n 8 ./fft_test_distributed --sizeExp 16 --seed 49
```

Install fftw3
```
sudo apt install libfftw3-dev
```

If vscode shows warnings, add "/usr/include" to includePath in c_cpp_properties.json
Add openMPI to includePath "/usr/lib/x86_64-linux-gnu/openmpi/**"

If you run very large array (> 2^21), the result will start to differ from fftw as error add up

Resource:
IPC Ch.13
https://learning.oreilly.com/library/view/introduction-to-parallel/0201648652/

Cooley–Tukey FFT algorithm
https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm

Construction of a High-Performance FFT 
https://edp.org/work/Construction.pdf
