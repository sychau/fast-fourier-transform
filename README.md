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
You can choose any number of processes that is power of 2
```
mpirun -n 4 ./fft_test_distributed
```

Install fftw3
```
sudo apt install libfftw3-dev
```

If vscode shows warnings, add "/usr/include" to includePath in c_cpp_properties.json
Add openMPI to includePath "/usr/lib/x86_64-linux-gnu/openmpi/**"

As array size increase, error will increase, at around 2^20 it will exceed the EPSILON 1e-4

Resource:
IPC Ch.13
https://learning.oreilly.com/library/view/introduction-to-parallel/0201648652/

Cooley–Tukey FFT algorithm
https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm