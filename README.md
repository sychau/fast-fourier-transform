Build
```
mkdir build
cd build
cmake ..
cmake --build .
./fft
```


Install fftw3
```
sudo apt install libfftw3-dev
```

If vscode shows warnings, add "/usr/include" to includePath in c_cpp_properties.json

Resource:
IPC Ch.13
https://learning.oreilly.com/library/view/introduction-to-parallel/0201648652/

Cooley–Tukey FFT algorithm
https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm