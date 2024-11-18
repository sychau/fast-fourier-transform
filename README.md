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