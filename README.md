# Fast Fourier Transform (FFT) Implementation  

This project provides serial, parallel, and distributed implementations of the Fast Fourier Transform (FFT) algorithm.

## Dependencies  

The following dependencies are required to build and run the project:  
- **openMPI**  
- **cmake ver3.20**  
- **make**  
- **g++10**  
- **fftw3** (used to verify the results of our algorithm)  

### Installing fftw3 on Ubuntu  
If you're using Ubuntu, you can install fftw3 with the following command:  
```bash  
sudo apt install libfftw3-dev  
```
## Build instructiuon
1. Create a build directory and navigate to it:
```bash
mkdir build  
cd build
```
2. Run `cmake` to generate the build files:
```bash
cmake ..
```

3. Build the project using the generated files:
```bash
cmake --build .  
```

## Usage
The size of the input array is $2^{sizeExp}$, and a seed value is used to generate random samples.

### Run Serial Test
Test the serial implementation of the FFT algorithm
```bash
./test_serial --sizeExp 20 --seed 49
```

### Run Parallel Test
Test the parallel implementation with multithreading. Specify the number of threads with `--nThreads`
```bash
./test_parallel --nThreads 4 --sizeExp 20 --seed 49
```

### Run Distributed Test
Test the distributed implementation using MPI. Specify the number of processes with `-n`
```bash
mpirun -n 8 ./test_distributed --sizeExp 20 --seed 49
```


## Plotting


### Requirements for creating plots
Matplotlib must be installed to be able to create plots. This can be done with: 
```bash
pip3 install matplotlib
```

### Create Parallel Plots
To begin we can first navigate to the plotting directory from the build directory like so: 
```bash
cd ../plotting
```
From there we can gather all of the necessary data as follows: 
```bash
python3 get_data_parallel.py
```
Afterwhich we can generate plots by running the plot constructing scripts and selecting parallel plots as follows: 
```bash
python3 construct_plots.py
>>>Enter 1 for parallel plots, 2 for distributed plots: 1
```
After which there should be the parallel plots in the actual_plots folder. 

### Create Distributed Plots
To begin we can first navigate to the plotting directory from the build directory like so: 
```bash
cd ../plotting
```
From there we can gather all of the necessary data as follows: 
```bash
python3 get_data_distributed.py
```
Afterwhich we can generate plots by running the plot constructing scripts and selecting parallel plots as follows: 

```bash
python3 construct_plots.py
>>>Enter 1 for parallel plots, 2 for distributed plots: 2
```
After which there should be the distributed plots in the actual_plots folder. 


## Note
### Include Path Configuration in VSCode
If you see warnings in VSCode, add the following paths to your c_cpp_properties.json:
```json
"/usr/include",  
"/usr/lib/x86_64-linux-gnu/openmpi/**"  
```
### Numerical Precision
For very large input arrays (e.g. $2^{21}$ or larger), results may differ slightly from fftw due to accumulated floating-point errors

## References  
- **Introduction to Parallel Computing, Ch. 13 (IPC)**  
  [Link to O'Reilly](https://learning.oreilly.com/library/view/introduction-to-parallel/0201648652/)  

- **Cooley–Tukey FFT Algorithm**:  
  [Wikipedia](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm)  

- **Construction of a High-Performance FFT**:  
  [EDP Paper](https://edp.org/work/Construction.pdf)  