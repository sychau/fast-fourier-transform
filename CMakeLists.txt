# Specify the minimum CMake version
cmake_minimum_required(VERSION 3.20)

# Define the project name and version
project(cmpt431_fft VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Use mpic++ wrapper
set(MPI_COMPILER mpic++)
set(CMAKE_CXX_FLAGS "-O3")

# Find MPI (ensures that MPI is configured correctly)
find_package(MPI REQUIRED)

# Add source files
add_executable(fft_test testing.cpp fft_serial.cpp fft_multiThreaded.cpp utillity.cpp)
add_executable(fft_test_distributed testing_distributed.cpp fft_distributed.cpp fft_serial.cpp utillity.cpp)

# Link the FFTW library
target_link_libraries(fft_test fftw3)
target_link_libraries(fft_test_distributed fftw3)
target_link_libraries(fft_test_distributed MPI::MPI_CXX)
