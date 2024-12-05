# Specify the minimum CMake version
cmake_minimum_required(VERSION 3.20)

# Define the project name and version
project(cmpt431_fft VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Use mpic++ wrapper
set(MPI_COMPILER mpic++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Find MPI (ensures that MPI is configured correctly)
find_package(MPI REQUIRED)

# Add source files
add_executable(test_serial test_serial.cpp fft_serial.cpp utillity.cpp)
add_executable(test_parallel test_parallel.cpp fft_multiThreaded.cpp fft_serial.cpp utillity.cpp)
add_executable(test_distributed test_distributed.cpp fft_distributed.cpp fft_serial.cpp utillity.cpp)

# Link the FFTW library
target_link_libraries(test_serial fftw3)
target_link_libraries(test_parallel fftw3)
target_link_libraries(test_distributed fftw3)
target_link_libraries(test_distributed MPI::MPI_CXX)
