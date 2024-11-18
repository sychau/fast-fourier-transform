# Specify the minimum CMake version
cmake_minimum_required(VERSION 3.20)

# Define the project name and version
project(cmpt431_fft VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add source files
add_executable(fft fft_serial.cpp)

# Link the FFTW library
target_link_libraries(fft fftw3)