# Specify the minimum CMake version
cmake_minimum_required(VERSION 3.20)

# Define the project name and version
project(cmpt431_fft VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set to how close two complex numbers must be to be considered equal
set(VALIDATE_EPSILON 1e-12)
configure_file(utillity.cpp.in ../utillity.cpp)

# Add source files
add_executable(fft_test testing.cpp fft_serial.cpp fft_multiThreaded.cpp utillity.cpp)

# Link the FFTW library
target_link_libraries(fft_test fftw3)