cmake_minimum_required(VERSION 3.14)
project(hw2_1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(OpenMP REQUIRED)

add_executable(serial main.cpp serial.cpp common.h)

add_executable(openmp main.cpp openmp.cpp common.h)
target_link_libraries(openmp PRIVATE OpenMP::OpenMP_CXX)
