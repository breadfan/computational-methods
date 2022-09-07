cmake_minimum_required(VERSION 3.14)
project(stiff_systems)

set(CMAKE_CXX_STANDARD 14)

add_executable(stiff_systems euler_method.cpp)