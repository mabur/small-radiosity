# Top-Level CMakeList.txt
cmake_minimum_required(VERSION 3.4)
enable_language(CXX)
add_definitions(-std=c++11)

PROJECT(small_radiosity)
ADD_EXECUTABLE(small_radiosity main.cpp)

