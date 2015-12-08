cmake_minimum_required(VERSION 3.3)
project(CSN6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(CSN6 ${SOURCE_FILES})