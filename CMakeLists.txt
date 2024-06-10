cmake_minimum_required(VERSION 3.22.1)
project(openmp_examples)

set(CMAKE_CXX_STANDARD 17)

# for macos required https://gist.github.com/ijleesw/4f863543a50294e3ba54acf588a4a421
#set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
#set(CMAKE_OSX_ARCHITECTURES "arm64")
# use for x86_64
# set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(openmp_examples main.cpp)
