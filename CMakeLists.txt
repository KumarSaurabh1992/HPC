cmake_minimum_required(VERSION 3.8)
project(tutorials)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")

set (INCLUDES
        include/DGEMM.h)
set(SRC
        src/DGEMM.cpp)


add_executable(tutorials  ${INCLUDES} ${SRC} main.cpp)