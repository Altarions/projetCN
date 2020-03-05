cmake_minimum_required(VERSION 3.15)
project(projetC_NA)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(projetC_NA
        main.cpp
        Makefile
        matrix.cpp
        matrix.h
        )
