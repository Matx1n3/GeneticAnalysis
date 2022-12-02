cmake_minimum_required(VERSION 3.24)
project(GeneticAnalysis C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(GeneticAnalysis
        definegg.h
        fungg.h
        fungg_s.c
        gengroups_s.c)
