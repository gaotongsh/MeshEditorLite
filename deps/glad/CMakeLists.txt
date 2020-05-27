cmake_minimum_required(VERSION 3.15)
project(glad)
set(CMAKE_CXX_STANDARD 14)

add_library(glad src/glad.c)
target_include_directories(glad PUBLIC include)
