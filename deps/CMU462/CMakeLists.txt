cmake_minimum_required(VERSION 2.8)

set(CMU462_SOURCE
        src/vector2D.cpp
        src/vector3D.cpp
        src/vector4D.cpp
        src/matrix3x3.cpp
        src/matrix4x4.cpp
        )

add_library(CMU462 ${CMU462_SOURCE})
target_include_directories(CMU462 PUBLIC include)
