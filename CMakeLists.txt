cmake_minimum_required(VERSION 3.15)
project(MeshEditorLite)
set(CMAKE_CXX_STANDARD 14)

include_directories(include)

# OpenGL
find_package(OpenGL REQUIRED)

# glad
add_subdirectory(deps/glad)

# glfw
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(deps/glfw)

# glm
add_subdirectory(deps/glm)

# CMU462
add_subdirectory(deps/CMU462)

# executable
add_executable(mesh_editor
        include/camera.h
        include/shader.h
        include/obj_loader_halfedge.h
        include/halfEdgeMesh.h
        include/meshEdit.h
        include/mutablePriorityQueue.h
        src/main.cpp
        src/halfEdgeMesh.cpp
        src/meshEdit.cpp
        )
target_link_libraries(mesh_editor OpenGL::GL glad glfw glm CMU462)
