cmake_minimum_required(VERSION 3.22)
project(TraceMeSenpai)

set(CMAKE_CXX_STANDARD 14)

add_executable(TraceMeSenpai 1705077_main.cpp)
#target_link_libraries(TraceMeSenpai libopengl32.A libglu32.A libfreeglut.A)
target_link_libraries(TraceMeSenpai -lOpenGL32 -lglu32  -lfreeGLUT)
