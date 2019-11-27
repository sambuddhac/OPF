# Specify the minimum version for CMake

cmake_minimum_required(VERSION 2.8)

# Set the output folder where your program will be created
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-std=c++0x")
#set(CMAKE_BINARY_DIR $/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/glpk-4.65/examples/UC/cmake-3.11.1/bin)
#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
# Project's name
project(SerialCode)
# The following folder will be included
#include_directories("$/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/glpk-4.65/examples/UC")
add_executable(
SerialCode
node.cpp 
generator.cpp 
load.cpp 
transl.cpp 
network.cpp 
maintwoserialOPF.cpp 
gensolver.cpp
loadsolver.cpp
transolver.cpp
)
