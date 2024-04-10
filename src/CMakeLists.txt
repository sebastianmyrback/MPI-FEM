cmake_minimum_required(VERSION 3.16)
project(TensorFEM)
set(CMAKE_CXX_STANDARD 17)

# Step 1: Add the directory that contains the header files
include_directories(${CMAKE_SOURCE_DIR}/mesh)

# Step 2: Find all the source files in the directory
file(GLOB SOURCES "mesh/*.cpp")

# Step 3: Add these source files to your executable
add_executable(main main.cpp ${SOURCES})