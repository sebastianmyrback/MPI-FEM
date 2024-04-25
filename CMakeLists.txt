# Create a Cmake list that build the library. The header files are in the include directory and the source files are in the src directory
cmake_minimum_required(VERSION 3.16)
project(TensorFEM)
set(CMAKE_CXX_STANDARD 17)

# Step 1: Add the directory that contains the header files
include_directories(${CMAKE_SOURCE_DIR}/include)

# Step 2: Find all the source files in the directory
file(GLOB SOURCES "src/*.cpp")

# Step 3: Add these source files to your library
add_library(TensorFEM SHARED ${SOURCES})

set(MAIN_DIR ${CMAKE_SOURCE_DIR}/examples)
# Step 4: Add the executable that uses the library
add_executable(main ${MAIN_DIR}/main.cpp)
target_link_libraries(main TensorFEM)

