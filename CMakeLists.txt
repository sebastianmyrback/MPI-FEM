cmake_minimum_required(VERSION 3.16)
project(TensorFEM)

find_package(Python3 COMPONENTS Interpreter Development NumPy)
include_directories(${Python3_INCLUDE_DIRS} include)

set(CMAKE_CXX_STANDARD 17)

# Find all the source files in the directory
file(GLOB SOURCES "src/*.cpp")

# Add source files to library
add_library(TensorFEM SHARED ${SOURCES})

add_executable(main examples/main.cpp)

target_link_libraries(main TensorFEM ${Python3_LIBRARIES})

