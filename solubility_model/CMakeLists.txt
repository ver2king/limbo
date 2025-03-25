cmake_minimum_required(VERSION 3.22.1)

#TODO: ? Requirement to set the compilers here
#SET(CMAKE_C_COMPILER C:/msys64/mingw64/bin/gcc.exe)
#SET(CMAKE_CXX_COMPILER C:/msys64/mingw64/bin/g++.exe)
#SET(CMAKE_MAKE_PROGRAM C:/msys64/mingw64/bin/cmake.exe)

project(SOL_OPT VERSION 1.0
        DESCRIPTION "SOLUBILITY MODEL OPTIMIZATION"
        LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)

find_package (Eigen3 3.3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

include(GoogleTest)
find_package(GTest REQUIRED)

add_executable(Test test.cpp)

add_library(Utils utils/utils.hpp utils/utils.cpp)
add_library(Model src/model.hpp src/model.cpp)
add_library(OptFunc opt_func/opt_func.hpp opt_func/opt_func.cpp)
add_library(Data data/data.hpp data/data.cpp)

target_link_libraries(Utils PUBLIC Eigen3::Eigen)
target_link_libraries(Data PUBLIC Utils)

target_link_libraries(Model PUBLIC Utils)
target_link_libraries(Model PUBLIC Data)

target_link_libraries(OptFunc PUBLIC Utils)
target_link_libraries(OptFunc PUBLIC Data)
target_link_libraries(OptFunc PUBLIC Model)

target_link_libraries(Test PUBLIC Utils)
target_link_libraries(Test PUBLIC Data)
target_link_libraries(Test PUBLIC Model)
target_link_libraries(Test PUBLIC OptFunc)
target_link_libraries(Test PUBLIC ${Boost_LIBRARIES} )