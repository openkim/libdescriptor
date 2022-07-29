cmake_minimum_required(VERSION 3.16)
project(libdescriptor)

set(CMAKE_CXX_STANDARD 14)
find_package(EnvModules REQUIRED)
env_module(load python39)

#----Pybind11
find_package(PythonInterp REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
find_package(pybind11 REQUIRED)
#------------

#----module dir
include_directories(./)
include_directories(SymFun)
file(GLOB SYMFUNFILES SymFun/*)
#include_directories(CM)
#file(GLOB CMFILES CM/*)
#--------------

pybind11_add_module(descriptors  descriptors.cpp descriptors.hpp python_bindings.cpp)
add_library(descriptor SHARED ${SYMFUNFILES} descriptors.cpp descriptors.hpp)
add_library(descriptor_static STATIC descriptors.cpp descriptors.hpp)

target_compile_options(descriptors PRIVATE "SHELL: -Xclang -load -Xclang /opt/enzyme/enzyme/build/Enzyme/ClangEnzyme-12.so")
target_compile_options(descriptor PRIVATE "SHELL: -Xclang -load -Xclang /opt/enzyme/enzyme/build/Enzyme/ClangEnzyme-12.so")
target_compile_options(descriptor_static PRIVATE "SHELL: -Xclang -load -Xclang /opt/enzyme/enzyme/build/Enzyme/ClangEnzyme-12.so")
