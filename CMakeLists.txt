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
include_directories(SymmetryFunctions)
file(GLOB SYMFUNFILES SymmetryFunctions/*)
include_directories(Bispectrum)
file(GLOB BISFILES Bispectrum/*)

#include_directories(CM)
#file(GLOB CMFILES CM/*)
#--------------
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fuse-ld=lld -flto")

#pybind11_add_module(descriptors  descriptors.cpp descriptors.hpp python_bindings.cpp)
add_library(descriptor SHARED ${SYMFUNFILES} Descriptors.cpp helper.cpp)
#add_library(descriptor_static STATIC descriptors.cpp descriptors.hpp)

target_link_options(descriptor PRIVATE "LINKER:--lto-legacy-pass-manager")
target_link_options(descriptor PRIVATE "LINKER:-mllvm=-load=/opt/enzyme/enzyme/build/Enzyme/LLDEnzyme-12.so")
