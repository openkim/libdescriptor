cmake_minimum_required(VERSION 3.16)
project(libdescriptor)

set(CMAKE_CXX_STANDARD 14)
find_package(Doxygen OPTIONAL_COMPONENTS dot)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.out)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # Note: do not put "ALL" - this builds docs together with application EVERY TIME!
    add_custom_target(docs COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

#find_package(EnvModules REQUIRED)
#env_module(load python39)
#
##----Pybind11
#find_package(PythonInterp REQUIRED)
#include_directories(${PYTHON_INCLUDE_DIRS})
#find_package(pybind11 REQUIRED)
##------------

#----module dir
include_directories(./)
include_directories(SymmetryFunctions)
file(GLOB SYMFUNFILES SymmetryFunctions/*)
include_directories(Bispectrum)
file(GLOB BISFILES Bispectrum/*)

#--------------
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fuse-ld=lld -flto -Xclang -new-struct-path-tbaa")

#pybind11_add_module(descriptors  descriptors.cpp descriptors.hpp python_bindings.cpp)
add_library(descriptor SHARED ${BISFILES} ${SYMFUNFILES} Descriptors.cpp helper.cpp)
#add_library(descriptor_static STATIC descriptors.cpp descriptors.hpp)

target_link_options(descriptor PRIVATE "LINKER:--lto-legacy-pass-manager")
target_link_options(descriptor PRIVATE "LINKER:-mllvm=-load=/opt/enzyme/enzyme/build/Enzyme/LLDEnzyme-12.so")
