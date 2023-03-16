cmake_minimum_required(VERSION 3.16)
project(libdescriptor)

set(CMAKE_CXX_STANDARD 17)
#find_package(Doxygen OPTIONAL_COMPONENTS dot)
#if (DOXYGEN_FOUND)
#    # set input and output files
#    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
#    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile.out)
#
#    # request to configure the file
#    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
#    message("Doxygen build started")
#
#    # Note: do not put "ALL" - this builds docs together with application EVERY TIME!
#    add_custom_target(docs COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating API documentation with Doxygen"
#        VERBATIM )
#else (DOXYGEN_FOUND)
#  message("Doxygen need to be installed to generate the doxygen documentation")
#endif (DOXYGEN_FOUND)


#----module dir --------------------
# Define all source files here
include_directories(./)
include_directories(SymmetryFunctions)
file(GLOB SYMFUNFILES SymmetryFunctions/*)
include_directories(Bispectrum)
file(GLOB BISFILES Bispectrum/*)

#--------------

include_directories(autodiff)

add_library(descriptor SHARED ${BISFILES} ${SYMFUNFILES} Descriptors.cpp helper.cpp)
#add_library(descriptor SHARED ${SYMFUNFILES} Descriptors.cpp helper.cpp)

set_target_properties(descriptor PROPERTIES PUBLIC_HEADER "Descriptors.hpp")

install(TARGETS descriptor
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

# ---------------- Python module -------------------------
## Python Interp
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# PYBIND11
find_package(pybind11 REQUIRED)

pybind11_add_module(libdescriptor python_bindings.cpp)
target_link_libraries(libdescriptor PUBLIC descriptor)

#TODO add install location


# ----------------- Tests ---------------------
enable_testing()
add_test(NAME SymmetryFunctions COMMAND ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/symmetry_functions_test.py)
add_test(NAME Bispectrum COMMAND ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/bispectrum_test.py)