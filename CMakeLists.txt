cmake_minimum_required(VERSION 3.16)
project(libdescriptor)

set(CMAKE_CXX_STANDARD 17)
find_package(Doxygen OPTIONAL_COMPONENTS dot)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile.out)

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
# Symmetry functions calc test
add_test(NAME SymmetryFunctionsCalc COMMAND ${PYTHON_EXECUTABLE}  symmetry_function_calc_test.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/autodiff_tests/symmetry_function_test)
set_tests_properties(SymmetryFunctionsCalc PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$ENV{LD_LIBRARY_PATH};PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH}")

# Symmetry functions load test
add_test(NAME SymmetryFunctionsLoad COMMAND ${PYTHON_EXECUTABLE}  symmetry_function_load_test.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/autodiff_tests/symmetry_function_test)
set_tests_properties(SymmetryFunctionsCalc PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$ENV{LD_LIBRARY_PATH};PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH}")

# Bispectrum calc test
add_test(NAME BispectrumCalc COMMAND ${PYTHON_EXECUTABLE}  bispectrum_calc_test.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/autodiff_tests/bispectrum_test)
set_tests_properties(BispectrumCalc PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$ENV{LD_LIBRARY_PATH};PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH}")

# Descriptor load test
add_test(NAME BispectrumLoad COMMAND ${PYTHON_EXECUTABLE}  bispectrum_load_test.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/autodiff_tests/bispectrum_test)
set_tests_properties(BispectrumLoad PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$ENV{LD_LIBRARY_PATH};PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH}")
