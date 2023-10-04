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

#----------------
# Find enzyme at location ENZYME_LIB
if(NOT DEFINED ENZYME_LIB)
    message(WARNING "ENZYME_LIB variable not defined, switching to default: /usr/local/lib" )
    set(ENZYME_LIB "/usr/local/lib")
endif()
message("Searching for ENZYME at ${ENZYME_LIB}")
if(EXISTS "${ENZYME_LIB}/LLDEnzyme-13.so")
    message("Found LLDEnzyme plugin.")
else()
    message(FATAL_ERROR "LLDEnzyme plugin not found, please check your Enzyme installation.")
endif()
if(EXISTS "${ENZYME_LIB}/ClangEnzyme-13.so")
    message("Found ClangEnzyme plugin.")
else()
    message(WARNING "ClangEnzyme plugin not found, it is NOT needed for libdescriptor, but usually indicates improper installation of Enzyme.")
endif()
if(EXISTS "${ENZYME_LIB}/LLVMEnzyme-13.so")
    message("Found LLVMEnzyme plugin.")
else()
    message(WARNING "LLVMEnzyme plugin not found, it is NOT needed for libdescriptor, but usually indicates improper installation of Enzyme.")
endif()


#----module dir --------------------
# Define all source files here
include_directories(./)
include_directories(SymmetryFunctions)
file(GLOB SYMFUNFILES SymmetryFunctions/*)
include_directories(Bispectrum)
file(GLOB BISFILES Bispectrum/*)
include_directories(SOAP)
file(GLOB SOAPFILES SOAP/*)
#--------------

# Set LLD, just in case... sometimes default LLD the name of LLD file is lld-12
if (NOT DEFINED LLD)
    set(LLD "lld")
endif()

# Compile files, but use lld for linking. lld is patched with enzyme to generate AD code
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fuse-ld=${LLD} -flto -Xclang -new-struct-path-tbaa")

add_library(descriptor SHARED ${BISFILES} ${SYMFUNFILES} ${SOAPFILES} Descriptors.cpp helper.cpp)
#add_library(descriptor_static STATIC ${BISFILES} ${SYMFUNFILES} ${SOAPFILES} Descriptors.cpp helper.cpp)


# pass linking to lld, with enzyme plugin
target_link_options(descriptor PRIVATE "LINKER:--lto-legacy-pass-manager")
target_link_options(descriptor PRIVATE "LINKER:-mllvm=-load=${ENZYME_LIB}/LLDEnzyme-13.so")

# TODO REMOVE EXTRA HEADER FILES ONCE ENZYME PROVIDES SOLUTION
#set_target_properties(descriptor PROPERTIES PUBLIC_HEADER "Descriptors.hpp;SymmetryFunctions/SymmetryFunctions.hpp;Bispectrum/Bispectrum.hpp;helper.hpp")
# Further testing indicates that library might be working fine with virtual desctructors. So everting the
# changes for now. But leaving above as reminder
set_target_properties(descriptor PROPERTIES PUBLIC_HEADER "Descriptors.hpp")
#set_target_properties(descriptor_static PROPERTIES OUTPUT_NAME descriptor)

install(TARGETS descriptor
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

# ---------------- Python module -------------------------
# Python Interp
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# PYBIND11
find_package(pybind11 REQUIRED)

pybind11_add_module(libdescriptor python_bindings.cpp)
target_link_libraries(libdescriptor PUBLIC descriptor)

#TODO add install location