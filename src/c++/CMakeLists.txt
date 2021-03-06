#
# Manta
# Copyright (c) 2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

################################################################################
##
## Configuration file for the c++ subfolder
##
## author Come Raczy
##
################################################################################

set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set (MANTA_CXX_EXECUTABLE_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxExecutable.cmake")
set (MANTA_CXX_LIBRARY_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxLibrary.cmake")
set (MANTA_CXX_CONFIGURE_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxConfigure.cmake")
set (MANTA_CXX_TEST_LIBRARY_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxTestLibrary.cmake")

include ("${MANTA_CXX_CONFIGURE_CMAKE}")

add_custom_target(MANTA_CPPCHECK
    ALL
    COMMAND ${MANTA_SOURCE_QC_DIR}/run_cppcheck.bash ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_dependencies(MANTA_OPT MANTA_CPPCHECK)

##
## Build all the libraries for the project
##
add_subdirectory (lib)

##
## build all the applications for the project
##
add_subdirectory (bin)


##
## build the documentation when available
##
include  (FindDoxygen)
message (STATUS "Doxygen: ${DOXYGEN_EXECUTABLE}. Dot: ${DOXYGEN_DOT_EXECUTABLE}.")
if (DOXYGEN_FOUND)
    set (DOXYFILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    message (STATUS "Creating Doxygen config file: ${DOXYFILE}")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${DOXYFILE} @ONLY)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOXYFILE})
endif ()
