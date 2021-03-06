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
## Configuration file for the opt subfolder (external tools)
##
## author Come Raczy
##
################################################################################

message(STATUS "Creating external tools in subdirectory opt")

include ("${MANTA_GLOBALS_CMAKE}")

# convenience macro to set in both current and parent scope:
macro(superset symbol value)
    set(${symbol} "${value}")
    set(${symbol} "${value}" PARENT_SCOPE)
endmacro()

#
# samtools
#
set(SAMTOOLS_PREFIX "samtools-0.1.18_no_tview")
superset(SAMTOOLS_DIR "${CMAKE_CURRENT_BINARY_DIR}/${SAMTOOLS_PREFIX}")
superset(SAMTOOLS_LIBRARY "${SAMTOOLS_DIR}/libbam.a")

add_custom_command(OUTPUT ${SAMTOOLS_LIBRARY}
                   COMMAND rm -rf "${SAMTOOLS_DIR}"
                   COMMAND tar -xjf "${MANTA_REDIST_DIR}/${SAMTOOLS_PREFIX}.tar.bz2"
                   COMMAND $(MAKE) -C "${SAMTOOLS_DIR}" all bgzip
                   COMMENT "Building the SAM/BAM library and tools")
add_custom_target(MANTA_SAMTOOLS DEPENDS "${SAMTOOLS_LIBRARY}")

install(PROGRAMS "${SAMTOOLS_DIR}/samtools" DESTINATION "${MANTA_LIBEXECDIR}")
install(PROGRAMS "${SAMTOOLS_DIR}/bgzip"    DESTINATION "${MANTA_LIBEXECDIR}")

#
# tabix
#
set(TABIX_PREFIX "tabix-0.2.6")
superset(TABIX_DIR "${CMAKE_CURRENT_BINARY_DIR}/${TABIX_PREFIX}")
set(TABIX_PROG "${TABIX_DIR}/tabix")
superset(TABIX_LIBRARY "${TABIX_DIR}/libtabix.a")

add_custom_command(OUTPUT ${TABIX_PROG}
                   COMMAND rm -rf "${TABIX_DIR}"
                   COMMAND tar -xjf "${MANTA_REDIST_DIR}/${TABIX_PREFIX}.tar.bz2"
                   COMMAND $(MAKE) -C "${TABIX_DIR}" all
                   COMMENT "Building tabix")
add_custom_target(MANTA_TABIX DEPENDS "${TABIX_PROG}")

install(PROGRAMS "${TABIX_PROG}" DESTINATION "${MANTA_LIBEXECDIR}")


#
# pyflow
#

#make sure unpacking happens regardless of whether ${VCF_TOOLS_DIR} exists
set(PYFLOW_PREFIX "pyflow-1.0.0")
set(PYFLOW_DIR "${CMAKE_CURRENT_BINARY_DIR}/${PYFLOW_PREFIX}")
set(PYFLOW_SCRIPT "${PYFLOW_DIR}/src/pyflow.py")

add_custom_command(OUTPUT ${PYFLOW_SCRIPT}
                   COMMAND rm -rf "${PYFLOW_DIR}"
                   COMMAND tar -xjf "${MANTA_REDIST_DIR}/${PYFLOW_PREFIX}.tar.bz2"
                   COMMENT "Extracting pyflow")
add_custom_target(MANTA_PYFLOW DEPENDS "${PYFLOW_SCRIPT}")

install (DIRECTORY ${PYFLOW_DIR}/src/ DESTINATION "${MANTA_PYTHON_LIBDIR}/pyflow" FILES_MATCHING PATTERN "*.py")

# tie results back to parent:
#
add_dependencies(MANTA_OPT MANTA_SAMTOOLS MANTA_TABIX MANTA_PYFLOW)


