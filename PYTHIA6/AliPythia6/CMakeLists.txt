# **************************************************************************
# * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
# *                                                                        *
# * Author: The ALICE Off-line Project.                                    *
# * Contributors are mentioned in the code where appropriate.              *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************

# Module
set(MODULE AliPythia6)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PYTHIA6/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/EVGEN
                    ${AliRoot_SOURCE_DIR}/FASTSIM
                    ${AliRoot_SOURCE_DIR}/LHAPDF/base
                    ${AliRoot_SOURCE_DIR}/PYTHIA6/pythia6
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    AliPythiaRndm.cxx
    AliDecayerPythia.cxx
    AliGenPythia.cxx
    AliGenPythiaPlus.cxx
    AliPythia6.cxx
    AliPythiaBase.cxx
    AliPythia.cxx
   )

# Fortran sources
set(FSRCS
    openDecayTable.F
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Additional header files
set(HDRS ${HDRS}
    PythiaProcesses.h
    PyquenCommon.h
   )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER EVGEN FASTSIM EGPythia6 microcern lhapdfbase)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")


# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx ${FSRCS})

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)
