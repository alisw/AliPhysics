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
set(MODULE qpythia)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PYTHIA6/QPYTHIA)

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/FASTSIM
                    ${AliRoot_SOURCE_DIR}/PYTHIA6/pythia6
                   )

# Sources in alphabetical order
set(SRCS
    AliQPythiaWrapper.cxx
   )

# Fortran sources
set(FSRCS
    pythia6_common_block_address.F
    pythia6_called_from_cc.F
    pythia-6.4.14.f
    pyquen1_5.F
    q-pyshow.1.2.F
    ${AliRoot_SOURCE_DIR}/PYTHIA6/pythia6/pyevwt.f
   )

# C sources
set(CSRCS
    main.c
    pythia6_common_address.c
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "qpythiaLinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS FASTSIM microcern)
# TO DO: temporary fix: removing the rootmaps to enable the manual loading. Fix the dependecies!
#generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")


# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx ${FSRCS} ${CSRCS})

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
install(FILES qgrid DESTINATION PYTHIA6/QPYTHIA)
