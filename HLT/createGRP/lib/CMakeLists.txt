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
set(MODULE AliHLTCreateGRP)
set(CREATEGRPDIR ${AliRoot_SOURCE_DIR}/HLT/createGRP/lib)

if(NOT DIM_FOUND)
    message( FATAL_ERROR "createGRP needs DIM! Please provide correct -DDIMDIR to cmake!" )
endif()

message( STATUS "DIMDIR set, building HLT createGRP library" )

# Module include folder
include_directories(${CREATEGRPDIR}
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${DIMDIR}/dim
                   )

# Sources in alphabetical order
set(SRCS
    ${CREATEGRPDIR}/AliHLTCreateGRP.cxx
   )

# Headers from sources
set(CINTHDRS
    ${CREATEGRPDIR}/AliHLTCreateGRP.h
    )
  
set(HDRS
    ${CINTHDRS}
  )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "" "${CINTHDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
generate_rootmap("${MODULE}" "" "")

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
target_link_libraries(${MODULE} ${DIM_LIBRARIES})

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)
install(FILES ../AliHLTCreateGRP.C DESTINATION HLT/createGRP)