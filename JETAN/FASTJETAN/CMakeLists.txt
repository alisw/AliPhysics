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
set(MODULE FASTJETAN)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/JETAN/${MODULE})

# Additional include folders in alphabetical order
include_directories(${AliPhysics_SOURCE_DIR}/JETAN/JETAN
                   )

# Sources in alphabetical order
set(SRCS
    AliAnalysisTaskJetBackgroundSubtract.cxx
    AliAnalysisTaskJetCluster.cxx
    AliAnalysisTaskJetClusterKine.cxx
    AliFastJetBkg.cxx
    AliFastJetFinder.cxx
    AliFastJetHeaderV1.cxx
    AliFastJetInput.cxx
    AliSISConeJetFinder.cxx
    AliSISConeJetHeader.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

if(FASTJET_FOUND)
    include_directories(SYSTEM ${FASTJET_INCLUDE_DIR})
    link_directories(${FASTJET_LIBS_DIR})
    add_definitions(${FASTJET_DEFINITIONS})
endif(FASTJET_FOUND)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}" "${FASTJET_ROOTDICT_OPTS}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase ESD AOD ANALYSIS ANALYSISalice JETAN EG Net)
if(FASTJET_FOUND)
  set(LIBDEPS
    ${LIBDEPS}
    fastjetcontribfragile
    fastjetplugins
    siscone_spherical
    siscone
    fastjettools
    fastjet
    )
endif(FASTJET_FOUND)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Create an object to be reused in case of static libraries
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)
# Add a library to the project using the object
add_library_tested(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)
target_link_libraries(${MODULE} ${LIBDEPS})


# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)
