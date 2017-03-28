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
# **************************************************************************/

#Module
set(MODULE PWGCFCorrelationsThreePart)
add_definitions(-D_MODULE_="${MODULE}")


if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang OR CMAKE_COMPILER_IS_GNUCXX)
  if (NOT ${CMAKE_SYSTEM_VERSION} MATCHES el5)
    message(STATUS "Treating warnings as errors under ${CMAKE_CURRENT_SOURCE_DIR}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
  endif()
endif()

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGCF/Correlations/ThreePart)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
  )

# Sources - alphabetical order
set(SRCS
  AliAnalysisTaskBuildCorrTree.cxx
  AliAnalysisTaskCorrelation3p.cxx
  AliAnalysisTaskCorrelation3p_lightefficiency.cxx
  AliCFPI0.cxx
  AliCorrelation3p.cxx
  AliCorrelation3p_noQA.cxx
  AliFilteredEvent.cxx
  AliFilteredEventInputHandler.cxx
  AliFilteredTrack.cxx
  AliThreeParticleCorrelator.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSIS ANALYSISalice STEERBase ESD PWGLFforward2 AOD)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library. Note: the extra path has a variable that is NOT
# substituted by CMake (it has rounded and not curly braces). The $(VARIABLE) in this format will be
# properly substituted by make with the correct path in the remote location. So, rounded braces
# here are not an error!
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}" "$(ALICE_ROOT)/include/pythia")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)

# Install the macros
install(DIRECTORY macros DESTINATION PWGCF/Correlations/ThreePart)