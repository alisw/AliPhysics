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
set(MODULE PWGmuonaaf)
#add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWG/muon/aaf)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/PWG/muon
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    )

# Sources - alphabetical order
set(SRCS
  StageAndFilter.cxx
  FILTER_AODMUONWITHTRACKLETS.cxx
  FILTER_RAWMUON.cxx
  FILTER_ESDMUON.cxx
  FILTER_AODMUONONLY_PP2015.cxx
  FILTER_AODMUONONLY_PBPB2015.cxx
  FILTER_AODMUONONLY_PBPB2015_WITHMULTSEL.cxx
  FILTER_ESD2AODMUON.cxx
  )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

set(EXE aaf-stage-and-filter)

add_executable(${EXE} ${EXE}.cxx)
add_executable(aliraw-generate-entrylist aliraw-generate-entrylist.cxx)
add_executable(aliraw-apply-entrylist aliraw-apply-entrylist.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS PWGmuon RAWDatabase)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})
target_link_libraries(${EXE} ${MODULE} Core)
target_link_libraries(aliraw-generate-entrylist ${MODULE} Core RAWDatarec)
target_link_libraries(aliraw-apply-entrylist ${MODULE} Core RIO )

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

install(TARGETS ${EXE} RUNTIME DESTINATION bin)
install(TARGETS aliraw-generate-entrylist RUNTIME DESTINATION bin)
install(TARGETS aliraw-apply-entrylist RUNTIME DESTINATION bin)
