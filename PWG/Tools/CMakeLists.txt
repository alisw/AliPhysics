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
set(MODULE PWGTools)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWG/Tools)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools/yaml-cpp/include
  )

# Sources - alphabetical order
set(SRCS
  AliAnalysisHelperJetTasks.cxx
  AliBasicParticle.cxx
  AliTHn.cxx
  AliPWGHistoTools.cxx
  AliPWGFunc.cxx
  AliLatexTable.cxx
  AliFigure.cxx
  AliCanvas.cxx
  AliHelperPID.cxx
  AliNamedArrayI.cxx
  AliNamedString.cxx
  TCustomBinning.cxx
  TLinearBinning.cxx
  TVariableBinning.cxx
  THistManager.cxx
  AliYAMLConfiguration.cxx
  AliJSONReader.cxx
  AliJSONData.cxx
  AliAnalysisTaskDummy.cxx
  AliTLorentzVector.cxx
  )

if(${ROOT_VERSION} GREATER_EQUAL 6.0)
  list(APPEND SRCS AliMCSpectraWeights.cxx)
endif()

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
set(HDRS
  "${HDRS}"
  TBinning.h
  )

# Statically build YAML
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Disable additional yaml tools (which are not used)")
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Disable yaml-cpp tests (which we are not running)")
add_subdirectory(yaml-cpp)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSIS AOD CORRFW ESD STEERBase ASImage OADB)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Link against yaml-cpp. It must be included _after_ the ROOT map because it is static rather than shared!
set(LIBDEPS ${LIBDEPS} yaml-cpp)

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

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

# Installing the macros
install (DIRECTORY macros DESTINATION PWG/Tools)

# Tests
install (DIRECTORY test DESTINATION PWG/Tools)

# Histmanager test
set(HISTMGRTESTS
    build_simple
    build_grouped
    fill_simple
    fill_grouped
    )
foreach(TEST_HMGR ${HISTMGRTESTS})
    add_test (histmgr_${TEST_HMGR}
        env
        LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}
        DYLD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{DYLD_LIBRARY_PATH}
        root -l -b -q "${CMAKE_INSTALL_PREFIX}/PWG/tools/test/histmgr/runtest.C(\"${TEST_HMGR}\")")
endforeach()
