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
set(MODULE PWGHMTF)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWG/HMTF)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
  )

# Sources - alphabetical order
set(SRCS
  AliAnalysisTaskHMTFMCMultEst.cxx
  AliAnalysisTrackingUncertaintiesHMTF.cxx
  AliEventClassifierBase.cxx
  AliEventClassifierMult.cxx
  AliEventClassifierMPI.cxx
  AliEventClassifierSphericity.cxx
  AliEventClassifierSpherocity.cxx
  AliEventClassifierQ2.cxx
  AliIsPi0PhysicalPrimary.cxx
  AliObservableBase.cxx
  AliObservableClassifierpTPID.cxx
  AliObservableCorrelationsOfClassifiers.cxx
  AliObservableEtaNch.cxx
  MonteCarlo/AliAnalysisTaskHMTFMC.cxx
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
set(LIBDEPS ANALYSIS ANALYSISalice)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

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

## Installing the macros
install (FILES
  macros/AddMCGenDPMJET.C
  macros/AddTaskHMTFMCMultEst.C
  macros/AddTaskTrackingUncertHMTF.C
  DESTINATION PWG/HMTF/macros)

set (INPUT_FILES
  MonteCarlo/template.jdl
  MonteCarlo/templateEPOS3111.jdl
  MonteCarlo/DIPSYpp_HepMC.in
  MonteCarlo/crmc_template.param
  MonteCarlo/pp.optns
  MonteCarlo/pPb.optns
  MonteCarlo/PbPb.optns
  MonteCarlo/LoadLibs.C
  MonteCarlo/rungen.C
  MonteCarlo/runTask.C
  )

install (FILES
  ${INPUT_FILES}
  MonteCarlo/parse_opts.sh
  MonteCarlo/prepare_runscript.sh
  DESTINATION PWG/HMTF/MonteCarlo)

# Install executable scripts
install (FILES
  MonteCarlo/gen_dipsy.sh
  MonteCarlo/runmc
  MonteCarlo/run.sh
  MonteCarlo/validation.sh
  DESTINATION PWG/HMTF/MonteCarlo
  PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
