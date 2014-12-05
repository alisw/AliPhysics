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
set(MODULE PWGDQdielectron)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PWGDQ/dielectron)

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/CORRFW/
                    ${AliRoot_SOURCE_DIR}/PWGCF/Correlations # what deps here?
                    ${AliRoot_SOURCE_DIR}/PWG/FLOW/Base
                    ${AliRoot_SOURCE_DIR}/PWG/FLOW/Tasks
                    ${AliRoot_SOURCE_DIR}/PWGLF/FORWARD
                    ${AliRoot_SOURCE_DIR}/STEER/AOD
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TRD/TRDbase
                   )

# Sources in alphabetical order
set(SRCS
      AliAnalysisTaskDielectronEfficiency.cxx
      AliAnalysisTaskDielectronFilter.cxx
      AliAnalysisTaskDielectronME.cxx
      AliAnalysisTaskDielectronReadAODBranch.cxx
      AliAnalysisTaskDielectronSE.cxx
      AliAnalysisTaskMultiDielectron.cxx
      AliAnalysisTaskMultiDielectronTG.cxx
      AliAnalysisTaskReducedTree.cxx
      AliAnalysisTaskSingleParticle.cxx
      AliDielectronBtoJPSItoEleCDFfitFCN.cxx
      AliDielectronBtoJPSItoEleCDFfitFCNfitter.cxx
      AliDielectronBtoJPSItoEleCDFfitHandler.cxx
      AliDielectronBtoJPSItoEle.cxx
      AliDielectronCF.cxx
      AliDielectronCFdraw.cxx
      AliDielectronCutGroup.cxx
      AliDielectronCutQA.cxx
      AliDielectron.cxx
      AliDielectronDebugTree.cxx
      AliDielectronEventCuts.cxx
      AliDielectronEvent.cxx
      AliDielectronHelper.cxx
      AliDielectronHF.cxx
      AliDielectronHFhelper.cxx
      AliDielectronHistos.cxx
      AliDielectronMC.cxx
      AliDielectronMixingHandler.cxx
      AliDielectronPair.cxx
      AliDielectronPairLegCuts.cxx
      AliDielectronPID.cxx
      AliDielectronSignalBase.cxx
      AliDielectronSignalExt.cxx
      AliDielectronSignalFunc.cxx
      AliDielectronSignalMC.cxx
      AliDielectronSpectrum.cxx
      AliDielectronTrackCuts.cxx
      AliDielectronTrackRotator.cxx
      AliDielectronTRDpidCut.cxx
      AliDielectronV0Cuts.cxx
      AliDielectronVarCuts.cxx
      AliDielectronVarManager.cxx
      AliReducedEvent.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES Core EG Gpad Graf Hist MathCore Matrix Minuit Net Physics RIO Tree)
set(ALIROOT_DEPENDENCIES ANALYSIS ANALYSISalice AOD ESD PWGflowTasks PWGflowBase STEERBase TRDbase PWGLFforward2)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)
# Add a library to the project using the object
add_library(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)
target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

if(DATE_FOUND)
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS}")
endif(DATE_FOUND)

# Additional compilation and linking flags
set(MODULE_COMPILE_FLAGS " ${MODULE_COMPILE_FLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Setting compilation flags for the object
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
# Setting the linking flags for the library
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)