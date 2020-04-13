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
set(MODULE PWGCFFemtoDream)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGCF/FEMTOSCOPY/FemtoDream
                    ${AliPhysics_SOURCE_DIR}/PWGGA/GammaConvBase
                    ${AliPhysics_SOURCE_DIR}/PWGGA/Hyperon
                    ${AliPhysics_SOURCE_DIR}/PWG/DevNanoAOD)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY/
	)
		  

# Sources - alphabetical order
set(SRCS
  AliFemtoDreamEvent.cxx  
  AliFemtoDreamEventHist.cxx
  AliFemtoDreamEventCuts.cxx
  AliFemtoDreamBasePart.cxx
  AliFemtoDreamTrack.cxx
  AliFemtoDreamTrackHist.cxx
  AliFemtoDreamTrackMCHist.cxx 
  AliFemtoDreamTrackCuts.cxx
  AliFemtoDreamv0.cxx
  AliFemtoDreamv0Hist.cxx
  AliFemtoDreamv0MCHist.cxx
  AliFemtoDreamv0Cuts.cxx
  AliFemtoDreamCascade.cxx
  AliFemtoDreamCascadeHist.cxx
  AliFemtoDreamCascadeCuts.cxx
  AliFemtoDreamControlSample.cxx
  AliFemtoDreamPairCleanerHists.cxx
  AliFemtoDreamPairCleaner.cxx 
  AliFemtoDreamCollConfig.cxx 
  AliFemtoDreamCorrHists.cxx 
  AliFemtoDreamPartContainer.cxx 
  AliFemtoDreamZVtxMultContainer.cxx 
  AliFemtoDreamPartCollection.cxx 
  AliFemtoDreamAnalysis.cxx 
  AliAnalysisTaskFemtoDream.cxx
  AliAnalysisTaskFemtoDreamDeuteron.cxx
  AliAnalysisTaskNanoPt.cxx
  AliAnalysisTaskNanoXioton.cxx
  AliAnalysisTaskAODXioton.cxx
  AliAnalysisTaskNanoXioton1530.cxx
  AliAnalysisTaskNanoLoton.cxx
  AliAnalysisTaskLeuteronAOD.cxx
  AliAnalysisTaskLeuteronNanoAOD.cxx
  AliAnalysisTaskAODLoton.cxx
  AliAnalysisTaskNanoBBar.cxx  
  AliAnalysisTaskFemtoDreamPhi.cxx
  AliAnalysisTaskNanoAODFemtoDreamPhi.cxx
  AliAnalysisTaskFemtoDreamSigPi.cxx
  AliAnalysisTaskNanoAODFemtoDreamSigPi.cxx
  AliAnalysisTaskFemtoDreamPion.cxx
  AliAnalysisTaskNanoAODFemtoDreamPion.cxx
  AliAnalysisTaskNanoAODSigma0Femto.cxx
  AliAnalysisTaskNanoSigmaPlus.cxx
  AliAnalysisTaskNanoLD.cxx
  AliAnalysisTaskGrandma.cxx
  AliAnalysisTaskOtonOmega.cxx
  AliOtonOmegaAnalysis.cxx
  AliOtonOmegaCascadeCuts.cxx
  AliOtonOmegaCascade.cxx
  AliFemtoDreamHigherPairMath.cxx
  AliAnalysisTaskNanoLX.cxx
  AliAnalysisTaskOtonOmegaNanoAOD.cxx
  AliAnalysisTaskGeorgiosNTuple.cxx
  AliAnalysisTaskPOmegaPenne.cxx
  AliFemtoDreamBaseDump.cxx
  
  )


# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES Core EG GenVector Geom Gpad Hist MathCore Matrix Net Physics RIO Tree)
set(ALIROOT_DEPENDENCIES ANALYSIS ANALYSISalice AOD PWGGAGammaConvBase PWGGAHyperon)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Create an object to be reused in case of static libraries
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)

# Add a library to the project using the object
add_library_tested(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)

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


# generate_par("Hyperon")

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)

