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
set(MODULE ESD)

# Module include folders
include_directories(${CMAKE_SOURCE_DIR}/STEER/${MODULE})

# Additional include directories - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${CMAKE_SOURCE_DIR}/STEER/STEERBase)

# Sources - alphabetical order
set(SRCS
    AliCascadeVertexer.cxx
    AliESDACORDE.cxx
    AliESDAD.cxx
    AliESDCaloCells.cxx
    AliESDCaloCluster.cxx
    AliESDCaloTrigger.cxx
    AliESDcascade.cxx
    AliESDCosmicTrack.cxx
    AliESD.cxx
    AliESDEvent.cxx
    AliESDFIT.cxx
    AliESDFMD.cxx
    AliESDfriend.cxx
    AliESDfriendTrack.cxx
    AliESDHandler.cxx
    AliESDHeader.cxx
    AliESDHLTDecision.cxx
    AliESDHLTtrack.cxx
    AliESDInputHandler.cxx
    AliESDInputHandlerRP.cxx
    AliESDkink.cxx
    AliESDMultITS.cxx
    AliESDMuonCluster.cxx
    AliESDMuonGlobalTrack.cxx
    AliESDMuonPad.cxx
    AliESDMuonTrack.cxx
    AliESDpid.cxx
    AliESDPmdTrack.cxx
    AliESDRun.cxx
    AliESDTOFCluster.cxx
    AliESDTOFHit.cxx
    AliESDTOFMatch.cxx
    AliESDtrack.cxx
    AliESDTrdTrack.cxx
    AliESDTrdTracklet.cxx
    AliESDTrdTrigger.cxx
    AliESDTZERO.cxx
    AliESDTZEROfriend.cxx
    AliESDUtils.cxx
    AliESDv0.cxx
    AliESDV0Params.cxx
    AliESDVertex.cxx
    AliESDVZERO.cxx
    AliESDVZEROfriend.cxx
    AliESDZDC.cxx
    AliExpression.cxx
    AliFMDFloatMap.cxx
    AliFMDMap.cxx
    AliKalmanTrack.cxx
    AliKFParticleBase.cxx
    AliKFParticle.cxx
    AliKFVertex.cxx
    AliMeanVertex.cxx
    AliMultiplicity.cxx
    AliRawDataErrorLog.cxx
    AliSelector.cxx
    AliStrLine.cxx
    AliTrackerBase.cxx
    AliTrackPointArray.cxx
    AliTriggerBCMask.cxx
    AliTriggerClass.cxx
    AliTriggerCluster.cxx
    AliTriggerConfiguration.cxx
    AliTriggerDescriptor.cxx
    AliTriggerInput.cxx
    AliTriggerInteraction.cxx
    AliTriggerIR.cxx
    AliTriggerPFProtection.cxx
    AliTriggerScalersESD.cxx
    AliTriggerScalersRecordESD.cxx
    AliV0vertexer.cxx
    AliVertex.cxx
    AliVertexerTracks.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "-O -g")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)
