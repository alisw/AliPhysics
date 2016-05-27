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
set(MODULE AliHLTGlobal)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/HLT/global
                    ${AliRoot_SOURCE_DIR}/HLT/global/physics
                   )

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/ITS/ITSrec
                    ${AliRoot_SOURCE_DIR}/HLT/BASE
                    ${AliRoot_SOURCE_DIR}/HLT/BASE/util
                    ${AliRoot_SOURCE_DIR}/HLT/ITS
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib
                    ${AliRoot_SOURCE_DIR}/HLT/CALO
                    ${AliRoot_SOURCE_DIR}/ITS/ITSbase
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCbase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCrec
                    ${AliRoot_SOURCE_DIR}/TPC/TPCcalib
                   )

# Sources in alphabetical order
set(SRCS
    AliHLTGlobalAgent.cxx
    AliHLTGlobalDCSPublisherComponent.cxx
    AliHLTGlobalEsdConverterComponent.cxx
    AliHLTGlobalHistoCollector.cxx
    AliHLTGlobalHistoComponent.cxx
    AliHLTGlobalOfflineVertexerComponent.cxx
    AliHLTGlobalPreprocessor.cxx
    AliHLTGlobalTrackMatcherComponent.cxx
    AliHLTGlobalTrackMatcher.cxx
    AliHLTGlobalTrackMergerComponent.cxx
    AliHLTGlobalTrackMerger.cxx
    AliHLTGlobalVertexerComponent.cxx
    AliHLTGlobalVertexerHistoComponent.cxx
    AliHLTPrimaryVertexFinderComponent.cxx
    AliHLTV0FinderComponent.cxx
    AliHLTVertexFinderBase.cxx
    AliHLTAsyncTestComponent.cxx
    AliHLTAsyncCalibrationComponent.cxx
    AliHLTZeroComponent.cxx
    AliHLTGlobalCompareFlatComponent.cxx
    AliHLTGlobalFlatEsdConverterComponent.cxx
    AliHLTGlobalEsdToFlatConverterComponent.cxx
    AliHLTGlobalFlatEsdTestComponent.cxx
    AliHLTRootObjectMergerComponent.cxx
    AliFlatESDEvent.cxx
    AliFlatESDEvent.h
    AliFlatESDFriend.cxx
    AliFlatESDFriend.h
    AliFlatESDFriendTrack.cxx
    AliFlatESDFriendTrack.h
    AliFlatESDTrack.cxx
    AliFlatESDTrack.h
    AliFlatESDTrigger.cxx
    AliFlatESDTrigger.h
    AliFlatESDVZERO.cxx
    AliFlatMultiplicity.cxx
    AliFlatESDVZEROFriend.cxx
    AliFlatESDV0.cxx
    AliFlatESDV0.h
    AliFlatESDVertex.cxx
    AliFlatESDVertex.h
    AliFlatExternalTrackParam.cxx
    AliFlatExternalTrackParam.h
    AliFlatTPCCluster.h
    AliFlatTPCseed.cxx
    AliFlatTPCseed.h
    AliHLTLumiRegComponent.cxx
    AliHLTGlobalPromptRecoQAComponent.cxx
    physics/AliHLTCaloHistoCellEnergy.cxx
    physics/AliHLTCaloHistoClusterEnergy.cxx
    physics/AliHLTCaloHistoComponent.cxx
    physics/AliHLTCaloHistoInvMass.cxx
    physics/AliHLTCaloHistoMatchedTracks.cxx
    physics/AliHLTCaloHistoProducer.cxx
    physics/AliHLTMultiplicityCorrelationsComponent.cxx
    physics/AliHLTMultiplicityCorrelations.cxx
    physics/AliHLTV0HistoComponent.cxx
    physics/AliHLTAnalysisManager.cxx
    physics/AliHLTAnalysisManagerComponent.cxx
    physics/examples/AliAnalysisTaskExampleV.cxx
   )

if(ZEROMQ_FOUND)
  include_directories(${ZEROMQ_INCLUDE_DIR})
  set(SRCS 
      ${SRCS}
      AliHLTZMQsink.cxx
      AliHLTZMQsource.cxx
     )
endif(ZEROMQ_FOUND)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER ESD CDB ANALYSIS ANALYSISalice HLTbase AliHLTUtil AliHLTITS AliHLTTPC TPCrec TPCbase ITSbase ITSrec TPCcalib AliHLTCalo)
if(ZEROMQ_FOUND)
  set(LIBDEPS ${LIBDEPS} ${ZEROMQ_LIBRARIES})
endif(ZEROMQ_FOUND)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
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

#install macros
install(DIRECTORY macros DESTINATION HLT/global)
install(DIRECTORY physics/macros DESTINATION HLT/global/physics/)

if(ZEROMQ_FOUND) 

  # Additional include folders in alphabetical order except ROOT
  include_directories(
                      ${AliRoot_SOURCE_DIR}/HLT/BASE/HOMER
                      ${AliRoot_SOURCE_DIR}/HLT/BASE/util
                      ${ZEROMQ_INCLUDE_DIR}
                      )

  add_executable(ZMQDIMlumiregServer ZMQDIMlumiregServer.cxx)
  target_link_libraries(ZMQDIMlumiregServer HLTbase AliHLTGlobal AliHLTUtil AliHLTHOMER RAWDatabase Geom Graf MathCore Net Tree EG Gpad Matrix Minuit Physics VMC Thread STEERBase XMLParser Graf3d RIO Hist Core zmq)
  
  # Installation of 
  install(TARGETS ZMQDIMlumiregServer RUNTIME DESTINATION bin)
endif(ZEROMQ_FOUND)
