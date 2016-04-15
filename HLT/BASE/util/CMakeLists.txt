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
set(MODULE AliHLTUtil)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/HLT/BASE/util
                   )

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/HLT/BASE
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/AOD
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    AliHLTAgentUtil.cxx
    AliHLTBlockFilterComponent.cxx
    AliHLTCaloClusterReader.cxx
    AliHLTCompStatCollector.cxx
    AliHLTCorruptorComponent.cxx
    AliHLTDataGenerator.cxx
    AliHLTESDCaloClusterMaker.cxx
    AliHLTEsdCollectorComponent.cxx
    AliHLTESDMCEventPublisherComponent.cxx
    AliHLTESDTrackCuts.cxx
    AliHLTFilePublisher.cxx
    AliHLTFileWriter.cxx
    AliHLTFXSFileWriter.cxx
    AliHLTFXSWriterComponent.cxx
    AliHLTGlobalBarrelTrack.cxx
    AliHLTJets.cxx
    AliHLTLoaderPublisherComponent.cxx
    AliHLTMCEvent.cxx
    AliHLTMonitoringRelay.cxx
    AliHLTObjectCompressionComponent.cxx
    AliHLTOUTPublisherComponent.cxx
    AliHLTRawReaderPublisherComponent.cxx
    AliHLTReadoutListDumpComponent.cxx
    AliHLTRecoParamComponent.cxx
    AliHLTRootFilePublisherComponent.cxx
    AliHLTRootFileStreamerComponent.cxx
    AliHLTRootFileWriterComponent.cxx
    AliHLTRootSchemaEvolutionComponent.cxx
    AliHLTObjectRenamer.cxx
   )

if(ZEROMQ_FOUND) 
  include_directories(
                      ${AliRoot_SOURCE_DIR}/HLT/BASE/HOMER
                      ${AliRoot_SOURCE_DIR}/HLT/global
                      ${ZEROMQ_INCLUDE_DIR}
                      )

  set(SRCS
      ${SRCS}
      AliZMQhistViewer.cxx
      AliZMQMTviewerGUI.cxx
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
set(LIBDEPS STEERBase STEER AOD CDB ESD RAWDatabase RAWDatarec ANALYSISalice HLTbase Graf)
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

if(ZEROMQ_FOUND) 

  # Additional include folders in alphabetical order except ROOT
  include_directories(
                      ${AliRoot_SOURCE_DIR}/HLT/BASE/HOMER
                      ${AliRoot_SOURCE_DIR}/HLT/global
                      ${ZEROMQ_INCLUDE_DIR}
                      )

  add_executable(ZMQfileProxy ZMQfileProxy.cxx)
  target_link_libraries(ZMQfileProxy HLTbase AliHLTHOMER AliHLTGlobal RAWDatabase Geom Graf MathCore Net Tree EG Gpad Matrix Minuit Physics VMC Thread STEERBase XMLParser Graf3d RIO Hist Core)
  add_executable(ZMQDIMlumiregServer ZMQDIMlumiregServer.cxx)
  target_link_libraries(ZMQDIMlumiregServer HLTbase AliHLTHOMER AliHLTGlobal RAWDatabase Geom Graf MathCore Net Tree EG Gpad Matrix Minuit Physics VMC Thread STEERBase XMLParser Graf3d RIO Hist Core)
  add_executable(ZMQHOMERproxy ZMQHOMERproxy.cxx)
  target_link_libraries(ZMQHOMERproxy HLTbase AliHLTHOMER AliHLTGlobal RAWDatabase Geom Graf MathCore Net Tree EG Gpad Matrix Minuit Physics VMC Thread STEERBase XMLParser Graf3d RIO Hist Core)
  add_executable(ZMQROOTmerger ZMQROOTmerger.cxx)
  target_link_libraries(ZMQROOTmerger HLTbase AliHLTGlobal RAWDatabase Geom Graf MathCore Net Tree EG Gpad Matrix Minuit Physics VMC Thread STEERBase XMLParser Graf3d RIO Hist Core)
  add_executable(ZMQhistSource ZMQhistSource.cxx)
  target_link_libraries(ZMQhistSource HLTbase AliHLTGlobal RAWDatabase Geom Graf MathCore Net Tree EG Gpad Matrix Minuit Physics VMC Thread STEERBase XMLParser Graf3d RIO Hist Core)
  add_executable(ZMQhistViewer ZMQhistViewer.cxx)
  target_link_libraries(ZMQhistViewer HLTbase AliHLTGlobal AliHLTUtil RAWDatabase Geom Graf MathCore Net Tree EG Gpad Matrix Minuit Physics VMC Thread STEERBase XMLParser Graf3d RIO Hist Core Cint Hist Rint)
  add_executable(ZMQproxy ZMQproxy.cxx)
  target_link_libraries(ZMQproxy HLTbase AliHLTGlobal RAWDatabase Geom Graf MathCore Net Tree EG Gpad Matrix Minuit Physics VMC Thread STEERBase XMLParser Graf3d RIO Hist Core)
  add_executable(ZMQHLTchain ZMQHLTchain.cxx)
  target_link_libraries(ZMQHLTchain HLTbase AliHLTGlobal RAWDatabase Geom Graf MathCore Net Tree EG Gpad Matrix Minuit Physics VMC Thread STEERBase XMLParser Graf3d RIO Hist Core)
  
  # Installation of 
  install(TARGETS ZMQfileProxy ZMQHOMERproxy ZMQROOTmerger ZMQhistSource ZMQhistViewer ZMQproxy ZMQHLTchain ZMQDIMlumiregServer RUNTIME DESTINATION bin)
  install(PROGRAMS ZMQmonitor.py DESTINATION bin)
endif(ZEROMQ_FOUND)
