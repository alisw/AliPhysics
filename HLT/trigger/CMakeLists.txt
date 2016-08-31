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
set(MODULE AliHLTTrigger)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/HLT/trigger)

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/HLT/BASE
                    ${AliRoot_SOURCE_DIR}/HLT/BASE/util
                    ${AliRoot_SOURCE_DIR}/HLT/ITS
                    ${AliRoot_SOURCE_DIR}/HLT/MUON
                    ${AliRoot_SOURCE_DIR}/HLT/TRD
                    ${AliRoot_SOURCE_DIR}/MUON/MUONraw
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/AOD
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCcalib
                    ${AliRoot_SOURCE_DIR}/TRD/TRDbase
                    ${AliRoot_SOURCE_DIR}/TRD/TRDrec
                   )

# Sources in alphabetical order
set(SRCS
    AliHLTD0Candidate.cxx
    AliHLTD0toKpi.cxx
    AliHLTD0Trigger.cxx
    AliHLTEmcalElectronMonitorComponent.cxx
    AliHLTEmcalElectronMonitor.cxx
    AliHLTFastJetMonitorComponent.cxx
    AliHLTFastJetMonitor.cxx
    AliHLTGlobalTriggerComponent.cxx
    AliHLTGlobalTriggerConfig.cxx
    AliHLTGlobalTrigger.cxx
    AliHLTGlobalTriggerWrapper.cxx
    AliHLTMuonSpectroScalars.cxx
    AliHLTMuonSpectroTriggerComponent.cxx
    AliHLTTRDTriggerComponent.cxx
    AliHLTTriggerAgent.cxx
    AliHLTTriggerBarrelGeomMultiplicity.cxx
    AliHLTTriggerBarrelMultiplicity.cxx
    AliHLTTriggerCaloClusterEnergy.cxx
    AliHLTTriggerCaloClusterEnergy.cxx
    AliHLTTriggerCosmics.cxx
    AliHLTTriggerCounterComponent.cxx
    AliHLTTriggerCounters.cxx
    AliHLTTriggerDecisionParameters.cxx
    AliHLTTriggerDetectorGeom.cxx
    AliHLTTriggerDetectorGeomRectangle.cxx
    AliHLTTriggerEmcalClusterEnergy.cxx
    AliHLTTriggerEmcalElectron.cxx
    AliHLTTriggerGammaConversion.cxx
    AliHLTTrigger.cxx
    AliHLTTriggerITSMultiplicity.cxx
    AliHLTTriggerPhosClusterEnergy.cxx
    AliHLTTriggerPhosMip.cxx
    AliHLTTriggerTrdClusterMultiplicity.cxx
    AliHLTUpcTriggerComponent.cxx
    AliHLTTriggerCTPRE.cxx
   )

# TO DO : Check if this is really needed
#if(FASTJET_FOUND)
#    include_directories(SYSTEM ${FASTJET_INCLUDE_DIR})
#    link_directories(${FASTJET_LIBS_DIR})
#    add_definitions(${FASTJET_DEFINITIONS})
#
#    set(SRCS ${SRCS}
#        AliHLTTriggerFastJet.cxx
#       )
#endif()

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER ESD CDB AOD RAWDatabase RAWDatarec ANALYSISalice HLTbase AliHLTUtil AliHLTITS AliHLTMUON AliHLTTRD TRDbase TRDrec TPCcalib MUONraw GenVector Gpad Graf Smatrix)

# !! TO FIX !!
# Temporary removing the linking to fastjet libraries as it adds system libraries to the rootmaps
# Rootmaps should contain only libraries that can be found in the LD_LIBRARY_PATH
#if(FASTJET_FOUND)
#    set(LIBDEPS ${LIBDEPS} ${FASTJET_LIBS})
#endif(FASTJET_FOUND)

generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
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
