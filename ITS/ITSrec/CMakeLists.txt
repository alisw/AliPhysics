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
set(MODULE ITSrec)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/ITS/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliRoot_SOURCE_DIR}/ITS/ITSbase
                    ${AliRoot_SOURCE_DIR}/ITS/ITSsim
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                   )

# Sources in alphabetical order
set(SRCS
    AliITSAlignMille2Constraint.cxx
    AliITSAlignMille2ConstrArray.cxx
    AliITSAlignMille2.cxx
    AliITSAlignMille2Module.cxx
    AliITSAlignMille.cxx
    AliITSAlignMilleData.cxx
    AliITSAlignMilleModule.cxx
    AliITSBadChannelsAuxSPD.cxx
    AliITSBadChannelsSPD.cxx
    AliITSChannelSPD.cxx
    AliITSclustererV2.cxx
    AliITSCorrectSDDPoints.cxx
    AliITSdEdxAnalyzer.cxx
    AliITSIOTrack.cxx
    AliITSMeanVertexer.cxx
    AliITSOnlineCalibrationSPD.cxx
    AliITSOnlineCalibrationSPDhandler.cxx
    AliITSOnlineSDDBase.cxx
    AliITSOnlineSDDCMN.cxx
    AliITSOnlineSDD.cxx
    AliITSOnlineSDDInjectors.cxx
    AliITSOnlineSDDTP.cxx
    AliITSOnlineSPDfoAnalyzer.cxx
    AliITSOnlineSPDfoChipConfig.cxx
    AliITSOnlineSPDfoChip.cxx
    AliITSOnlineSPDfo.cxx
    AliITSOnlineSPDfoInfo.cxx
    AliITSOnlineSPDHitArray.cxx
    AliITSOnlineSPDHitEvent.cxx
    AliITSOnlineSPDphysAnalyzer.cxx
    AliITSOnlineSPDphys.cxx
    AliITSOnlineSPDphysInfo.cxx
    AliITSOnlineSPDscanAnalyzer.cxx
    AliITSOnlineSPDscan.cxx
    AliITSOnlineSPDscanInfo.cxx
    AliITSOnlineSPDscanInfoMeanTh.cxx
    AliITSOnlineSPDscanInfoMultiple.cxx
    AliITSOnlineSPDscanMeanTh.cxx
    AliITSOnlineSPDscanMultiple.cxx
    AliITSOnlineSPDscanSingle.cxx
    AliITSPreprocessorSDD.cxx
    AliITSPreprocessorSPD.cxx
    AliITSPreprocessorSSD.cxx
    AliITSQASSDRefData.cxx
    AliITSRad.cxx
    AliITSRealignTracks.cxx
    AliITSreconstruction.cxx
    AliITSSumTP.cxx
    AliITSTPArrayFit.cxx
    AliITSVertexer3DTapan.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDECIES Core EG Geom Gpad Graf Hist MathCore Matrix Net RIO Tree)
set(ALIROOT_DEPENDECIES CDB ESD STEERBase STEER ITSbase ITSsim RAWDatabase RAWDatarec)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ALIROOT_DEPENDECIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ALIROOT_DEPENDENCIES})

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