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
set(MODULE AliHLTEMCAL)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/HLT/EMCAL
                   )

# Additional include folders in alphabetical order
include_directories(${AliRoot_SOURCE_DIR}/EMCAL/EMCALbase
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALraw
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALUtils
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALTriggerBase
                    ${AliRoot_SOURCE_DIR}/HLT/CALO
                    ${AliRoot_SOURCE_DIR}/HLT/BASE
                    ${AliRoot_SOURCE_DIR}/HLT/BASE/util
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    AliHLTEMCALAgent.cxx
    AliHLTEMCALCaloCells.cxx
    AliHLTEMCALClusterizerComponent.cxx
    AliHLTEMCALClusterizerComponentNbyN.cxx
    AliHLTEMCALClusterMonitorComponent.cxx
    AliHLTEMCALClusterMonitor.cxx
    AliHLTEMCALConstants.cxx
    AliHLTEMCALDigitMakerComponent.cxx
    AliHLTEMCALMapper.cxx
    AliHLTEMCALRawAnalyzerComponent.cxx
    AliHLTEMCALRawAnalyzerStandardComponent.cxx
    AliHLTEMCALRawAnalyzerCrudeComponent.cxx
    AliHLTEMCALRawAnalyzerFastFitComponent.cxx
    AliHLTEMCALRawAnalyzerLMSComponent.cxx
    AliHLTEMCALRawAnalyzerNNComponent.cxx
    AliHLTEMCALRawAnalyzerPeakFinderComponent.cxx
    AliHLTEMCALRawHistoMaker.cxx
    AliHLTEMCALRawHistoMakerComponent.cxx
    AliHLTEMCALRecoParamHandler.cxx
    AliHLTEMCALTriggerDataMakerComponent.cxx
    AliHLTEMCALTriggerMaker.cxx
    AliHLTEMCALTriggerMakerComponent.cxx
    AliHLTEMCALTRURawDigitMaker.cxx
    AliHLTEMCALRawAnalyzerComponentTRU.cxx
    AliHLTEMCALSTURawDigitMaker.cxx
    AliHLTEMCALRawAnalyzerComponentSTU.cxx
    AliHLTEMCALTriggerQAComponent.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
set(HDRS ${HDRS}
    AliHLTEMCALSTUHeaderStruct.h
)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER CDB HLTbase AliHLTUtil AliHLTCalo EMCALbase EMCALUtils EMCALraw EMCALTriggerBase)
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
