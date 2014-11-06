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
set(MODULE TOFbase)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TOF/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatasim
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TOF/TOFrec
                    ${AliRoot_SOURCE_DIR}/TOF/TOFsim
                   )

# Sources in alphabetical order
set(SRCS
    AliTOFAlignment.cxx
    AliTOFArray.cxx
    AliTOFCableLengthMap.cxx
    AliTOFcalib.cxx
    AliTOFcalibHisto.cxx
    AliTOFChainSummaryData.cxx
    AliTOFChannelOffline.cxx
    AliTOFChannelOnlineArray.cxx
    AliTOFChannelOnline.cxx
    AliTOFChannelOnlineStatusArray.cxx
    AliTOFChannelOnlineStatus.cxx
    AliTOFCTPLatency.cxx
    AliTOFDaConfigHandler.cxx
    AliTOFDataDCS.cxx
    AliTOFDCSmaps.cxx
    AliTOFDecoder.cxx
    AliTOFDecoderSummaryData.cxx
    AliTOFDecoderV2.cxx
    AliTOFDeltaBCOffset.cxx
    AliTOFdigit.cxx
    AliTOFDigitMap.cxx
    AliTOFDRMSummaryData.cxx
    AliTOFFEEDump.cxx
    AliTOFFEEReader.cxx
    AliTOFFormatDCS.cxx
    AliTOFGeometry.cxx
    AliTOFHitDataBuffer.cxx
    AliTOFHitData.cxx
    AliTOFHitField.cxx
    AliTOFLTMSummaryData.cxx
    AliTOFLvHvDataPoints.cxx
    AliTOFNoiseConfigHandler.cxx
    AliTOFPreprocessor.cxx
    AliTOFPreprocessorFDR.cxx
    AliTOFrawData.cxx
    AliTOFRawMap.cxx
    AliTOFRawStream.cxx
    AliTOFReadoutInfo.cxx
    AliTOFResponseParams.cxx
    AliTOFRunParams.cxx
    AliTOFT0Fill.cxx
    AliTOFT0FillOnlineCalib.cxx
    AliTOFTDCErrorBuffer.cxx
    AliTOFTDCError.cxx
    AliTOFTDCHitBuffer.cxx
    AliTOFTDCHit.cxx
    AliTOFTRMSummaryData.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER ESD CDB RAWDatabase RAWDatarec RAWDatasim TOFrec TOFsim)
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