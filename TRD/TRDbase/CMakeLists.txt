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
set(MODULE TRDbase)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TRD/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TRD/TRDrec
                   )

# Sources in alphabetical order
set(SRCS
    AliTRDalignment.cxx
    AliTRDarrayADC.cxx
    AliTRDarrayDictionary.cxx
    AliTRDarraySignal.cxx
    AliTRDCalChamberStatus.cxx
    AliTRDCalDCS.cxx
    AliTRDCalDCSFEE.cxx
    AliTRDCalDCSFEEv2.cxx
    AliTRDCalDCSGTUBoardInfo.cxx
    AliTRDCalDCSGTUCtpOpc.cxx
    AliTRDCalDCSGTU.cxx
    AliTRDCalDCSGTUSegment.cxx
    AliTRDCalDCSGTUTgu.cxx
    AliTRDCalDCSGTUTmu.cxx
    AliTRDCalDCSPTR.cxx
    AliTRDCalDCSv2.cxx
    AliTRDCalDet.cxx
    AliTRDCalibChamberStatus.cxx
    AliTRDcalibDB.cxx
    AliTRDCalibPadStatus.cxx
    AliTRDCalibraExbAltFit.cxx
    AliTRDCalibraFit.cxx
    AliTRDCalibraMode.cxx
    AliTRDCalibraVdriftLinearFit.cxx
    AliTRDCalibraVector.cxx
    AliTRDCalibViewer.cxx
    AliTRDCalibViewerGUI.cxx
    AliTRDCalMonitoring.cxx
    AliTRDCalOnlineGainTable.cxx
    AliTRDCalOnlineGainTableMCM.cxx
    AliTRDCalOnlineGainTableROC.cxx
    AliTRDCalPad.cxx
    AliTRDCalPadStatus.cxx
    AliTRDCalPID.cxx
    AliTRDCalROC.cxx
    AliTRDCalSingleChamberStatus.cxx
    AliTRDCalTrapConfig.cxx
    AliTRDCalTrkAttach.cxx
    AliTRDCommonParam.cxx
    AliTRDDataDCS.cxx
    AliTRDdigit.cxx
    AliTRDdigitsManager.cxx
    AliTRDdigitsParam.cxx
    AliTRDEntriesInfo.cxx
    AliTRDfeeParam.cxx
    AliTRDgeometry.cxx
    AliTRDgtuParam.cxx
    AliTRDgtuSim.cxx
    AliTRDgtuTMU.cxx
    AliTRDhit.cxx
    AliTRDltuParam.cxx
    AliTRDmcmSim.cxx
    AliTRDonlineTrackMatching.cxx
    AliTRDpadPlane.cxx
    AliTRDPhInfo.cxx
    AliTRDPreprocessor.cxx
    AliTRDPrfInfo.cxx
    AliTRDQAChecker.cxx
    AliTRDrawData.cxx
    AliTRDrawStream.cxx
    AliTRDrecoParam.cxx
    AliTRDSaxHandler.cxx
    AliTRDSensorArray.cxx
    AliTRDSensor.cxx
    AliTRDSignalIndex.cxx
    AliTRDSimParam.cxx
    AliTRDtrackGTU.cxx
    AliTRDtrackletBase.cxx
    AliTRDtrackletGTU.cxx
    AliTRDtrackletMCM.cxx
    AliTRDtrackletWord.cxx
    AliTRDtrackOnline.cxx
    AliTRDtrapConfig.cxx
    AliTRDtrapConfigHandler.cxx
    AliTRDUshortInfo.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER ESD CDB RAWDatabase RAWDatarec TRDrec Gpad Gui)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")


# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "-O -g")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Check if CRS support was enabled
if(TRD_RAW_CRC)
    message(STATUS "CRC support for AliTRDrawStream enabled")
    add_definitions(-DTRD_RAW_CRC)
endif(TRD_RAW_CRC)

# Check if RAW debug enabled
if(TRD_RAW_DEBUG)
    message(STATUS "Debug support for AliTRDrawStream enabled")
    add_definitions(-DTRD_RAW_DEBUG)
endif(TRD_RAW_DEBUG)


# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)