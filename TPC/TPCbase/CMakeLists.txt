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
set(MODULE TPCbase)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TPC/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS} 
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/dateStream
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCrec
                    ${CMAKE_BINARY_DIR}/version
                   )

# Sources in alphabetical order
set(SRCS
    AliDetectorParam.cxx
    AliDigitsArray.cxx
    AliDigits.cxx
    AliH2F.cxx
    AliSegmentArray.cxx
    AliSegmentID.cxx
    AliSimDigits.cxx
    AliTPCAltroMapping.cxx
    AliTPCBoundaryVoltError.cxx
    AliTPCCalibCE.cxx
    AliTPCcalibDB.cxx
    AliTPCcalibDButil.cxx
    AliTPCCalibGlobalMisalignment.cxx
    AliTPCCalibPedestal.cxx
    AliTPCCalibPulser.cxx
    AliTPCCalibRawBase.cxx
    AliTPCCalibRaw.cxx
    AliTPCCalibVdrift.cxx
    AliTPCCalibViewer.cxx
    AliTPCCalibViewerGUI.cxx
    AliTPCCalibViewerGUItime.cxx
    AliTPCCalPad.cxx
    AliTPCCalROC.cxx
    AliTPCComposedCorrection.cxx
    AliTPCConfigDA.cxx
    AliTPCConfigParser.cxx
    AliTPCCorrection.cxx
    AliTPCCorrectionDrift.cxx
    AliTPCCorrectionLookupTable.cxx
    AliTPCdataQA.cxx
    AliTPCDigitsArray.cxx
    AliTPCExBBShape.cxx
    AliTPCExB.cxx
    AliTPCExBEffective.cxx
    AliTPCExBEffectiveSector.cxx
    AliTPCExBExact.cxx
    AliTPCExBFirst.cxx
    AliTPCExBTwist.cxx
    AliTPCFCVoltError3D.cxx
    AliTPCGGVoltError.cxx
    AliTPCInverseCorrection.cxx
    AliTPCLaserTrack.cxx
    AliTPCLoader.cxx
    AliTPCmapper.cxx
    AliTPCParam.cxx
    AliTPCParamSR.cxx
    AliTPCPointCorrection.cxx
    AliTPCPreprocessor.cxx
    AliTPCPreprocessorOnline.cxx
    AliTPCPRF2D.cxx
    AliTPCQAChecker.cxx
    AliTPCRawStreamV3.cxx
    AliTPCRF1D.cxx
    AliTPCROC.cxx
    AliTPCROCVoltError3D.cxx
    AliTPCSensorTempArray.cxx
    AliTPCSensorTemp.cxx
    AliTPCSpaceCharge3D.cxx
    AliTPCSpaceCharge.cxx
    AliTPCTempMap.cxx
    AliTPCTransform.cxx
    AliTransform.cxx
    AliXRDPROOFtoolkit.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER CDB ESD RAWDatabase RAWDatarec Minuit Gui)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")


# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
target_link_libraries(${MODULE} STEERBase STEER CDB ESD RAWDatabase RAWDatarec Minuit Gui)

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