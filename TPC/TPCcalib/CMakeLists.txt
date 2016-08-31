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
set(MODULE TPCcalib)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TPC/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/RAW/dateStream
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STAT
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCbase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCrec
                   )

# Sources in alphabetical order
set(SRCS
    AliAnalysisTaskPt.cxx
    AliTPCAnalysisTaskcalib.cxx
    AliTPCcalibAlign.cxx
    AliTPCcalibBase.cxx
    AliTPCcalibCalib.cxx
    AliTPCcalibCosmic.cxx
    AliTPCcalibGainMult.cxx
    AliTPCCalibKr.cxx
    AliTPCCalibKrTask.cxx
    AliTPCcalibLaser.cxx
    AliTPCcalibSummary.cxx
    AliTPCcalibTime.cxx
    AliTPCcalibTimeGain.cxx
    AliTPCcalibTracksCuts.cxx
    AliTPCcalibTracks.cxx
    AliTPCcalibTracksGain.cxx
    AliTPCcalibV0.cxx
    AliTPCCalPadRegion.cxx
    AliTPCCorrectionFit.cxx
    AliTPCDcalibRes.cxx
    AliTPCFitPad.cxx
    AliTPCkalmanAlign.cxx
    AliTPCMisAligner.cxx
    AliTPCcalibAlignInterpolation.cxx
    AliTPCPreprocessorOffline.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase CDB ESD STEER ANALYSIS ANALYSISalice STAT RAWDatarec RAWDatabase TPCbase TPCrec Core EG Geom Gpad Graf Hist MathCore Matrix Minuit Postscript Proof RIO Tree)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
target_link_libraries(${MODULE} ${LIBDEPS})

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# Calib maps
install(DIRECTORY maps DESTINATION TPC/${MODULE})

install(FILES ${HDRS} DESTINATION include)
