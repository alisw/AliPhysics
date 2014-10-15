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
set(MODULE MUONgraphics)

# Module include folder
include_directories(${CMAKE_SOURCE_DIR}/MUON/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${CMAKE_SOURCE_DIR}/MUON/MUONbase
                    ${CMAKE_SOURCE_DIR}/MUON/MUONcalib
                    ${CMAKE_SOURCE_DIR}/MUON/MUONcore
                    ${CMAKE_SOURCE_DIR}/MUON/MUONgeometry
                    ${CMAKE_SOURCE_DIR}/MUON/MUONmapping
                    ${CMAKE_SOURCE_DIR}/MUON/MUONrec
                    ${CMAKE_SOURCE_DIR}/RAW/RAWDatabase
                    ${CMAKE_SOURCE_DIR}/RAW/RAWDatarec
                    ${CMAKE_SOURCE_DIR}/STEER/CDB
                    ${CMAKE_SOURCE_DIR}/STEER/ESD
                    ${CMAKE_SOURCE_DIR}/STEER/STEER
                    ${CMAKE_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    AliMUONAlignmentCompareDialog.cxx
    AliMUONAttPainter.cxx
    AliMUONAttPainterSelectorFrame.cxx
    AliMUONBusPatchPainter.cxx
    AliMUONChamberPainter.cxx
    AliMUONContour.cxx
    AliMUONContourHandler.cxx
    AliMUONContourMaker.cxx
    AliMUONContourMakerTest.cxx
    AliMUONContourPainter.cxx
    AliMUONDEPainter.cxx
    AliMUONManuContourMaker.cxx
    AliMUONManuPadPainter.cxx
    AliMUONManuPainter.cxx
    AliMUONMchViewApplication.cxx
    AliMUONNode.cxx
    AliMUONPainterColorSlider.cxx
    AliMUONPainterDataRegistry.cxx
    AliMUONPainterDataSourceFrame.cxx
    AliMUONPainterDataSourceItem.cxx
    AliMUONPainterEnv.cxx
    AliMUONPainterGroup.cxx
    AliMUONPainterHelper.cxx
    AliMUONPainterHighlighter.cxx
    AliMUONPainterInterfaceHelper.cxx
    AliMUONPainterMasterFrame.cxx
    AliMUONPainterMatrix.cxx
    AliMUONPainterMatrixFrame.cxx
    AliMUONPainterPadStore.cxx
    AliMUONPainterPlotSelector.cxx
    AliMUONPainterRegistry.cxx
    AliMUONPCBPainter.cxx
    AliMUONPointWithRef.cxx
    AliMUONPolygon.cxx
    AliMUONSegment.cxx
    AliMUONSegmentTree.cxx
    AliMUONTrackerDataCompareDialog.cxx
    AliMUONTrackerDataHistogrammer.cxx
    AliMUONTrackerDataWrapper.cxx
    AliMUONTrackerHV.cxx
    AliMUONVPainter.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER CDB ESD RAWDatabase RAWDatarec MUONbase MUONcore MUONmapping MUONcalib MUONrec MUONgeometry)
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