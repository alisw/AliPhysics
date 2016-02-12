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
set(MODULE JETAN)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/JETAN/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                   )

# Sources in alphabetical order
set(SRCS
    AliAnalysisTaskDiJets.cxx
    AliAnalysisTaskFastEmbedding.cxx
    AliAnalysisTaskJetsFinder.cxx
    AliAnalysisTaskJetsReader.cxx
    AliAnalysisTaskKMeans.cxx
    AliCdfJetFinder.cxx
    AliCdfJetHeader.cxx
    AliDAJetFinder.cxx
    AliDAJetHeader.cxx
    AliEventShape.cxx
    AliJetBkg.cxx
    AliJetCalTrk.cxx
    AliJetFillCalTrkEvent.cxx
    AliJetFillCalTrkTrack.cxx
    AliJetFillCalTrkTrackMC.cxx
    AliJetFinder.cxx
    AliJetHeader.cxx
    AliJetHistos.cxx
    AliJetKineReaderHeader.cxx
    AliJetReader.cxx
    AliJetReaderHeader.cxx
    AliKMeansClustering.cxx
    AliUA1JetFinder.cxx
    AliUA1JetHeaderV1.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase ESD AOD ANALYSIS ANALYSISalice EMCALUtils EG Eve Net Smatrix)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)
