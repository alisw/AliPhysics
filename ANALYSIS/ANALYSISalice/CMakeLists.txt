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
set(MODULE ANALYSISalice)

# Module include folders in alphabetical order except ROOT
include_directories(${AliRoot_SOURCE_DIR}/ANALYSIS/${MODULE})

# Additional include directories
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISaliceBase
                    ${AliRoot_SOURCE_DIR}/STEER/AOD
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order 
set(SRCS 
    AliAnalysisCuts.cxx
    AliAnalysisFilter.cxx
    AliAnalysisTaskBadChunkID.cxx
    AliAnalysisTaskBaseLine.cxx
    AliAnalysisTaskMCParticleFilter.cxx
    AliAnalysisTaskME.cxx
    AliAnalysisTaskPIDCombined.cxx
    AliAnalysisTaskPIDqa.cxx
    AliAnalysisTaskPIDResponse.cxx
    AliAnalysisTaskSE.cxx
    AliAnalysisTaskStat.cxx
    AliAnalysisTaskTagCreator.cxx
    AliTrackSelectionFactory.cxx
    AliAODv0KineCuts.cxx
    AliAODITSsaTrackCuts.cxx
    AliAODTrackSelection.cxx
    AliESDpidCuts.cxx
    AliESDtrackCuts.cxx
    AliESDv0Cuts.cxx
    AliESDv0KineCuts.cxx
    AliESDTrackSelection.cxx
    AliEventPoolLoop.cxx
    AliEventPoolManager.cxx
    AliEventPoolOTF.cxx
    AliEventPoolSparse.cxx
    AliHEPDataParser.cxx
    AliKineTrackCuts.cxx
    AliMultiEventInputHandler.cxx
    AliPIDResponseInputHandler.cxx
    AliTagAnalysis.cxx
    AliUnfolding.cxx
    AliVZEROEPSelectionTask.cxx
    AliVTrackSelection.cxx
    AliXMLCollection.cxx
  )


# If ROOT has no xml support desabling collections and AliEn support
if(ROOT_HASXML)
    set(SRCS ${SRCS}
         AliTagAnalysis.cxx
         AliXMLCollection.cxx
       )
else()
    message(STATUS "${MODULE} - ROOT was built with no XML support. Disabling XML collection and AliEn support")
endif()

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSIS ANALYSISaliceBase AOD ESD STEERBase Core EG Gpad Gui Hist MathCore Matrix Net Physics RIO TreePlayer Tree XMLIO)
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
