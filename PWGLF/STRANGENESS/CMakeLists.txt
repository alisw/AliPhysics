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
# **************************************************************************/

#Module
set(MODULE PWGLFSTRANGENESS)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGLF/STRANGENESS)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/EVENTMIX
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWGCF/Correlations/Base
                    ${AliPhysics_SOURCE_DIR}/PWGLF/STRANGENESS/Cascades/lightvertexers
                    ${AliPhysics_SOURCE_DIR}/PWGUD/base
  )

# Sources - alphabetical order
set(SRCS
  Cascades/AliAnalysisTaskCheckCascade.cxx
  Cascades/AliAnalysisTaskCheckCascadePbPb.cxx
  Cascades/AliAnalysisTaskCheckPerformanceCascade.cxx
  Cascades/AliAnalysisTaskCheckPerformanceCascadePbPb.cxx
  LambdaK0PbPb/AliAnalysisCentralitySelector.cxx
  LambdaK0PbPb/AliAnalysisTaskCTauPbPb.cxx
  LambdaK0PbPb/AliAnalysisTaskCTauPbPbaod.cxx
  LambdaK0PbPb/AliAnalysisTaskLukeAOD.cxx
  LambdaK0PbPb/AliAnalysisTaskLukeV0.cxx
  LambdaK0PbPb/AliAnalysisTaskPerformanceStrange.cxx
  LambdaK0PbPb/AliAnalysisTaskV0ForRAA.cxx
  LambdaK0/AliAnalysisTaskExtractV0.cxx
  LambdaK0/AliAnalysisTaskExtractV0pPb.cxx
  LambdaK0/AliAnalysisTaskExtractV0AOD.cxx
  LambdaK0/AliAnalysisTaskExtractPerformanceV0.cxx
  LambdaK0/AliAnalysisTaskExtractPerformanceV0pPb.cxx
  LambdaK0PbPb/AliAnalysisTaskLK0Spectra.cxx
  Cascades/AliAnalysisTaskExtractCascade.cxx
  Cascades/AliAnalysisTaskExtractPerformanceCascade.cxx
  Cascades/AliAnalysisTaskStrangenessVsMultiplicity.cxx
  Cascades/AliAnalysisTaskStrangenessVsMultiplicityMC.cxx
  Cascades/AliAnalysisTaskPPVsMultCrossCheckMC.cxx
  Cascades/AliAnalysisTaskpANormalizationCheckMC.cxx
  Correlations/AliAnalysisTaskLambdaOverK0sJets.cxx
  Correlations/AliAnalysisTaskV0ChCorrelations.cxx
  Cascades/AliAnalysisTaskCheckCascadepp.cxx
  Cascades/AliAnalysisTaskCheckPerformanceCascadepp.cxx
  Cascades/lightvertexers/AliLightV0vertexer.cxx
  Cascades/lightvertexers/AliLightCascadeVertexer.cxx
  )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSISalice CORRFW OADB PWGUDbase)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)

# install macros
install(DIRECTORY Cascades/macros DESTINATION PWGLF/STRANGENESS/Cascades)
install(DIRECTORY Correlations/macros DESTINATION PWGLF/STRANGENESS/Correlations)
install(DIRECTORY LambdaK0/macros DESTINATION PWGLF/STRANGENESS/LambdaK0)
install(DIRECTORY LambdaK0PbPb DESTINATION  PWGLF/STRANGENESS)
install(DIRECTORY Cascades/corrections DESTINATION PWGLF/STRANGENESS/Cascades)
