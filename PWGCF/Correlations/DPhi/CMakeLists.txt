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
set (MODULE PWGCFCorrelationsDPhi)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PWGCF/Correlations/DPhi)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ESDfilter
                    ${AliRoot_SOURCE_DIR}/CORRFW
                    ${AliRoot_SOURCE_DIR}/EVGEN
                    ${AliRoot_SOURCE_DIR}/LHAPDF
                    ${AliRoot_SOURCE_DIR}/PWG/muon
                    ${AliRoot_SOURCE_DIR}/PWG/Tools
                    ${AliRoot_SOURCE_DIR}/PWGCF/Correlations/Base
                    ${AliRoot_SOURCE_DIR}/PWGLF/SPECTRA/PiKaPr/TestAOD
                    ${AliRoot_SOURCE_DIR}/PYTHIA6/AliPythia6
                    ${AliRoot_SOURCE_DIR}/PYTHIA8/AliPythia8
                    ${AliRoot_SOURCE_DIR}/PYTHIA8/pythia8175/include
                    ${AliRoot_SOURCE_DIR}/STEER/AOD
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
  )

# Sources - alphabetical order
set(SRCS
    AliAnalysisTask3PCorrelations.cxx
    AliAnalysisTaskContMC.cxx
    AliAnalysisTaskDiHadron.cxx
    AliAnalysisTaskDiJetCorrelations.cxx
    AliAnalysisTaskDptDptCorrelations.cxx
    AliAnalysisTaskDptDptQA.cxx
    AliAnalysisTaskLeadingTrackUE.cxx
    AliAnalysisTaskLongRangeCorrelations.cxx
    AliAnalysisTaskMinijet.cxx
    AliAnalysisTaskPhiCorrelations.cxx
    AliAnalysisTaskpxpx.cxx
    AliAnalysisTaskpxpy.cxx
    AliAnalysisTaskpxpy.cxx
    AliAnalysisTaskpypy.cxx
    AliAnalysisTaskpypy.cxx
    AliAnalysisTaskPythiaMpi.cxx
    AliAnalysisTaskpzpz.cxx
    AliAnalysisTaskpzpz.cxx
    AliAnalysisTaskTwoPlusOne.cxx
    AliDptDptInMC.cxx
    AliLeadingV0Correlation.cxx
    AliPhiCorrelationsQATask.cxx
    DiHadronPID/AliAnalysisTaskCompareAODTrackCuts.cxx
    DiHadronPID/AliAnalysisTaskDiHadronPID.cxx
    DiHadronPID/AliAODEventCutsDiHadronPID.cxx
    DiHadronPID/AliAODTrackCutsDiHadronPID.cxx
    DiHadronPID/AliFunctionsDiHadronPID.cxx
    DiHadronPID/AliHistToolsDiHadronPID.cxx
    DiHadronPID/AliTrackDiHadronPID.cxx
    FourierDecomposition/AliDhcTask.cxx
    FourierDecomposition/AliMCTruthCent.cxx
    FourierDecomposition/AliMCTruthTrackMaker.cxx
    FourierDecomposition/AliMuonEffMC.cxx
    FourierDecomposition/AliPool.cxx
    MuonHadron/AliAnalysisTaskDiMuonCorrelations.cxx
    MuonHadron/AliAnalysisTaskMuonHadronCorrelations.cxx
    PidPid/AliAnalysisTaskPidPidCorrelations.cxx
    TriggerPID/AliAnalysisTaskPIDCORR.cxx
    TriggerPID/AliTwoParticlePIDCorr.cxx
    TriggerPID/AliTwoParticlePIDCorrKine.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS AliPythia8 ANALYSISalice CORRFW PWGmuon PWGCFCorrelationsBase PWGLFspectra)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

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

# Macros
install(FILES MuonHadron/AddAnalysisTaskMuonHadronCorrelations.C DESTINATION PWGCF/Correlations/DPhi/MuonHadron/)