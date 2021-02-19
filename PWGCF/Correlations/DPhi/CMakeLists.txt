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
set(MODULE PWGCFCorrelationsDPhi)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGCF/Correlations/DPhi)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/PWG/muon
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWG/JETFW
                    ${AliPhysics_SOURCE_DIR}/PWG/DevNanoAOD
                    ${AliPhysics_SOURCE_DIR}/PWGCF/Correlations/Base
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWGLF/FORWARD/analysis2
		    ${AliPhysics_SOURCE_DIR}/EVENTMIX
  )

# Sources - alphabetical order
set(SRCS
    AliAnalysisTask3PCorrelations.cxx
    AliAnalysisTaskContMC.cxx
    AliAnalysisTaskCorPIDTOFQA.cxx
    AliAnalysisTaskCorPIDTOFdeut.cxx
    AliAnalysisTaskCorPIDTOFprot.cxx
    AliAnalysisTaskCorPIDTOFhadr.cxx
    AliAnalysisTaskCorPIDTOFdiprot.cxx
    AliAnalysisTaskDiHadCorrelHighPt.cxx
    AliAnalysisTaskDiHadron.cxx
    AliAnalysisTaskDptDptCorrelations.cxx
    AlianalysisTaskDptPIDpPb.cxx    
    AliAnalysisTaskDptDptQA.cxx
    AlianalysisTaskDptPID.cxx
    AliAnalysisTaskHadronPhiCorr.cxx
    AliAnalysisTaskLeadingTrackUE.cxx
    AliAnalysisTaskLongRangeCorrelations.cxx
    AliAnalysisTaskLongFluctuations2PC.cxx
    AliAnalysisTaskMinijet.cxx
    AliAnalysisTaskMinijetV2.cxx
    AliAnalysisTaskPhiCorrelations.cxx
    AliAnalysisTaskCorrelationsDev.cxx
    AliAnalysisTaskpxpx.cxx
    AliAnalysisTaskpxpy.cxx
    AliAnalysisTaskpxpy.cxx
    AliAnalysisTaskpypy.cxx
    AliAnalysisTaskpypy.cxx
    AliAnalysisTaskPythiaMpi.cxx
    AliAnalysisTaskpzpz.cxx
    AliAnalysisTaskpzpz.cxx
    AliAnalysisTaskR2P2.cxx
    AliAnalysisTaskSEpPbCorrelationsYS.cxx
    AliAnalysisTaskSEpPbCorrelationsMCYS.cxx
    AliAnalysisTaskSEpPbCorrelationsMCLEGOYS.cxx	
    AliTrigAssoPairST.cxx
    AliAnalysisTaskSEpPbCorrelationsJetV2.cxx
    AliAnalysisTaskSEPbPbCorrelationsJetV2.cxx
    AliAnaTaskMomScan.cxx
    AliAnalysisTaskTwoPlusOne.cxx
    AliDptDptInMC.cxx
    AliLeadingV0Correlation.cxx
    AliAnalysisTaskV0ChCorrelation.cxx
    AliAnalysisTaskV0ChCorrelationpp.cxx
    AliAnalysisTaskV0ChCorrelationppsys.cxx
    AliAnalysisTaskV0ChCorrelationsys.cxx
    AliPhiCorrelationsQATask.cxx
    AliAnalysisTaskDiJetCorrelationsAllb2b.cxx
    AliAnalysisTaskEPCorrPP.cxx
    AliAnalysisTaskEPCorrAA.cxx
    AliAnalysisTaskEPCorrAAMC.cxx
    AliAnalysisTaskPIDBFDptDpt.cxx
    AliAnalysisTaskGeneralBF.cxx
    AliCorrelationAnalysis.cxx
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
    Unfoldedhistos/AliAnalysisTaskCorrelationsStudies.cxx
    Unfoldedhistos/AliCSAnalysisCutsBase.cxx
    Unfoldedhistos/AliCSEventCuts.cxx
    Unfoldedhistos/AliCSPairAnalysis.cxx
    Unfoldedhistos/AliCSPIDCuts.cxx
    Unfoldedhistos/AliCSTrackCutsBase.cxx
    Unfoldedhistos/AliCSTrackCuts.cxx
    Unfoldedhistos/AliCSTrackMaps.cxx
    Unfoldedhistos/AliCSTrackSelection.cxx
    Unfoldedhistos/AliDptDptCorrelations.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES)
set(ALIROOT_DEPENDENCIES ANALYSISalice CORRFW PWGmuon PWGCFCorrelationsBase PWGJETFW OADB PWGDevNanoAOD EventMixing)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library. Note the extra includes ("Tracks UserTasks")
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}" "Tracks UserTasks")

# Create an object to be reused in case of static libraries
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)

# Add a library to the project using the object
add_library_tested(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)
target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

# Additional compilation and linking flags
set(MODULE_COMPILE_FLAGS " ${MODULE_COMPILE_FLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Setting compilation flags for the object
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
# Setting the linking flags for the library
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)

set(LIBDEPS )

# Macros
install(FILES MuonHadron/AddAnalysisTaskMuonHadronCorrelations.C DESTINATION PWGCF/Correlations/DPhi/MuonHadron/)
