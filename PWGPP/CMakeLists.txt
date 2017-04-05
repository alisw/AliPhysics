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

# Libraries from subdirectories
add_subdirectory(pid)
add_subdirectory(EVCHAR)
add_subdirectory(EVCHAR/FlowVectorCorrections/QnCorrections)
add_subdirectory(EVCHAR/FlowVectorCorrections/QnCorrectionsInterface)
add_subdirectory(MUON)
add_subdirectory(EMCAL)
add_subdirectory(Alignment)

#Module
set(MODULE PWGPP)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGPP)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/TENDER/Tender
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/PWG/CaloTrackCorrBase
                    ${AliPhysics_SOURCE_DIR}/OADB
  )

# Sources as introduced by the authors
#file ( GLOB SRCS1 "*.cxx" )
set ( SRCS1
  AliAnaFwdDetsQA.cxx
  AliAnalysisTaskFilteredTree.cxx
  AliAnalysisTaskIPInfo.cxx
  AliAnalysisTaskITSTPCalignment.cxx
  AliAnalysisTaskQASym.cxx
  AliAnalysisTaskV0QA.cxx
  AliAnalysisTaskVtXY.cxx
  AliAnaVZEROQA.cxx
  AliFilteredTreeAcceptanceCuts.cxx
  AliFilteredTreeEventCuts.cxx
  AliIntSpotEstimator.cxx
  AliRelAlignerKalmanArray.cxx
  AliTaskCDBconnect.cxx
  AliTrackComparison.cxx
  AliTrackComparisonESD.cxx
  AliOfflineTrigger.cxx
  )
#file ( GLOB SRCS2 "global/*.cxx" )
set ( SRCS2
  global/AliAnalysisTaskGlobalQA.cxx
  global/AliAnalysisTaskVertexESD.cxx
  global/AliAnalysisTaskCheckPileup.cxx
  global/AliAnalysisTaskCheckESDTracks.cxx
  global/AliAnalysisTaskCheckAODTracks.cxx
  global/AliGlobalFBFqa.cxx
  )
#file ( GLOB SRCS3 "cosmic/*.cxx" )
set ( SRCS3
  cosmic/AliAnalysisTaskCosmic.cxx
  )
#file ( GLOB SRCS4 "tender/*.cxx" )
set ( SRCS4
  tender/AliTenderSupplyTRD.cxx
  )
#file ( GLOB SRCS5 "HMPID/*.cxx" )
set ( SRCS5
  HMPID/AliHMPIDTaskQA.cxx
  )
#file ( GLOB SRCS6 "ITS/*.cxx" )
set ( SRCS6
  ITS/AliAlignmentDataFilterITS.cxx
  ITS/AliAnalysisTaskdEdxSSDQA.cxx
  ITS/AliAnalysisTaskITSAlignQA.cxx
  ITS/AliAnalysisTaskITSsaTracks.cxx
  ITS/AliAnalysisTaskITSTrackingCheck.cxx
  ITS/AliAnalysisTaskSDDRP.cxx
  ITS/AliAnalysisTaskSEImpParRes.cxx
  ITS/AliAnalysisTaskSEImpParResSparse.cxx
  ITS/AliAnalysisTaskSPD.cxx
  ITS/AliMeanVertexCalibTask.cxx
  ITS/AliMeanVertexPreprocessorOffline.cxx
  ITS/AliSPDUtils.cxx
  ITS/AliTrackMatchingTPCITSCosmics.cxx
  )
#file ( GLOB SRCS7 "TRD/*.cxx" )
set ( SRCS7
  TRD/AliTRDalignmentTask.cxx
  TRD/AliTRDcalibration.cxx
  TRD/AliTRDcheckDET.cxx
  TRD/AliTRDcheckESD.cxx
  TRD/AliTRDcheckPID.cxx
  TRD/AliTRDcheckTRK.cxx
  TRD/AliTRDclusterResolution.cxx
  TRD/AliTRDdigitsFilter.cxx
  TRD/AliTRDefficiency.cxx
  TRD/AliTRDefficiencyMC.cxx
  TRD/AliTRDinfoGen.cxx
  TRD/AliTRDmultiplicity.cxx
  TRD/AliTRDonlineTrackletFilter.cxx
  TRD/AliTRDonlineTrackletQA.cxx
  TRD/AliTRDpidRefMaker.cxx
  TRD/AliTRDpidRefMakerLQ.cxx
  TRD/AliTRDpidRefMakerNN.cxx
  TRD/AliTRDpwgppHelper.cxx
  TRD/AliTRDrecoTask.cxx
  TRD/AliTRDresolution.cxx
  TRD/AliTRDv0Monitor.cxx
  )
#file ( GLOB SRCS8 "TPC/*.cxx" )
set ( SRCS8
  TPC/AliAnalysisNoiseTPC.cxx
  TPC/AliComparisonDraw.cxx
  TPC/AliComparisonObject.cxx
  TPC/AliESDRecInfo.cxx
  TPC/AliESDRecKinkInfo.cxx
  TPC/AliESDRecV0Info.cxx
  TPC/AliESDresolMakerFast.cxx
  TPC/AliESDresolParams.cxx
  TPC/AliGenInfoMaker.cxx
  TPC/AliGenInfoTask.cxx
  TPC/AliGenKinkInfo.cxx
  TPC/AliGenV0Info.cxx
  TPC/AliMaterialBudget.cxx
  TPC/AliMCInfoCuts.cxx
  TPC/AliMCInfo.cxx
  TPC/AliMCTrackingTestTask.cxx
  TPC/AliPerfAnalyzeInvPt.cxx
  TPC/AliPerformanceDCA.cxx
  TPC/AliPerformanceDEdx.cxx
  TPC/AliPerformanceEff.cxx
  TPC/AliPerformanceMatch.cxx
  TPC/AliPerformanceMC.cxx
  TPC/AliPerformanceObject.cxx
  TPC/AliPerformancePtCalib.cxx
  TPC/AliPerformancePtCalibMC.cxx
  TPC/AliPerformanceRes.cxx
  TPC/AliPerformanceTask.cxx
  TPC/AliPerformanceTPC.cxx
  TPC/AliRecInfoCuts.cxx
  TPC/AliRecInfoMaker.cxx
  TPC/AliTaskConfigOCDB.cxx
  TPC/AliTPCComparisonPID.cxx
  TPC/AliTPCPerformanceSummary.cxx
  TPC/AliTPCtaskPID.cxx
  TPC/AliTPCtaskQA.cxx
  TPC/AliTreeDraw.cxx
  TPC/AliTPCPIDBase.cxx
  TPC/AliTPCPIDEtaTree.cxx
  TPC/AliTPCPIDEtaQA.cxx
  TPC/AliTPCcalibResidualPID.cxx
  )
#file ( GLOB SRCS9 "TOF/*.cxx" )
set ( SRCS9
  TOF/AliAnalysisTaskTOFqa.cxx
  TOF/AliAnalysisTaskTOFqaID.cxx
  TOF/AliTOFAnalysisTaskCalibPass0.cxx
  TOF/AliTOFAnalysisTaskCalibTree.cxx
  )
#file ( GLOB SRCS10 "TRD/info/*.cxx" )
set ( SRCS10
  TRD/info/AliTRDchmbInfo.cxx
  TRD/info/AliTRDclusterInfo.cxx
  TRD/info/AliTRDeventCuts.cxx
  TRD/info/AliTRDeventInfo.cxx
  TRD/info/AliTRDpidInfo.cxx
  TRD/info/AliTRDtrackInfo.cxx
  TRD/info/AliTRDtrendingManager.cxx
  TRD/info/AliTRDtrendValue.cxx
  TRD/info/AliTRDtriggerInfo.cxx
  TRD/info/AliTRDv0Info.cxx
  )
#file ( GLOB SRCS11 "background/*.cxx" )
set ( SRCS11
  background/AliHistoListWrapper.cxx
  )
#file ( GLOB SRCS12 "Centrality/*.cxx" )
set ( SRCS12
  Centrality/AliAnalysisTaskHIMultCorr.cxx
  Centrality/AliMultiplicityCorrelations.cxx
  )
#file ( GLOB SRCS13 "ZDC/*.cxx" )
set ( SRCS13
  ZDC/AliAnalysisTaskZDCpp.cxx
  ZDC/AliAnalysisTaskZDCpA.cxx
  ZDC/AliAnalysisTaskZDCPbPb.cxx
  )
#file ( GLOB SRCS14 "T0/*.cxx" )
set ( SRCS14
  T0/AliT0AnalysisTaskQA.cxx
  )
#file ( GLOB SRCS15 "VZERO/*PbPb.cxx" )
set (SRCS15
  VZERO/AliAnaVZEROPbPb.cxx
  )
#file ( GLOB SRCS16 "EvTrkSelection/*.cxx" )
set ( SRCS16
  EvTrkSelection/AliAnalysisTrackingUncertainties.cxx
  EvTrkSelection/AliAnalysisTrackingUncertaintiesAOT.cxx
  EvTrkSelection/AliCFSingleTrackEfficiencyTask.cxx
  EvTrkSelection/AliSingleTrackEffCuts.cxx
  EvTrkSelection/AliAnalysisTaskEventCutsValidation.cxx
  )
#file ( GLOB SRCS17 "TRD/TRDPID/*.cxx" )
set ( SRCS17
  TRD/TRDPID/AliTRDPIDTree.cxx
  TRD/TRDPID/AliTRDPIDmatching.cxx
  )
#file ( GLOB SRCS18 "AD/*.cxx" )
set ( SRCS18
  AD/AliAnalysisTaskADQA.cxx
  AD/AliAnalysisTaskADVVQA.cxx
  AD/AliAnalysisTaskADPilot.cxx
  AD/ADESDFriendUtils.cxx
  AD/AliAnalysisTaskADCalib.cxx
  AD/AliAnalysisTaskADChargeMonitoring.cxx
  )
#file ( GLOB SRCS19 "BeamGasMonitoring/*.cxx" )
set ( SRCS19
#  BeamGasMonitoring/AliAnalysisMBVeto.cxx temporarily commented:sjena
  BeamGasMonitoring/AliAnalysisBGMonitorQAHLT.cxx
  BeamGasMonitoring/AliAnalysisBGMonitorQATrain.cxx
  )

#file ( GLOB SRCS19 "BeamGasMonitoring/*.cxx" )
set ( SRCS20
  HMTFqa/AliAnalysisMultCorrTaskQA.cxx
  HMTFqa/AliAnalysisTaskLeadingPt.cxx
  HMTFqa/AliAnalysisMultCorrTaskQAPF.cxx
  )

list ( APPEND SRCS
  ${SRCS1}
  ${SRCS2}
  ${SRCS3}
  ${SRCS4}
  ${SRCS5}
  ${SRCS6}
  ${SRCS7}
  ${SRCS8}
  ${SRCS9}
  ${SRCS10}
  ${SRCS11}
  ${SRCS12}
  ${SRCS13}
  ${SRCS14}
  ${SRCS15}
  ${SRCS16}
  ${SRCS17}
  ${SRCS18}
  ${SRCS19}
  ${SRCS20}
  )

###add_subdirectory(Alignment)

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
set(LIBDEPS ANALYSISalice OADB CDB CORRFW EMCALUtils ITSrec HMPIDbase PWGCaloTrackCorrBase STEER Tender TOFrec TPCrec TRDrec Postscript ADbase ADrec)
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

# install macros as needed
install(DIRECTORY macros DESTINATION PWGPP)
install(DIRECTORY CalibMacros DESTINATION PWGPP)
install(DIRECTORY PilotTrain DESTINATION PWGPP)
install(DIRECTORY CalibMacros DESTINATION PWGPP)
install(DIRECTORY EMCAL/QA/macros DESTINATION PWGPP/EMCAL/QA USE_SOURCE_PERMISSIONS)
install(DIRECTORY EMCAL/macros DESTINATION PWGPP/EMCAL)

install(DIRECTORY HMTFqa/macros DESTINATION PWGPP/HMTFqa)
install(FILES HMTFqa/AddTaskMultCorrTaskQAPF.C DESTINATION PWGPP/HMTFqa)

install(FILES HMPID/AddTaskHmpidQA.C DESTINATION PWGPP/HMPID)
install(FILES HMPID/macros/copyHMPQAfiles.C
 	      HMPID/macros/makeHMPQA.C DESTINATION PWGPP/HMPID/macros)
install(FILES T0/AddTaskT0QA.C
	      T0/MakeTrendT0.C
	      T0/drawPerformanceT0QATrends.C
	      DESTINATION PWGPP/T0)
install(FILES ITS/QA/MakeTrendingITSQA.C
	      ITS/QA/DrawTrendingITSQA.C
	      DESTINATION PWGPP/ITS/QA)
install(FILES TOF/AddTaskTOFQA.C
              TOF/AddTaskTOFqaID.C
              DESTINATION PWGPP/TOF)
install(FILES TOF/trending/MakeTrendingTOFQA.C
	      TOF/trending/DrawTrendingTOFQA.C
	      TOF/trending/MakeTrendingTOFQAv2.C
              DESTINATION PWGPP/TOF/trending)
install(FILES TOF/AddTaskTOFQA.C DESTINATION PWGPP/TOF)
install(DIRECTORY TPC/macros DESTINATION PWGPP/TPC PATTERN ".C")
install(DIRECTORY TPC/CalibMacros DESTINATION PWGPP/TPC PATTERN ".C")
install(DIRECTORY TPC/CalibMacros USE_SOURCE_PERMISSIONS DESTINATION PWGPP/TPC PATTERN ".sh")
install(DIRECTORY rawmerge  USE_SOURCE_PERMISSIONS DESTINATION PWGPP PATTERN ".sh")
install(DIRECTORY TRD/macros DESTINATION PWGPP/TRD FILES_MATCHING PATTERN "*.C")
install(FILES TRD/run.C TRD/runQAesds.C TRD/cacheOCDB.C DESTINATION PWGPP/TRD)
install(FILES TRD/TRDPID/macros/AddTaskTRDPIDTree.C
              TRD/TRDPID/macros/AddTaskTRDPIDmatching.C
DESTINATION PWGPP/TRD/TRDPID/macros)
install(DIRECTORY TRD/data DESTINATION PWGPP/TRD FILES_MATCHING PATTERN "*.root")
install(FILES ZDC/AddTaskZDCQA.C
              ZDC/AddTaskZDCQApp.C
              ZDC/AddTaskZDCQAPbPb.C
              ZDC/AddTaskZDCQApA.C
              DESTINATION PWGPP/ZDC)
install(FILES ZDC/macros/DrawPerformanceZDCQAMatch.C
	      ZDC/macros/DrawPerformanceZDCQAMatchTrends.C
	      ZDC/macros/MakeTrendZDC.C
	      DESTINATION PWGPP/ZDC/macros)
install(FILES AD/AddTaskADQA.C
	      AD/AddTaskADVVQA.C
	      AD/AddTaskADPilot.C
	      AD/AddTaskADCalib.C
	      AD/AddTaskADChargeMonitoring.C
	      DESTINATION PWGPP/AD)
install(FILES AD/trending/DrawTrendingADQA.C
 	      AD/trending/MakeTrendingADQA.C
	      DESTINATION PWGPP/AD/trending)
install(FILES PilotTrain/AddTaskMuonQA.C DESTINATION PWGPP/PilotTrain/ )
install(FILES VZERO/AddTaskVZEROPbPb.C DESTINATION PWGPP/VZERO/)
install(FILES VZERO/trending/DrawTrendingV0QA.C
 	      VZERO/trending/MakeTrendingV0QA.C
	      DESTINATION PWGPP/VZERO/trending)
install(DIRECTORY MUON/lite DESTINATION PWGPP/MUON FILES_MATCHING PATTERN "*.C")
install(DIRECTORY MUON/dep DESTINATION PWGPP/MUON FILES_MATCHING PATTERN "*.C")
install(DIRECTORY MUON/lite DESTINATION PWGPP/MUON FILES_MATCHING PATTERN "*.config")
install(DIRECTORY EVS DESTINATION PWGPP)
install(DIRECTORY scripts DESTINATION PWGPP USE_SOURCE_PERMISSIONS)
install(DIRECTORY QA/scripts DESTINATION PWGPP/QA USE_SOURCE_PERMISSIONS)
install(DIRECTORY QA/detectorQAscripts DESTINATION PWGPP/QA USE_SOURCE_PERMISSIONS)
install(DIRECTORY QA/Tracking  DESTINATION PWGPP/QA USE_SOURCE_PERMISSIONS)
install(DIRECTORY QA/Tracking/ExpertQA DESTINATION PWGPP/QA/Tracking USE_SOURCE_PERMISSIONS)
install(DIRECTORY benchmark DESTINATION PWGPP USE_SOURCE_PERMISSIONS)
install(DIRECTORY EvTrkSelection/macros DESTINATION PWGPP/EvTrkSelection FILES_MATCHING PATTERN "*.C")
install(FILES EvTrkSelection/macros/AddTaskTrackingUncert.C
	      EvTrkSelection/macros/AddSingleTrackEfficiencyTaskForAutomaticQA.C
	      DESTINATION PWGPP/EvTrkSelection)
install(FILES CalibMacros/CPass0/LoadLibraries.C
	      CalibMacros/CPass0/ConfigCalibTrain.C
	      CalibMacros/CPass0/AddTaskTPCCalib.C
	      CalibMacros/CPass0/AddTOFAnalysisTaskCalibPass0.C
	      CalibMacros/CPass0/AddTaskT0Calib.C
	      CalibMacros/CPass0/AddTaskMeanVertexCalib.C
	      CalibMacros/CPass0/AddTaskSDDCalib.C
	      DESTINATION PWGPP/CalibMacros/CPass0)
install(FILES CalibMacros/CPass1/LoadLibraries.C
	      CalibMacros/CPass1/ConfigCalibTrain.C
	      CalibMacros/CPass1/AddTaskTPCCalib.C
	      CalibMacros/CPass1/AddTOFAnalysisTaskCalibPass0.C
	      CalibMacros/CPass1/AddTOFAnalysisTaskCalibTree.C
	      CalibMacros/CPass1/AddTaskT0Calib.C
	      CalibMacros/CPass1/AddTaskT0Analysis.C
	      CalibMacros/CPass1/AddTaskSDDCalib.C
	      DESTINATION PWGPP/CalibMacros/CPass1)
install(FILES CalibMacros/commonMacros/CleanGeom.C
	      DESTINATION PWGPP/CalibMacros/commonMacros)
install(FILES BeamGasMonitoring/macros/AddTaskMBVeto.C
              BeamGasMonitoring/macros/AddTaskBGMonitorQAHLT.C
	      BeamGasMonitoring/macros/AddTaskBGMonitorQATrain.C
	      DESTINATION PWGPP/BeamGasMonitoring/macros)

message(STATUS "PWGPP enabled")
