#ifdef __CINT__

#pragma link off all glols;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliTenderSupplyTRD+;

#pragma link C++ class AliAnaFwdDetsQA+;
#pragma link C++ class AliAnalysisTaskVtXY+;
#pragma link C++ class AliESDresolParams+;
#pragma link C++ class AliESDresolMakerFast+;
#pragma link C++ class AliOfflineTrigger+;


#pragma link C++ class AliAnalysisTaskGlobalQA+;

#pragma link C++ class AliTreeDraw+;
//
#pragma link C++ class AliTPCdigitRow+;
#pragma link C++ class AliTPCComparisonPID+;
#pragma link C++ class AliMCInfo+;
#pragma link C++ class AliGenV0Info+;
#pragma link C++ class AliGenKinkInfo+;
#pragma link C++ class AliGenInfoMaker+;
//
#pragma link C++ class AliESDRecInfo+;
#pragma link C++ class AliESDRecV0Info+;
#pragma link C++ class AliESDRecKinkInfo+;
//
#pragma link C++ class AliRecInfoMaker+;
#pragma link C++ class AliComparisonDraw+;

#pragma link C++ class AliRecInfoCuts+;
#pragma link C++ class AliMCInfoCuts+;
#pragma link C++ class AliComparisonObject+;

#pragma link C++ class AliGenInfoTask+;
#pragma link C++ class AliMCTrackingTestTask+;
#pragma link C++ class AliTPCtaskPID+;
#pragma link C++ class AliTPCtaskQA+;

#pragma link C++ class AliPerformanceTask+;
#pragma link C++ class AliPerformanceObject+;
#pragma link C++ class AliPerformanceRes+;
#pragma link C++ class AliPerformanceEff+;
#pragma link C++ class AliPerformanceDEdx+;
#pragma link C++ class AliPerformanceDCA+;
#pragma link C++ class AliPerformanceTPC+;
#pragma link C++ class AliPerformanceMC+;
#pragma link C++ class AliPerformanceMatch+;
#pragma link C++ class AliPerformancePtCalib+;
#pragma link C++ class AliPerformancePtCalibMC+;
#pragma link C++ class AliPerfAnalyzeInvPt+;
#pragma link C++ class AliTPCPerformanceSummary+;
#pragma link C++ class AliAnalysisNoiseTPC+;

#pragma link C++ class AliIntSpotEstimator+;
#pragma link C++ class AliAnalysisTaskIPInfo+;

#pragma link C++ class AliAnalysisTaskVertexESD+;
#pragma link C++ class AliAnalysisTaskCheckPileup+;
#pragma link C++ class AliAnalysisTaskCheckESDTracks+;
#pragma link C++ class AliAnalysisTaskCheckAODTracks+;
#pragma link C++ class AliAlignmentDataFilterITS+;
#pragma link C++ class AliAnalysisTaskITSTrackingCheck+;
#pragma link C++ class AliAnalysisTaskITSsaTracks+;
#pragma link C++ class AliAnalysisTaskITSAlignQA+;
#pragma link C++ class AliAnalysisTaskSEImpParRes+;
#pragma link C++ class AliAnalysisTaskSEImpParResSparse+;
#pragma link C++ class AliTrackMatchingTPCITSCosmics+;
#pragma link C++ class AliAnalysisTaskV0QA+;
#pragma link C++ class AliMaterialBudget+;
#pragma link C++ class AliAnalysisTaskSPD+;
#pragma link C++ class AliAnalysisTaskSDDRP+;
#pragma link C++ class AliSPDUtils+;
#pragma link C++ class AliAnalysisTaskdEdxSSDQA+;
#pragma link C++ class AliMeanVertexCalibTask+;
#pragma link C++ class AliMeanVertexPreprocessorOffline+;

#pragma link C++ class AliRelAlignerKalmanArray+;
#pragma link C++ class AliAnalysisTaskITSTPCalignment+;

#pragma link C++ class AliAnalysisTaskQASym+;
#pragma link C++ class AliAnaVZEROQA+;
#pragma link C++ class AliAnaVZEROPbPb+;

#pragma link C++ class AliAnalysisTaskFilteredTree+;
#pragma link C++ class AliFilteredTreeEventCuts+;
#pragma link C++ class AliFilteredTreeAcceptanceCuts+;

#pragma link C++ class AliTaskConfigOCDB+;


// TRD performance classes
#pragma link C++ class  AliTenderSupplyTRD+;
#pragma link C++ class  AliTRDpwgppHelper+;
#pragma link C++ class  AliTRDtrendValue+;
#pragma link C++ class  AliTRDtrendingManager+;
#pragma link C++ class  AliTRDclusterInfo+;
#pragma link C++ class  AliTRDtriggerInfo+;
#pragma link C++ class  AliTRDchmbInfo+;
#pragma link C++ class  AliTRDv0Info+;
#pragma link C++ class  AliTRDtrackInfo+;
#pragma link C++ class  AliTRDtrackInfo::AliESDinfo+;
#pragma link C++ class  AliTRDtrackInfo::AliMCinfo+;
#pragma link C++ class  AliTRDeventCuts+;
#pragma link C++ class  AliTRDeventInfo+;
#pragma link C++ class  AliTRDpidInfo+;
#pragma link C++ class  AliTRDpidInfo::AliTRDpidData+;
#pragma link C++ class  AliTRDinfoGen+;
#pragma link C++ class  AliTRDrecoTask+;
#pragma link C++ class  AliTRDrecoTask::AliTRDrecoProjection+;
#pragma link C++ class  AliTRDcheckESD+;
#pragma link C++ class  AliTRDcheckDET+;
#pragma link C++ class  AliTRDcheckPID+;
#pragma link C++ class  AliTRDcheckTRK+;
#pragma link C++ class  AliTRDresolution+;
#pragma link C++ class  AliTRDefficiency+;
#pragma link C++ class  AliTRDefficiencyMC+;
#pragma link C++ class  AliTRDv0Monitor+;
#pragma link C++ class  AliTRDonlineTrackletFilter;
#pragma link C++ class  AliTRDonlineTrackletQA;
#pragma link C++ class  AliTRDPIDTree+;
#pragma link C++ class  AliTRDPIDmatching+;
#pragma link C++ class  AliTRDdigitsFilter;

// TRD offline calibration classes
#pragma link C++ class  AliTRDmultiplicity+;
#pragma link C++ class  AliTRDclusterResolution+;
#pragma link C++ class  AliTRDalignmentTask+;
#pragma link C++ class  AliTRDcalibration+;
#pragma link C++ class  AliTRDpidRefMaker+;
#pragma link C++ class  AliTRDpidRefMakerLQ+;
#pragma link C++ class  AliTRDpidRefMakerNN+;
//TPC PID calibration and QA
#pragma link C++ class  AliTPCPIDBase+;
#pragma link C++ class  AliTPCPIDEtaTree+;
#pragma link C++ class  AliTPCPIDEtaQA+;
#pragma link C++ class  AliTPCcalibResidualPID+;
// TOF QA
#pragma link C++ class  AliAnalysisTaskTOFqa+;
#pragma link C++ class  AliAnalysisTaskTOFqaID+;
#pragma link C++ class  AliTOFAnalysisTaskCalibPass0+;
#pragma link C++ class  AliTOFAnalysisTaskCalibTree+;
// HMPID QA
#pragma link C++ class  AliHMPIDTaskQA+;
// Cosmics QA
#pragma link C++ class  AliAnalysisTaskCosmic+;
// Background and luminosity studies. The obsolete task was removed. The list wrapper is still used in other classes.
#pragma link C++ class  AliHistoListWrapper+;
// CDB connect
#pragma link C++ class  AliTaskCDBconnect+;
// Centrality classes
#pragma link C++ class  AliMultiplicityCorrelations+;
#pragma link C++ class  AliAnalysisTaskHIMultCorr+;
// ZDC
#pragma link C++ class  AliAnalysisTaskZDCpp+;
#pragma link C++ class  AliAnalysisTaskZDCpA+;
#pragma link C++ class  AliAnalysisTaskZDCPbPb+;
// T0
#pragma link C++ class  AliT0AnalysisTaskQA+;
// AD
#pragma link C++ class  AliAnalysisTaskADQA+;
#pragma link C++ class  AliAnalysisTaskADVVQA+;
#pragma link C++ class  AliAnalysisTaskADPilot+;
#pragma link C++ class  ADESDFriendUtils+;
#pragma link C++ class  AliAnalysisTaskADCalib+;
#pragma link C++ class  AliAnalysisTaskADChargeMonitoring+;
//Beam Gas Monitoring HMTF
//#pragma link C++ class  AliAnalysisMBVeto+;
#pragma link C++ class  AliAnalysisMultCorrTaskQA+;
#pragma link C++ class  AliAnalysisTaskLeadingPt+;
#pragma link C++ class  AliAnalysisBGMonitorQATrain+;
#pragma link C++ class  AliAnalysisBGMonitorQAHLT+;
#pragma link C++ class  AliAnalysisMultCorrTaskQAPF+;

#pragma link C++ class  AliTrackComparison+;
#pragma link C++ class  AliTrackComparisonESD+;
#pragma link C++ class  AliAnalysisTaskGlobalQA+;
#pragma link C++ class  AliGlobalFBFqa+;

//EVENT AND TRACK SELECTION
#pragma link C++ class  AliCFSingleTrackEfficiencyTask+;
#pragma link C++ class  AliSingleTrackEffCuts+;
#pragma link C++ class  AliAnalysisTrackingUncertaintiesAOT+;
#pragma link C++ class  AliAnalysisTrackingUncertainties+;
#pragma link C++ class  AliAnalysisTaskEventCutsValidation+;

// vdM task
#pragma link C++ class  AliAnalysisTaskVdM::TreeData+;
#pragma link C++ class  AliAnalysisTaskVdM::EventInfo+;
#pragma link C++ class  AliAnalysisTaskVdM::ADV0+;
#pragma link C++ class  AliAnalysisTaskVdM+;
#pragma link C++ class  AliXMLEngine::NodeIterator+;
#pragma link C++ class  AliXMLEngine::Node+;
#pragma link C++ class  AliXMLEngine::Attr+;
#pragma link C++ class  AliXMLEngine+;

#endif
