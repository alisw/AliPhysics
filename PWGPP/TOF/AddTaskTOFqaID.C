/*  created by fbellini@cern.ch on 29/04/2013 */
/*  last modified by fbellini   on 29/04/2013 */

// UInt_t kTriggerMuonAll = AliVEvent::kMUL7 | AliVEvent::kMUSH7 | AliVEvent::kMUU7 | AliVEvent::kMUS7;
// UInt_t kTriggerMuonBarell = AliVEvent::kMUU7;
// UInt_t kTriggerEMC   = AliVEvent::kEMC7;
// UInt_t kTriggerHM   = AliVEvent::kHighMult;
UInt_t kTriggerInt = AliVEvent::kAnyINT;
UInt_t kTriggerMask = kTriggerInt;

AliAnalysisTaskSE * AddTaskTOFqaID(UInt_t triggerMask = kTriggerMask, 
				   Bool_t flagEnableAdvancedCheck=kFALSE, 
				   Bool_t useStdCuts2011 = kTRUE, 
				   Bool_t isMC = kFALSE, 
				   Short_t absPdgCode = 0) 
{
  // Task for checking TOF QA
 
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }   

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 
  if (isMC && (absPdgCode<0)) {
    ::Error("AddTask","Invalid selected PDG code for MC. Set absPdgCode=0 for no specie selection.");
  }
 
  // Create the task
  AliAnalysisTaskTOFqaID *task = new AliAnalysisTaskTOFqaID(Form("taskTOFqaID_%i",absPdgCode));
  task->EnableAdvancedCheck(flagEnableAdvancedCheck);
  task->SetSelectMCspecies(isMC, absPdgCode);
  task->SelectCollisionCandidates(triggerMask);
  //AliLog::SetClassDebugLevel("AliAnalysisTaskTOFqaID",4);
  mgr->AddTask(task);

  /* cuts used for QA in 2010 p-p */
  /*
  AliESDtrackCuts* esdTrackCutsLoose2010 = new AliESDtrackCuts("AliESDtrackCuts", "esdTrackCutsLoose2010");
  esdTrackCutsLoose2010->SetMinNClustersTPC(70); 
  esdTrackCutsLoose2010->SetMaxChi2PerClusterTPC(3.5); 
  esdTrackCutsLoose2010->SetMaxCovDiagonalElements(2, 2, 0.5, 0.5, 2);
  esdTrackCutsLoose2010->SetRequireTPCRefit(kTRUE);
  esdTrackCutsLoose2010->SetMaxDCAToVertexXY(3.0); 
  esdTrackCutsLoose2010->SetMaxDCAToVertexZ(3.0); 
  esdTrackCutsLoose2010->SetRequireSigmaToVertex(kTRUE); 
  esdTrackCutsLoose2010->SetAcceptKinkDaughters(kFALSE); 
  esdTrackCutsLoose2010->SetMaxNsigmaToVertex(4.0);
  */

  /* standard cuts ITS-TPC 2010 */
  AliESDtrackCuts* esdTrackCutsStd2010 = new AliESDtrackCuts("AliESDtrackCuts", "Standard2010");
  esdTrackCutsStd2010->GetStandardITSTPCTrackCuts2010(kTRUE,0);
  // TPC  
  // esdTrackCutsStd2010->SetMinNClustersTPC(70); 
  // esdTrackCutsStd2010->SetMaxChi2PerClusterTPC(4);
  // esdTrackCutsStd2010->SetAcceptKinkDaughters(kFALSE); 
  // esdTrackCutsStd2010->SetRequireTPCRefit(kTRUE);
  // // ITS
  // esdTrackCutsStd2010->SetRequireITSRefit(kTRUE);
  // esdTrackCutsStd2010->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
  // 						AliESDtrackCuts::kAny);
  // esdTrackCutsStd2010->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//selects primaries
  // esdTrackCutsStd2010->SetMaxDCAToVertexZ(2);
  // esdTrackCutsStd2010->SetDCAToVertex2D(kFALSE);
  // esdTrackCutsStd2010->SetRequireSigmaToVertex(kFALSE);

 
  /* standard cuts ITS-TPC 2011 */
  AliESDtrackCuts* esdTrackCutsStd2011 = new AliESDtrackCuts("AliESDtrackCuts", "Standard2011");
  esdTrackCutsStd2011->GetStandardITSTPCTrackCuts2011(kTRUE,0);
  
  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  if (useStdCuts2011) 
    trackFilter->AddCuts(esdTrackCutsStd2011);
  else
    trackFilter->AddCuts(esdTrackCutsStd2010);
  task->SetTrackFilter(trackFilter);
  
  TString partName(task->GetSpeciesName(absPdgCode));
  
  // Create containers for input/output
  AliAnalysisDataContainer *cInputTOFqa = mgr->CreateContainer("cInputTOFqa",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *cGeneralTOFqa = mgr->CreateContainer(Form("base_%s",partName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName(), partName.Data()));
  AliAnalysisDataContainer *cTimeZeroTOFqa = mgr->CreateContainer(Form("timeZero_%s",partName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName(),partName.Data()));
   AliAnalysisDataContainer *cPIDTOFqa = mgr->CreateContainer(Form("pid_%s",partName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName(),partName.Data()));
   AliAnalysisDataContainer *cTRDcheckTOFqa = mgr->CreateContainer(Form("trd_%s",partName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName(),partName.Data()));
   AliAnalysisDataContainer *cTriggerTOFqa = mgr->CreateContainer(Form("trigger_%s",partName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName(),partName.Data()));

  // Attach i/o
  mgr->ConnectInput(task, 0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cGeneralTOFqa);
  mgr->ConnectOutput(task, 2, cTimeZeroTOFqa);
  mgr->ConnectOutput(task, 3, cPIDTOFqa);
  mgr->ConnectOutput(task, 4, cTRDcheckTOFqa);
  mgr->ConnectOutput(task, 5, cTriggerTOFqa);
  return task;
}
