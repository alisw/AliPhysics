/*  created by fbellini@cern.ch on 14/09/2010 */
/*  last modified by fbellini   on 11/11/2011 */

AliAnalysisTaskSE * AddTaskTOFQA(Bool_t flagEnableAdvancedCheck=kFALSE, 
				 UInt_t triggerMask = AliVEvent::kAnyINT,
				 Int_t trackCutSetTOFqa = 0,
				 TString cutName = "") 
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
  
  // Create the task
  AliAnalysisTaskTOFqa *task = new AliAnalysisTaskTOFqa("taskTOFqa");
  task->EnableAdvancedCheck(flagEnableAdvancedCheck);
  task->SelectCollisionCandidates(triggerMask);
  //AliLog::SetClassDebugLevel("AliAnalysisTaskTOFqa",1);
  mgr->AddTask(task);

  /* cuts used for QA in 2010 p-p */
    //Define track cut set
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "TrackCutsTOFqa"); 
  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");

  if ( (trackCutSetTOFqa<0) || (trackCutSetTOFqa>=AliAnalysisTaskTOFqaID::kNCutSetTOFqa) ) trackCutSetTOFqa = 0;

  if ( trackCutSetTOFqa == AliAnalysisTaskTOFqaID::kRun1Cuts ) {
    //use track cuts used for QA during run1 (before July 2014)
    esdTrackCuts->SetMinNClustersTPC(70); 
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE); 
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					   AliESDtrackCuts::kAny);
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//selects primaries
    esdTrackCuts->SetMaxDCAToVertexZ(2);
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  } else {
    if ( trackCutSetTOFqa == AliAnalysisTaskTOFqaID::kStd2010 ) esdTrackCuts->GetStandardITSTPCTrackCuts2010(kTRUE,0);
    if ( trackCutSetTOFqa == AliAnalysisTaskTOFqaID::kStd2010crossedRows ) esdTrackCuts->GetStandardITSTPCTrackCuts2010(kTRUE,1);
    if ( trackCutSetTOFqa == AliAnalysisTaskTOFqaID::kStd2011) esdTrackCuts->GetStandardITSTPCTrackCuts2011(kTRUE,0); 
    if ( trackCutSetTOFqa == AliAnalysisTaskTOFqaID::kStd2011crossedRows ) esdTrackCuts->GetStandardITSTPCTrackCuts2011(kTRUE,1);
  }
  
  trackFilter->AddCuts(esdTrackCuts);
  task->SetTrackFilter(trackFilter); 
  
  // Create containers for input/output
  AliAnalysisDataContainer *cInputTOFqa = mgr->CreateContainer("cInputTOFqa",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *cGeneralTOFqa = mgr->CreateContainer(Form("cGeneralTOFqa%s",cutName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));
  AliAnalysisDataContainer *cTimeZeroTOFqa = mgr->CreateContainer(Form("cTimeZeroTOFqa%s",cutName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));
  AliAnalysisDataContainer *cPIDTOFqa = mgr->CreateContainer(Form("cPIDTOFqa%s",cutName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));
  AliAnalysisDataContainer *cPosTracksTOFqa = mgr->CreateContainer(Form("cPosTracksTOFqa%s",cutName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));
  AliAnalysisDataContainer *cNegTracksTOFqa = mgr->CreateContainer(Form("cNegTracksTOFqa%s",cutName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));

  // Attach i/o
  mgr->ConnectInput(task, 0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cGeneralTOFqa);
  mgr->ConnectOutput(task, 2, cTimeZeroTOFqa);
  mgr->ConnectOutput(task, 3, cPIDTOFqa);
  mgr->ConnectOutput(task, 4, cPosTracksTOFqa);
  mgr->ConnectOutput(task, 5, cNegTracksTOFqa);
  return task;
}
