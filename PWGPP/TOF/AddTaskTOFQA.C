/*  created by fbellini@cern.ch on 14/09/2010 */
/*  last modified by fbellini   on 11/11/2011 */

AliAnalysisTaskSE * AddTaskTOFQA(Bool_t flagEnableAdvancedCheck=kFALSE) 
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
  //AliLog::SetClassDebugLevel("AliAnalysisTaskTOFqa",1);
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
  // TPC  
  esdTrackCutsStd2010->SetMinNClustersTPC(70); 
  esdTrackCutsStd2010->SetMaxChi2PerClusterTPC(4);
  esdTrackCutsStd2010->SetAcceptKinkDaughters(kFALSE); 
  esdTrackCutsStd2010->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCutsStd2010->SetRequireITSRefit(kTRUE);
  esdTrackCutsStd2010->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  esdTrackCutsStd2010->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//selects primaries
  esdTrackCutsStd2010->SetMaxDCAToVertexZ(2);
  esdTrackCutsStd2010->SetDCAToVertex2D(kFALSE);
  esdTrackCutsStd2010->SetRequireSigmaToVertex(kFALSE);

  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  trackFilter->AddCuts(esdTrackCutsStd2010);
  task->SetTrackFilter(trackFilter);
   
  
  // Create containers for input/output
  AliAnalysisDataContainer *cInputTOFqa = mgr->CreateContainer("cInputTOFqa",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *cGeneralTOFqa = mgr->CreateContainer("cGeneralTOFqa",TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));
  AliAnalysisDataContainer *cTimeZeroTOFqa = mgr->CreateContainer("cTimeZeroTOFqa",TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));
   AliAnalysisDataContainer *cPIDTOFqa = mgr->CreateContainer("cPIDTOFqa",TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));
   AliAnalysisDataContainer *cPosTracksTOFqa = mgr->CreateContainer("cPosTracksTOFqa",TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));
   AliAnalysisDataContainer *cNegTracksTOFqa = mgr->CreateContainer("cNegTracksTOFqa",TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));

  // Attach i/o
  mgr->ConnectInput(task, 0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cGeneralTOFqa);
  mgr->ConnectOutput(task, 2, cTimeZeroTOFqa);
  mgr->ConnectOutput(task, 3, cPIDTOFqa);
  mgr->ConnectOutput(task, 4, cPosTracksTOFqa);
  mgr->ConnectOutput(task, 5, cNegTracksTOFqa);
  return task;
}
