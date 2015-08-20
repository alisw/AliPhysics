AliAnalysisTaskCheckPileup *AddTaskPileupChecks(TString suffix="",
						Int_t nContribSPD=3,
						Double_t zDiffSPD=0.8,
						Bool_t doNtuple=kTRUE,
						Int_t nContribMV=5,
						Double_t zDiffMV=15.,
						Double_t chi2MV=5.,
						Bool_t flagBCMV=kFALSE,
						Bool_t readMC=kFALSE
						){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPileup", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPileup", "This task requires an input event handler");
    return NULL;
  }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AddTaskPileup", "This task requires to run on ESD");
    return NULL;
  }
  
  //Bool_t isMC=kFALSE;
  //if (mgr->GetMCtruthEventHandler()) isMC=kTRUE;
  
  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }
  // Create and configure the task
  printf("CREATE PILEUP TASK\n");

  AliAnalysisTaskCheckPileup *taskpil = new AliAnalysisTaskCheckPileup();
  taskpil->SetFillNtuple(doNtuple);
  taskpil->SetCutOnContribToSPDPileupVert(nContribSPD);
  taskpil->SetCutOnSPDZDiff(zDiffSPD);
  taskpil->ConfigureMultiTrackVertexPileup(nContribMV,zDiffMV,chi2MV,flagBCMV);
  mgr->AddTask(taskpil);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += Form(":CheckPileup%s",suffix.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("clistPrimaryV%s",suffix.Data()),
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName );
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("clistPileupSPD%s",suffix.Data()),
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName );
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("clistPileupMV%s",suffix.Data()),
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName );
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(Form("cCounters%s",suffix.Data()),
							    AliCounterCollection::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName);

  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(Form("cNtuple%s",suffix.Data()),
							    TNtuple::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName);

  mgr->ConnectInput(taskpil, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskpil, 1, coutput1);
  mgr->ConnectOutput(taskpil, 2, coutput2);
  mgr->ConnectOutput(taskpil, 3, coutput3);
  mgr->ConnectOutput(taskpil, 4, coutput4);
  mgr->ConnectOutput(taskpil, 5, coutput5);
  return taskpil;
}   
