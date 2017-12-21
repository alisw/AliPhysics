
AliAnalysisTaskMLTreeMaker *AddTaskMLTreeMaker(TString taskname = "ESDExample", 
                                             Double_t etaMin = -0.8,
                                             Double_t etaMax = 0.8,
                                             Double_t ptMin = 0.2,
                                             Double_t ptMax = 10.0
					     ) {				


AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
if (!mgr) {
  ::Error("AddTaskBalancePsiCentralityTrain",  "No analysis manager to connect to.");
  return NULL;
}

// Check the analysis type using the event handlers connected to the analysis manager.
//===========================================================================
if (!mgr->GetInputEventHandler()) {
  ::Error("AddTaskMLTreeMaker",  "This task requires an input event handler");
  return NULL;
}
TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

 if (analysisType!="ESD"){
   ::Warning("AddTaskMLTreeMaker",  "analysis type NOT ESD --> some variables might not be filled");
 }
      

AliAnalysisTaskMLTreeMaker *taskESD = new AliAnalysisTaskMLTreeMaker(taskname);

   // ==========================================================================
  // user customization part

//  taskESD->SetEtaRange(etaMin, etaMax);
//  taskESD->SetPtRange(ptMin, ptMax);
  taskESD->SelectCollisionCandidates(AliVEvent::kINT7);
//  taskESD->SetLoCuts(kTRUE);
//  taskESD->SetFilterBit(AliDielectronTrackCuts::kTPCqualSPDany);// for AOD analyses: TPC cuts + any SPD hit
  // ==========================================================================

  mgr->AddTask(taskESD);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  
// AliAnalysisDataContainer *coutESD = mgr->CreateContainer("",  TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
// AliAnalysisDataContainer *coutESD1 = mgr->CreateContainer("QAHist",  TH1::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
// mgr->ConnectInput(taskESD,  0,  mgr->GetCommonInputContainer());
// mgr->ConnectOutput(taskESD,  1,  coutESD);
// mgr->ConnectOutput(taskESD,  2,  coutESD1);
 
  AliAnalysisDataContainer *coutESD = mgr->CreateContainer("output", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(taskESD, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskESD, 1, coutESD);
 
 return taskESD;
}
