AliAnalysisTask *AddTaskXtAnalysis() {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFluctuations", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFluctuations", "This task requires an input event handler");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  
  // parameter provides the location of the input card
  AliXtAnalysis *xtTask = new AliXtAnalysis("AliXtAnalysis",Form("%s%s",gSystem->Getenv("ALICE_ROOT"),"/PWGLF/SPECTRA/XtAnalysis/card_xT.input"));
  xtTask->SetDebugLevel(0);
  xtTask->SetDebugMode(0);  

  mgr->AddTask((AliAnalysisTask*) xtTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("xtAnalysis", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:xtAnalysis",AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("xtHistos", TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:xtHistos",AliAnalysisManager::GetCommonFileName()));
  // Connect input/output
  mgr->ConnectInput(xtTask, 0, cinput);
  mgr->ConnectOutput(xtTask, 1, coutput1);
  mgr->ConnectOutput(xtTask, 2, coutput2);

  return xtTask;
}

