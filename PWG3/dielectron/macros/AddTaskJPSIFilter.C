AliAnalysisTask *AddTaskJPSIFilter(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jpsi_DielectronFilter", "No analysis manager found.");
    return 0;
  }
  // currently don't accept AOD input
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()) {
    Warning("AddTask_jpsi_JPsi","No AOD input supported currently. Not adding the task!");
    return 0;
  }
  //check for output aod handler
  if (!mgr->GetOutputEventHandler()||mgr->GetOutputEventHandler()->IsA()!=AliAODHandler::Class()) {
    Warning("AddTask_jpsi_DielectronFilter","No AOD output handler available. Not adding the task!");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //Create task and add it to the analysis manager
  AliAnalysisTaskDielectronFilter *task=new AliAnalysisTaskDielectronFilter("jpsi_DielectronFilter");
  
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/dielectron/macros/ConfigJpsi2eeFilter.C");
  AliDielectron *jpsi=ConfigJpsi2eeFilter();
  if (!hasMC) task->UsePhysicsSelection();
  task->SetDielectron(jpsi);
  mgr->AddTask(task);

  //----------------------
  //create data containers
  //----------------------
  
  
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWG3_dielectronFilter";
  
  //create output container
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("QA", TList::Class(), AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("CF", TList::Class(), AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);

  return task;
}
