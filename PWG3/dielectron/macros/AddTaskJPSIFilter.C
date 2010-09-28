AliAnalysisTask *AddTaskJPSIFilter(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskJPSIFilter", "No analysis manager found.");
    return 0;
  }
  
  //check for output aod handler
  if (!mgr->GetOutputEventHandler()||mgr->GetOutputEventHandler()->IsA()!=AliAODHandler::Class()) {
    Error("AddTaskJPSIFilter","No AOD output handler available. Cannot Proceed!");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //Do we run on AOD?
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  //Create task and add it to the analysis manager
  AliAnalysisTaskDielectronFilter *task=new AliAnalysisTaskDielectronFilter("jpsi_DielectronFilter");
  
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/dielectron/macros/ConfigJpsi2eeFilter.C");
  AliDielectron *jpsi=ConfigJpsi2eeFilter(isAOD);
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
    mgr->CreateContainer("jpsi_FilterQA",
                         THashList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("jpsi_FilterEventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  
  return task;
}
