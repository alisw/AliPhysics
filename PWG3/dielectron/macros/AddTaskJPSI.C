AliAnalysisTask *AddTaskJPSI(const char* config=""){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTasJPSI", "No analysis manager found.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJPSI", "This task requires an input event handler");
    return NULL;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(mgr->GetMCtruthEventHandler()!=0x0);
  
  TString configFile("$ALICE_ROOT/PWG3/dielectron/macros/ConfigJpsi2eeData.C");
  if (hasMC) configFile="$ALICE_ROOT/PWG3/dielectron/macros/ConfigJpsi2eeEff.C";
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDie");
  mgr->AddTask(task);
  
  //load dielectron configuration file
  gROOT->LoadMacro(configFile.Data());
  
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigJpsi2ee(i,isAOD);
    if (jpsi) task->AddDielectron(jpsi);
  }
  
  //----------------------
  //create data containers
  //----------------------
  
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWG3_dielectron";
    
  //create output container
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("jpsi_QA", TList::Class(), AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("jpsi_CF", TList::Class(), AliAnalysisManager::kOutputContainer,
                         containerName.Data());

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("jpsi_EventStat", TH1D::Class(), AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
