AliAnalysisTask *AddTaskHFE(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFE", "No analysis manager found.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFE", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type=="AOD"){
    ::Error("AddTaskHFE", "The tasks exits because AODs are in input");
    return NULL;
  }

  //============= Set Task Name ===================
  //TString taskName=("AliAnalysisTaskHFE.cxx+");
  //===============================================

  AliHFEcuts *hfecuts = new AliHFEcuts;
  hfecuts->CreateStandardCuts();
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetMinNTrackletsTRD(1);
  //hfecuts->SetCheckITSLayerStatus(kFALSE);

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE("Heavy Flavour Electron Analysis");
  task->SetPIDStrategy(4);  
  task->SetHFECuts(hfecuts);
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);
  task->SetSecVtxOn();
  mgr->AddTask(task);

  //----------------------
  //create data containers
  //----------------------

  //find input container
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWG3_hfe";
  
  task->ConnectOutput(0, mgr->CreateContainer("nEvents", TH1I::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  task->ConnectOutput(1, mgr->CreateContainer("Results", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  task->ConnectOutput(2, mgr->CreateContainer("QA", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  mgr->ConnectInput  (task,  0, cinput );
  
  return task;
}
