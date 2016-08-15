AliAnalysisTaskEmcalTriggerPosition *AddTaskEmcalTriggerPosition(){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisTaskEmcalTriggerPosition *task = new AliAnalysisTaskEmcalTriggerPosition("triggerpos");
  mgr->AddTask(task);

  TString outputfile = mgr->GetCommonFileName();
  outputfile += ":emcaltriggerpos";

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("EMCALtriggerPosition", TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));
  return task;
}


