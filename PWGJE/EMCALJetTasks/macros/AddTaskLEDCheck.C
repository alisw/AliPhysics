AliAnalysisTaskLEDCheck *AddTaskLEDCheck(){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisTaskLEDCheck *ledtask = new AliAnalysisTaskLEDCheck("LEDTASK");
  mgr->AddTask(ledtask);

  TString outputfile = mgr->GetCommonFileName();
  outputfile += ":LEDCheck";

  mgr->ConnectOutput(ledtask, 0, mgr->CreateContainer("LEDCheck", TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));

  return ledtask;
}
