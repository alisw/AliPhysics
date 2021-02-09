AliTOFTriggerEfficiencyTask *AddTOFTriggerEfficiencyTask(
    const char* outputFileName = 0,
    const char* folderName = "TOFTriggerEfficiency")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTOFTriggerEfficiencyTask", "No analysis manager to connect to.");
    return NULL;
  }  

  AliTOFTriggerEfficiencyTask* ana = new AliTOFTriggerEfficiencyTask();
  mgr->AddTask(ana);
  
  if (!outputFileName) outputFileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histos", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("events", TTree::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1);
  mgr->ConnectOutput (ana, 2, coutput2);
  return ana;
}









