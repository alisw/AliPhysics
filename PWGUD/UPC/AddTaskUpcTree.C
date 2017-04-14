AliAnalysisTaskUpcTree *AddTaskUpcTree(
    const char* outputFileName = 0,
    const char* folderName = "UpcTree")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskUpcTree", "No analysis manager to connect to.");
    return NULL;
  }  

  AliAnalysisTaskUpcTree* ana = new AliAnalysisTaskUpcTree();
  mgr->AddTask(ana);
  
  if (!outputFileName) outputFileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histos", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("events", TTree::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1);
  mgr->ConnectOutput (ana, 2, coutput2);
  return ana;
}









