AliAnalysisTaskSE* AddTaskSimple(const char* containerName = "")
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCorrelationsDev", "No analysis manager to connect to.");
    return NULL;
  }  
  
  AliAnalysisTaskNanoSimple* ana = new AliAnalysisTaskNanoSimple(Form("AliAnalysisTaskNanoSimple%s", containerName));
  mgr->AddTask(ana);

  const char* outputFileName = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("histograms%s", containerName), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, "Simple"));
  
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1 );
   
  return ana;
}
