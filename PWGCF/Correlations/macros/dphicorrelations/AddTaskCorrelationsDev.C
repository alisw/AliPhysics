AliAnalysisTask* AddTaskCorrelationsDev(const char* containerName = "histograms")
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhiCorrelations", "No analysis manager to connect to.");
    return NULL;
  }  
  
  AliAnalysisTaskCorrelationsDev* ana = new  AliAnalysisTaskCorrelationsDev(containerName);

  ana->SetEventSelectionBit(0);
  ana->SetDebugLevel(0); 
  ana->SetFilterBit(0);
  
  ana->SetTrackEtaCut(0.9);
  ana->SetPtMin(0);

//   ana->SetZVertex(10);
  ana->SetCentralityMethod("");
  
  mgr->AddTask(ana);
  
  const char* outputFileName = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName, TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, "PWGCF_CorrelationsDev"));
  
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1 );
   
  return ana;
}
