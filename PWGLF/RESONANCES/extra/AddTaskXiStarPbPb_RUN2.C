AliXiStarPbPb_RUN2 *AddTaskXiStarPbPb_RUN2(bool AODcase=kFALSE, bool MCcase=kFALSE, int CutList=0) {
 
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskXiStar", "No analysis manager to connect to.");
    return NULL;
  }
  
    if (AODcase == kTRUE) {
        Printf("INFO! You are using AODs!");}

  //____________________________________________//
  // Create tasks
  AliXiStarPbPb_RUN2 *XiStarTask_RUN2 = new AliXiStarPbPb_RUN2("XiStarTask_RUN2", AODcase, MCcase);
  if(!XiStarTask_RUN2) exit(-1);
    
  mgr->AddTask(XiStarTask_RUN2);
    
   
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGLF.outputXiStarAnalysis.root";
  AliAnalysisDataContainer *coutXiStar = mgr->CreateContainer("XiStarOutput_RUN2", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(XiStarTask_RUN2, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(XiStarTask_RUN2, 1, coutXiStar);
  
  
  // Return the task pointer
  return XiStarTask_RUN2;
}
