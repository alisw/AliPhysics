AliXiStarpp13TeV *AddTaskXiStarpp13TeV(bool MCcase=kFALSE, bool AODcase=kFALSE, int CutList=0, bool DevelopmentMode=kFALSE, bool HMTrigger=kFALSE, bool PIDOption=kFALSE, bool SetSystematic=kTRUE){

  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }

  //____________________________________________//
  // Create tasks
  AliXiStarpp13TeV *XiStarTask = new AliXiStarpp13TeV("XiStarTask", AODcase, MCcase, CutList, DevelopmentMode, HMTrigger, PIDOption, SetSystematic);
  if(!XiStarTask) exit(-1);
  mgr->AddTask(XiStarTask);


  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGLF.outputXiStarAnalysis.root";
  AliAnalysisDataContainer *coutXiStar = mgr->CreateContainer(Form("MyList_%d", HMTrigger), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(XiStarTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(XiStarTask, 1, coutXiStar);


  // Return the task pointer
  return XiStarTask;
}
