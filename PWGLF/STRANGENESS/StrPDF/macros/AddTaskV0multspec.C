AliAnalysisTaskV0multspec *AddTaskV0multspec(int part = 0, bool ismc = kFALSE, TString suffix = "")
{
  // analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { 
    ::Error("AddTaskV0multspec", "No analysis manager to connect to."); 
    return NULL; 
  }
  if (!mgr->GetInputEventHandler()) { 
    ::Error("AddTaskV0multspec", "This task requires an input event handler"); 
    return NULL; 
  }

  // Create the task and add it to the manager
  TString tskname = Form("StrMultSpec_%s", suffix.Data());
  AliAnalysisTaskV0multspec *mytask = new AliAnalysisTaskV0multspec(tskname, part, "", ismc);
  mytask->SetIsMC(ismc);
  mgr->AddTask(mytask);

  // output file name
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  //create and link only used containers
  AliAnalysisDataContainer *coutput[2];

  mgr->ConnectInput(mytask, 0, mgr->GetCommonInputContainer());
  coutput[0] = mgr->CreateContainer(Form("fHistos_misc_%s",suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  coutput[1] = mgr->CreateContainer(Form("fTree_%s",suffix.Data()),      TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectOutput(mytask, 1, coutput[0]);
  mgr->ConnectOutput(mytask, 2, coutput[1]);
  
  return mytask;

}
