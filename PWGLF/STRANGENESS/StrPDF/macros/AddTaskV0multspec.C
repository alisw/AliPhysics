AliAnalysisTaskV0multspec *AddTaskV0multspec(int part = 0, bool addCasc = kFALSE, bool ismc = kFALSE, TString suffix = "")
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
  TString tskname = Form("V0multspec_Task_%s", suffix.Data());
  AliAnalysisTaskV0multspec *mytask = new AliAnalysisTaskV0multspec(tskname, part, addCasc, "", ismc);
  mytask->SetIsMC(ismc);
  mgr->AddTask(mytask);

  // output file name
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  //create and link only used containers
  AliAnalysisDataContainer *coutput[3];

  mgr->ConnectInput(mytask, 0, mgr->GetCommonInputContainer());
  coutput[0] = mgr->CreateContainer("fHistos_misc", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  coutput[1] = mgr->CreateContainer("fTreeV0",      TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  if(addCasc) coutput[2] = mgr->CreateContainer("fTreeCascade", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectOutput(mytask, 1, coutput[0]);
  mgr->ConnectOutput(mytask, 2, coutput[1]);
  if(addCasc) mgr->ConnectOutput(mytask, 3, coutput[2]);
  
  return mytask;

}
