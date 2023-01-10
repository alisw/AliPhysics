AliAnalysisTaskMultspec *AddTaskMultspec(int part = 0, bool ismc = kTRUE, bool removePythiaGen = kTRUE, TString suffix = "")
{
  // analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { 
    ::Error("AddTaskMultspec", "No analysis manager to connect to."); 
    return NULL; 
  }
  if (!mgr->GetInputEventHandler()) { 
    ::Error("AddTaskMultspec", "This task requires an input event handler"); 
    return NULL; 
  }

  // Create the task and add it to the manager
  TString tskname = Form("Multspec_Task_%s", suffix.Data());
  AliAnalysisTaskMultspec *mytask = new AliAnalysisTaskMultspec(tskname, part, "", ismc, removePythiaGen);
  mytask->SetIsMC(ismc);
  mytask->SetRemovePythiaGen(removePythiaGen);
  mgr->AddTask(mytask);

  // output file name
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  //create and link only used containers
  AliAnalysisDataContainer *coutput[2];

  mgr->ConnectInput(mytask, 0, mgr->GetCommonInputContainer());
  coutput[0] = mgr->CreateContainer("fHistos_misc", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  coutput[1] = mgr->CreateContainer("fTree", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectOutput(mytask, 1, coutput[0]);
  mgr->ConnectOutput(mytask, 2, coutput[1]);
  
  return mytask;

}
