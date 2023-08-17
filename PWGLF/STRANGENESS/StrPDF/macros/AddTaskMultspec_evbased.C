AliAnalysisTaskMultspec_evbased *AddTaskMultspec_evbased(int part = 0, bool ismc = kFALSE, bool removePythiaGen = kTRUE, TString suffix = "")
{

  printf("HERE -------------------\n");
  // analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { 
    ::Error("AliAnalysisTaskMultspec_evbased", "No analysis manager to connect to.");
    return NULL; 
  }
  if (!mgr->GetInputEventHandler()) { 
    ::Error("AliAnalysisTaskMultspec_evbased", "This task requires an input event handler");
    return NULL; 
  }

  // Create the task and add it to the manager
  TString tskname = Form("Multspec_evbased_%s", suffix.Data());
  AliAnalysisTaskMultspec_evbased *mytask = new AliAnalysisTaskMultspec_evbased(tskname, part, "", ismc, removePythiaGen);
  mytask->SetIsMC(ismc);
  mytask->SetRemovePythiaGen(removePythiaGen);
  mgr->AddTask(mytask);

  // output file name
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  //create and link only used containers
  AliAnalysisDataContainer *coutput[2];

  mgr->ConnectInput(mytask, 0, mgr->GetCommonInputContainer());
  coutput[0] = mgr->CreateContainer("fHistos", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  // coutput[1] = mgr->CreateContainer("fTree", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectOutput(mytask, 1, coutput[0]);
  // mgr->ConnectOutput(mytask, 2, coutput[1]);
  
  return mytask;

}
