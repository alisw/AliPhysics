class AliAnalysisDataContainer;

AliAnalysisTask *AddTaskNucleiv2SP(TString name="name"){
  
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNucleiv2", "No analysis manager found.");
    return 0;
  }
  
  //========= Add task to the ANALYSIS manager =====

  AliAnalysisTaskNucleiv2SP *task = new   AliAnalysisTaskNucleiv2SP(name);

  //================================================
  //              data containers
  //================================================
  //            find input container

  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clisthist", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
 
    //           connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput1);
 
  mgr->AddTask(task);
 
  return task;
}

