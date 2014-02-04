AliAnalysisTask *AddTask_Helium3Pi(){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_Helium3Pi", "No analysis manager found.");
    return 0;
  }
  
  //========= Add task to the ANALYSIS manager =====

  AliAnalysisTaskSE *taskHelium3Pi = new AliAnalysisTaskHelium3Pi("Helium3Pi_task");

  mgr->AddTask(taskHelium3Pi);

  //================================================
  //              data containers
  //================================================
  //            find input container
 
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Helium3Pi_tree", TTree::Class(), AliAnalysisManager::kOutputContainer, "He3Pi.Ntuple.root");  

  //           connect containers
  mgr->ConnectInput  (taskHelium3Pi,  0, cinput );
  mgr->ConnectOutput (taskHelium3Pi,  1, coutput1);
  
  return taskHelium3Pi;
}
