AliAnalysisTask *AddTask_Helium3Pi(TString name="name"){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_Helium3Pi", "No analysis manager found.");
    return 0;
  }
  
  //========= Add task to the ANALYSIS manager =====

  AliAnalysisTaskHelium3PiAOD *taskHelium3Pi = new AliAnalysisTaskHelium3PiAOD(name);
   
  mgr->AddTask(taskHelium3Pi);
  
  //================================================
  //              data containers
  //================================================
  //            find input container

  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Helium3Pi_tree", TTree::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");  
  //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Helium3Pi_tree", TTree::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");  
 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clisthistHyper", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("treeHyper", TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("treeHelium" , TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  //           connect containers
  mgr->ConnectInput  (taskHelium3Pi,  0, cinput );
  mgr->ConnectOutput (taskHelium3Pi,  1, coutput1);
  mgr->ConnectOutput (taskHelium3Pi,  2, coutput2);
  mgr->ConnectOutput (taskHelium3Pi,  3, coutput3);

  return taskHelium3Pi;
}
