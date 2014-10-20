AliAnalysisTask *AddTaskNucleiv2(TString name="name",TString eventtype="REAL", Bool_t saveTree = kFALSE){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_Helium3Pi", "No analysis manager found.");
    return 0;
  }
  
  //========= Add task to the ANALYSIS manager =====

  AliAnalysisTaskNucleiv2 *task = new   AliAnalysisTaskNucleiv2(name,eventtype,saveTree);
  task->SetDataType(eventtype);
  mgr->AddTask(task);
  
  //================================================
  //              data containers
  //================================================
  //            find input container

  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Helium3Pi_tree", TTree::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");  
  //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Helium3Pi_tree", TTree::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");  

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clisthist", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("treeNuclei", TTree::Class(),AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("treeMC" , TTree::Class(),AliAnalysisManager::kOutputContainer   , outputFileName);

    //           connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput1);
  mgr->ConnectOutput (task,  2, coutput2);
  mgr->ConnectOutput (task,  3, coutput3);

  return task;
}
