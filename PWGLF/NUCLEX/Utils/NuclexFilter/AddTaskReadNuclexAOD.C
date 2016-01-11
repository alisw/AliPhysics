AliAnalysisTask *AddTaskReadNuclexAOD(TString name="name"){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskReadNuclexAOD", "No analysis manager found.");
    return 0;
  }
  
  //========= Add task to the ANALYSIS manager =====

  AliAnalysisTaskReadNuclexAOD *taskReadNuclexAOD = new AliAnalysisTaskReadNuclexAOD(name);
   
  mgr->AddTask(taskReadNuclexAOD);
  
  //================================================
  //              data containers
  //================================================
  //            find input container

  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Helium3Pi_tree", TTree::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");  
  //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Helium3Pi_tree", TTree::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");  
 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clisthistHyper"  , TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("treeAODrecoDecay", TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("treeMySecVert"   , TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  //           connect containers
  mgr->ConnectInput  (taskReadNuclexAOD,  0, cinput );
  mgr->ConnectOutput (taskReadNuclexAOD,  1, coutput1);
  mgr->ConnectOutput (taskReadNuclexAOD,  2, coutput2);
  mgr->ConnectOutput (taskReadNuclexAOD,  3, coutput3);

  return taskReadNuclexAOD;
}
