AliAnalysisTask *AddTask_Helium3PiMC(){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_Helium3Pi", "No analysis manager found.");
    return 0;
  }
  
  // mc event handler                                                                                                                                                                                                           
  AliMCEventHandler* mchandler = new AliMCEventHandler();
  // Not reading track references                                                                                                                                                                                           
  mchandler->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(mchandler);
    
  //========= Add task to the ANALYSIS manager =====

  AliAnalysisTaskSE *taskHelium3PiMC = new AliAnalysisTaskHelium3PiMC("Helium3PiMC_task");

  mgr->AddTask(taskHelium3PiMC);

  //================================================
  //              data containers
  //================================================
  //            find input container
 
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Helium3PiMC_tree", TTree::Class(), AliAnalysisManager::kOutputContainer, "He3Pi.Ntuple.MC.root");  

  //           connect containers
  mgr->ConnectInput  (taskHelium3PiMC,  0, cinput );
  mgr->ConnectOutput (taskHelium3PiMC,  1, coutput1);
  
  return taskHelium3PiMC;
}
