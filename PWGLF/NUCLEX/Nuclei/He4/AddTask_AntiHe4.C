AliAnalysisTask *AddTask_AntiHe4(){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_AntiHe4", "No analysis manager found.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTask_AntiHe4", "This task requires an input event handler");
      return NULL;
   }  


  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskAntiHe4 *task = new AliAnalysisTaskAntiHe4("nmartinTaskAntiHe4");

  Int_t iResult = task->Initialize();
  if (!iResult)
    mgr->AddTask(task);
  else {
    //AliError("NO pt ranges specfied, not adding the task !!!");
    return -1;
  }

  //mgr->AddTask(task);
  
  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  //dumm output container
  AliAnalysisDataContainer *coutput0 =
      mgr->CreateContainer("nmartin_treeAntiHe4",
                           TTree::Class(),
                           AliAnalysisManager::kExchangeContainer,
                           "nmartin_default");

  //define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
      mgr->CreateContainer("nmartin_AntiHe4", TObjArray::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");

AliAnalysisDataContainer *coutput2 = 
      mgr->CreateContainer("treeAntiHe4", TTree::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");


  //connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  0, coutput0);
  mgr->ConnectOutput (task,  1, coutput1);
  mgr->ConnectOutput (task,  2, coutput2);

  return task;
}
