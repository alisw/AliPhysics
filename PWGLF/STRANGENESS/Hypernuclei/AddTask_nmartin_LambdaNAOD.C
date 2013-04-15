AliAnalysisTask *AddTask_nmartin_LambdaNAOD(){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_nmartin_LambdaNAOD", "No analysis manager found.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTask_nmartin_LambdaNAOD", "This task requires an input event handler");
      return NULL;
   }  


  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskLambdaNAOD *task = new AliAnalysisTaskLambdaNAOD("nmartinTaskLambdaNAOD");

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  task->SetAnalysisType               (type);


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
      mgr->CreateContainer("nmartin_treeLambdaN",
                           TTree::Class(),
                           AliAnalysisManager::kExchangeContainer,
                           "nmartin_default");

  //define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
      mgr->CreateContainer("nmartin_LambdaNAOD", TObjArray::Class(),AliAnalysisManager::kOutputContainer,"nmartin_LambdaNAOD.root");

AliAnalysisDataContainer *coutput2 = 
      mgr->CreateContainer("treeLambdaNAOD", TTree::Class(),AliAnalysisManager::kOutputContainer,"nmartin_TreeLambdaNAOD.root");


  //connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  0, coutput0);
  mgr->ConnectOutput (task,  1, coutput1);
  mgr->ConnectOutput (task,  2, coutput2);

  return task;
}
