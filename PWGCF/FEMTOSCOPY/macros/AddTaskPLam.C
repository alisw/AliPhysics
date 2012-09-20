AliAnalysisTaskSE *AddTaskPLam(){
  // Adds the Proton-Lambda Femtoscopy task to the manager

  // Get the manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {printf("E-AddTaskPLam: Couldn't get the manager!\n");return 0;}

  // Add the proton-lambda task
  AliAnalysisTaskSE *taskPLam = new AliAnalysisTaskProtonLambda("TaskProtonLambda");
  if (!taskPLam){printf("E-AddTaskPLam: Couldn't create the task!\n");return 0;}
  UInt_t triggerMask=AliVEvent::kMB;
  triggerMask|=AliVEvent::kCentral;
  triggerMask|=AliVEvent::kSemiCentral;
  taskPLam->SelectCollisionCandidates(triggerMask); // std is AliVEvent::kMB
  mgr->AddTask(taskPLam);
  
  // Create containers for output
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("cHistLambdas", TList::Class(),    AliAnalysisManager::kOutputContainer, "ProtonLambda.AOD.root");
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cHistProtons", TList::Class(),    AliAnalysisManager::kOutputContainer, "ProtonLambda.AOD.root");
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("cHist2Part", TList::Class(),    AliAnalysisManager::kOutputContainer, "ProtonLambda.AOD.root");

  // Connect input/output
  mgr->ConnectInput(taskPLam, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPLam, 1, coutput);
  mgr->ConnectOutput(taskPLam, 2, coutput2);
  mgr->ConnectOutput(taskPLam, 3, coutput3);

  // Return the task
  return taskPLam;
}
