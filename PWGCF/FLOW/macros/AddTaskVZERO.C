AliAnalysisTask *AddTaskVZERO(AliAnalysisManager *mgr,Bool_t kV2=kTRUE,Bool_t kV3=kTRUE){
  char fileout[100];
  sprintf(fileout,"outVZEROv2.root");
  char fileout2[100];
  sprintf(fileout2,"outVZEROv3.root");

  //get the current analysis manager
  //  AnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("No manager found in AddTaskVZERO. Why?");
    return 0;
  }
  // currently don't accept AOD input
  if (!mgr->GetInputEventHandler()->InheritsFrom(AliAODInputHandler::Class())) {
    Error("AddTaskVZERO","This task works only with AOD input!");
    return 0;
  }

  //========= Add tender to the ANALYSIS manager and set default storage =====
  char mytaskName[100];
  sprintf(mytaskName,"AliAnalysisTaskVnV0.cxx"); 

  AliAnalysisTaskVnV0 *task = new AliAnalysisTaskVnV0(mytaskName);
  task->SetV2(kV2);
  task->SetV3(kV3);

  mgr->AddTask(task);

  //Attach input to my tasks
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // Attach output to my tasks
  if(kV2){
    AliAnalysisDataContainer *cOutputL= mgr->CreateContainer("contVZEROv2",TList::Class(), AliAnalysisManager::kOutputContainer, fileout);
    mgr->ConnectOutput(task, 1, cOutputL);
  }
  if(kV3){
    AliAnalysisDataContainer *cOutputL2= mgr->CreateContainer("contVZEROv3",TList::Class(), AliAnalysisManager::kOutputContainer, fileout2);
    mgr->ConnectOutput(task, 2, cOutputL2);
  }
  printf("task really added\n");

  return task;
}
