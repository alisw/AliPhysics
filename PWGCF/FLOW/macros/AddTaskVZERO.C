AliAnalysisTask *AddTaskVZERO(Bool_t ismc=kFALSE,Bool_t kV2=kTRUE,Bool_t kV3=kTRUE,Bool_t qa=kTRUE,Bool_t modulationdEdx=kFALSE,Bool_t globalTrack=kFALSE){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
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
  snprintf(mytaskName,100,"AliAnalysisTaskVnV0.cxx"); 

  AliAnalysisTaskVnV0 *task = new AliAnalysisTaskVnV0(mytaskName);
  task->SetV2(kV2);
  task->SetV3(kV3);
  if(ismc) task->SetMC();
  if(qa) task->SetQA();

  task->SetFillDCAinfo(globalTrack); // 0 = TPC only track, 1 = global tracks
  task->SetModulationDEDx(modulationdEdx);

  mgr->AddTask(task);

  //Attach input to my tasks
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // Attach output to my tasks
  if(kV2){
    AliAnalysisDataContainer *cOutputL= mgr->CreateContainer("contVZEROv2",TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
    mgr->ConnectOutput(task, 1, cOutputL);
  }
  if(kV3){
    AliAnalysisDataContainer *cOutputL2= mgr->CreateContainer("contVZEROv3",TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
    mgr->ConnectOutput(task, 2, cOutputL2);
  }
  if(ismc){
    AliAnalysisDataContainer *cOutputL3= mgr->CreateContainer("contVZEROmc",TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
    mgr->ConnectOutput(task, 3, cOutputL3);
  }
  if(qa){
    AliAnalysisDataContainer *cOutputL4= mgr->CreateContainer("contVZEROqa",TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
    mgr->ConnectOutput(task, 4, cOutputL4);
  }
  printf("task really added\n");

  return task;
}

