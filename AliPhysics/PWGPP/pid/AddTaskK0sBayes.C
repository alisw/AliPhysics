AliAnalysisTask *AddTaskK0sBayes(Bool_t ismc=kFALSE,Bool_t qa=kTRUE,Int_t filterbit=4,Int_t typeCol=2,Bool_t toEP=kFALSE,Int_t species = 2){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("No manager found in AddTaskVZERO. Why?");
    return 0;
  }
  // currently don't accept AOD input
  if (!mgr->GetInputEventHandler()->InheritsFrom(AliAODInputHandler::Class())) {
    Error("AddTaskK0sBayes","This task works only with AOD input!");
    return 0;
  }

  //========= Add tender to the ANALYSIS manager and set default storage =====
  char mytaskName[100];
  snprintf(mytaskName,100,"AliAnalysisTaskK0sBayes.cxx"); 

  AliAnalysisTaskK0sBayes *task = new AliAnalysisTaskK0sBayes(mytaskName);
  if(ismc) task->SetMC();
  if(qa) task->SetQA();
  task->SetEtaCut(0.8);
  task->SetFilterBit(filterbit);
  task->SetTypeCollisions(typeCol);
  task->SetCorrEP(toEP);
  task->SetRefSpecies(species);

  AliPIDmaxProb *userCut = new AliPIDmaxProb("maxProbPion");
  task->SetPIDuserCut(userCut);
  userCut->RequireTPC();
  userCut->RequireTOF();
  mgr->AddTask(task);

  //Attach input to my tasks
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // Attach output to my tasks
  AliAnalysisDataContainer *cOutputL= mgr->CreateContainer("contK0sBayes1",TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 1, cOutputL);

  AliAnalysisDataContainer *cOutputL2= mgr->CreateContainer("contK0sBayes2",TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 2, cOutputL2);

  AliAnalysisDataContainer *cOutputL3= mgr->CreateContainer("contK0sBayes3",TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 3, cOutputL3);

  return task;
}

