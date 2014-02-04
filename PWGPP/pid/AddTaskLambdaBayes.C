AliAnalysisTask *AddTaskLambdaBayes(Bool_t ismc=kFALSE,Bool_t qa=kTRUE,Int_t filterbit=4,Int_t typeCol=2,Bool_t toEP=kFALSE,Int_t species=4){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("No manager found in AddTaskVZERO. Why?");
    return 0;
  }
  // currently don't accept AOD input
  if (!mgr->GetInputEventHandler()->InheritsFrom(AliAODInputHandler::Class())) {
    Error("AddTaskLambdaBayes","This task works only with AOD input!");
    return 0;
  }

  //========= Add tender to the ANALYSIS manager and set default storage =====
  char mytaskName[100];
  snprintf(mytaskName,100,"AliAnalysisTaskLambdaBayes.cxx"); 

  AliAnalysisTaskLambdaBayes *task = new AliAnalysisTaskLambdaBayes(mytaskName);
  if(ismc) task->SetMC();
  if(qa) task->SetQA();
  task->SetEtaCut(0.8);
  task->SetFilterBit(filterbit);
  task->SetTypeCollisions(typeCol);
  task->SetCorrEP(toEP);
  task->SetRefSpecies(species);

  AliPIDmaxProb *userCut = new AliPIDmaxProb("maxProbProton");
  userCut->RequireTPC();
  userCut->RequireTOF();
  task->SetPIDuserCut(userCut);

  mgr->AddTask(task);

  //Attach input to my tasks
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // Attach output to my tasks
  AliAnalysisDataContainer *cOutputL= mgr->CreateContainer("contLambdaBayes1",TList::Class(), AliAnalysisManager::kOutputContainer, 
AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 1, cOutputL);

  AliAnalysisDataContainer *cOutputL2= mgr->CreateContainer("contLambdaBayes2",TList::Class(), AliAnalysisManager::kOutputContainer, 
AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 2, cOutputL2);

  AliAnalysisDataContainer *cOutputL3= mgr->CreateContainer("contLambdaBayes3",TList::Class(), AliAnalysisManager::kOutputContainer, 
AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 3, cOutputL3);

  return task;
}

