AliAnalysisTaskSE *AddTaskNuclei(Bool_t kAOD=kFALSE,Bool_t kMC=kFALSE,Float_t *kCentrality,Int_t filterBit,Int_t nTPCminCluster,Float_t DCAzCut,Float_t DCAxyCut,Bool_t bTPCcut=kTRUE){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  if (!mgr) {
    Error("No manager found in AddTaskVZERO. Why?");
    return 0;
  }
  // currently don't accept AOD input
  if (kAOD && !mgr->GetInputEventHandler()->InheritsFrom(AliAODInputHandler::Class())) { // check AOD
    Error("AddTaskVZERO","This task works only with AOD input!");
    return 0;
  }
  else if(1){ // check ESD

  }

  //========= Add tender to the ANALYSIS manager and set default storage =====
  char mytaskName[100];
  
  snprintf(mytaskName,100,"AliAnalysisNucleiMass%02i%02i",kCentrality[0],kCentrality[1]);

  AliAnalysisNucleiMass *task = new AliAnalysisNucleiMass(mytaskName);

  task->SetCentrality(kCentrality);
  task->SetFilterBit(filterBit);
  task->SetNTPCcluster(nTPCminCluster);
  task->SetDCAzCut(DCAzCut);
  task->SetDCAxyCut(DCAxyCut);
  task->SetkTPCcut(bTPCcut);

  mgr->AddTask(task);

  //Attach input to my tasks
  char name[200];
  
  //sprintf(name,"cchain1%02i%02i",kCentrality[0],kCentrality[1]);

  sprintf(name,"cchain1%02i%02i_FilterBit=%02i_NminTPCclusters=%03i_DCAzCUT=%.1f_DCAxyCUT=%.1f_kTPCcut=%i",kCentrality[0],kCentrality[1],filterBit,nTPCminCluster,DCAzCut,DCAxyCut,bTPCcut);
  
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(name,TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // Attach output to my tasks
  
  //sprintf(name,"Results_CC%02i%02i",kCentrality[0],kCentrality[1]);
  
  sprintf(name,"Results_CC%02i%02i_FilterBit=%02i_NminTPCclusters=%03i_DCAzCUT=%.1f_DCAxyCUT=%.1f_kTPCcut=%i",kCentrality[0],kCentrality[1],filterBit,nTPCminCluster,DCAzCut,DCAxyCut,bTPCcut);

  AliAnalysisDataContainer *cOutputL= mgr->CreateContainer(name,TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 1, cOutputL);

  return task;
}
