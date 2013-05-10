AliAnalysisTaskSE *AddTaskNuclei(Bool_t kAOD=kFALSE,Bool_t kMC=kFALSE,Float_t fCentralityMin=0.0,Float_t fCentralityMax=100.0,Int_t filterBit,Int_t nTPCminCluster,Float_t DCAzCut,Float_t DCAxyCut,Bool_t bTPCcut=kTRUE,Float_t fNsigmaTpcCut,Bool_t bSignalCheck=kTRUE){

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
  
  snprintf(mytaskName,100,"AliAnalysisNucleiMass");

  AliAnalysisNucleiMass *task = new AliAnalysisNucleiMass(mytaskName);

  Float_t *kCentrality = new Float_t[2];
  kCentrality[0]=fCentralityMin;
  kCentrality[1]=fCentralityMax;  

  task->SetCentrality(kCentrality);
  task->SetFilterBit(filterBit);
  task->SetNTPCcluster(nTPCminCluster);
  task->SetDCAzCut(DCAzCut);
  task->SetDCAxyCut(DCAxyCut);
  task->SetkTPCcut(bTPCcut);
  task->SetNsigmaTPCCut(fNsigmaTpcCut);
  task->SetisSignalCheck(bSignalCheck);
  
  
  mgr->AddTask(task);

  //Attach input to my tasks
  char name[200];
  
  //sprintf(name,"cchain1%02i%02i",kCentrality[0],kCentrality[1]);

  sprintf(name,"cchain1%02i%02i_FilterBit=%02i_NminTPCclusters=%03i_DCAzCUT=%.1f_DCAxyCUT=%.1f_kTPCcut=%i_NsigTPCcut=%1.0f_bSignCheck=%i",kCentrality[0],kCentrality[1],filterBit,nTPCminCluster,DCAzCut,DCAxyCut,bTPCcut,fNsigmaTpcCut,bSignalCheck);
  
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(name,TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // Attach output to my tasks
  
  //sprintf(name,"Results_CC%02i%02i",kCentrality[0],kCentrality[1]);
  
  sprintf(name,"ResultsBmm_CC%02i%02i_FilterBit=%02i_NminTPCclusters=%03i_DCAzCUT=%.1f_DCAxyCUT=%.1f_kTPCcut=%i_NsigTPCcut=%1.0f_bSignCheck=%i",kCentrality[0],kCentrality[1],filterBit,nTPCminCluster,DCAzCut,DCAxyCut,bTPCcut,fNsigmaTpcCut,bSignalCheck);
  AliAnalysisDataContainer *cOutputL= mgr->CreateContainer(name,TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 1, cOutputL);

  sprintf(name,"ResultsBpp_%02i%02i_FilterBit=%02i_NminTPCclusters=%03i_DCAzCUT=%.1f_DCAxyCUT=%.1f_kTPCcut=%i_NsigTPCcut=%1.0f_bSignCheck=%i",kCentrality[0],kCentrality[1],filterBit,nTPCminCluster,DCAzCut,DCAxyCut,bTPCcut,fNsigmaTpcCut,bSignalCheck);
  AliAnalysisDataContainer *cOutputL2= mgr->CreateContainer(name,TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 2, cOutputL2);

  return task;
}
