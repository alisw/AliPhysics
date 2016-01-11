AliAnalysisTaskSE *AddTaskNuclei(Bool_t kAOD=kTRUE, Double_t CentralityMin=0.0, Double_t CentralityMax=100.0, Int_t filterBit=16, Double_t EtaMin=-0.8, Double_t EtaMax=0.8, Double_t DCAxyCut=0.1, Double_t DCAzCut=1000.0, Double_t fNsigmaTpcCut=2.0, Int_t NminTpcCluster=0, Int_t iTRDslices=0, Int_t kSignalCheck=1, Int_t iMtof=1, Int_t kPvtxCorr=1){

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
  task->SetCentrality(CentralityMin,CentralityMax);
  task->SetFilterBit(filterBit);
  task->SetEtaLimit(EtaMin,EtaMax);
  task->SetDCACut(DCAxyCut,DCAzCut);
  task->SetNsigmaTPCCut(fNsigmaTpcCut);
  task->SetNminTPCcluster(NminTpcCluster);
  task->SetTrdCut(iTRDslices);
  task->SetisSignalCheck(kSignalCheck);
  task->SetMtofMethod(iMtof);
  task->SetPvtxNucleiCorrection(kPvtxCorr);

  mgr->AddTask(task);

  //Attach input to my tasks
  char name[1000];

  snprintf(name,1000,"cchain1%.0f_%.0f_FilterBit=%02i_EtaMin=%.1f_EtaMax=%.1f_DCAxyCUT=%.2f_DCAzCUT=%.1f_NsigTPCcut=%1.0f_NminTpcClusters=%03i_iTrdCut=%i_kSignCheck=%i_iMtof=%i_kPvtxCorr=%i",CentralityMin,CentralityMax,filterBit,EtaMin,EtaMax,DCAxyCut,DCAzCut,fNsigmaTpcCut,NminTpcCluster,iTRDslices,kSignalCheck,iMtof,kPvtxCorr);

  AliAnalysisDataContainer *cinput = mgr->CreateContainer(name,TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // Attach output to my tasks
  
  snprintf(name,1000,"ResultsBmm_CC%.0f_%.0f_FilterBit=%02i_EtaMin=%.1f_EtaMax=%.1f_DCAxyCUT=%.2f_DCAzCUT=%.1f_NsigTPCcut=%1.0f_NminTpcClusters=%03i_iTrdCut=%i_kSignCheck=%i_iMtof=%i_kPvtxCorr=%i",CentralityMin,CentralityMax,filterBit,EtaMin,EtaMax,DCAxyCut,DCAzCut,fNsigmaTpcCut,NminTpcCluster,iTRDslices,kSignalCheck,iMtof,kPvtxCorr);
  AliAnalysisDataContainer *cOutputL= mgr->CreateContainer(name,TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 1, cOutputL);

  snprintf(name,1000,"ResultsBpp_CC%.0f_%.0f_FilterBit=%02i_EtaMin=%.1f_EtaMax=%.1f_DCAxyCUT=%.2f_DCAzCUT=%.1f_NsigTPCcut=%1.0f_NminTpcClusters=%03i_iTrdCut=%i_kSignCheck=%i_iMtof=%i_kPvtxCorr=%i",CentralityMin,CentralityMax,filterBit,EtaMin,EtaMax,DCAxyCut,DCAzCut,fNsigmaTpcCut,NminTpcCluster,iTRDslices,kSignalCheck,iMtof,kPvtxCorr);
  AliAnalysisDataContainer *cOutputL2= mgr->CreateContainer(name,TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 2, cOutputL2);

  return task;
}
