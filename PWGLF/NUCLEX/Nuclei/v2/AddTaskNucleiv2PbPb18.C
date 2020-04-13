class AliAnalysisDataContainer;

AliAnalysisTask *AddTaskNucleiv2PbPb18(TString name="name",Int_t ptctype =1, Float_t vzmax=10.,  Int_t harmonic = 2, Int_t nsigma = 3, Short_t filterbit = 4, Bool_t isLHC18r = 1 , Short_t CentEstim = 0, Short_t recPass = 0){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNucleiv2", "No analysis manager found.");
    return 0;
  }
 
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskNucleiv2PbPb18 *task = new AliAnalysisTaskNucleiv2PbPb18(name.Data());
  task->SetParticle(ptctype);
  task->SetVzMax(vzmax);
  task->SetCentralityEstimator(CentEstim);
  task->SetHarmonic(harmonic);
  task->SetRecPass(recPass);
  task->SetFilterBit(filterbit);
  task->SetPeriod(isLHC18r);
  task->SetTPCnsigma(nsigma);
  mgr->AddTask(task);

  //================================================
  //              data containers
  //================================================
  //            find input container

  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clisthist", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("treeNuclei", TTree::Class(),AliAnalysisManager::kOutputContainer, outputFileName);
  //           connect containers

  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput1);
  mgr->ConnectOutput (task,  2, coutput2);

  return task;
}
