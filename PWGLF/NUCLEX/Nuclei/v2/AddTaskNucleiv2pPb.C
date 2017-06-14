class AliAnalysisDataContainer;

AliAnalysisTask *AddTaskNucleiv2pPb(TString name="name",Int_t ptctype =1, Float_t vzmax=10., Short_t CentEstim = 0, Short_t recPass = 0, Int_t harmonic = 2, Short_t filterbit = 4){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNucleiv2", "No analysis manager found.");
    return 0;
  }
 
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskNucleiv2pPb *task = new AliAnalysisTaskNucleiv2pPb(name.Data());
  task->SetParticle(ptctype);
  task->SetVzMax(vzmax);
  task->SetCentralityEstimator(CentEstim);
  task->SetHarmonic(harmonic);
  task->SetRecPass(recPass);
  task->SetFilterBit(filterbit);
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
