class AliAnalysisDataContainer;

AliAnalysisTask *AddTaskHypv2PbPb18(TString name="name", Float_t vzmax=10., Short_t CentEstim = 0, Short_t recPass = 0, Int_t harmonic = 2, Bool_t isLHC18r = 1){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHypv2", "No analysis manager found.");
    return 0;
  }
 
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskHypv2PbPb18 *task = new AliAnalysisTaskHypv2PbPb18(name.Data());
  task->SetVzMax(vzmax);
  task->SetCentralityEstimator(CentEstim);
  task->SetHarmonic(harmonic);
  task->SetPeriod(isLHC18r);
  mgr->AddTask(task);

  //================================================
  //              data containers
  //================================================
  //            find input container

  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clisthist", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("treeHyp", TTree::Class(),AliAnalysisManager::kOutputContainer, outputFileName);
  //           connect containers

  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput1);
  mgr->ConnectOutput (task,  2, coutput2);

  return task;
}
