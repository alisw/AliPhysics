class AliAnalysisDataContainer;

AliAnalysisTask *AddTaskAllPtcv2(TString name="name",Int_t ptctype =1,Bool_t dcacut = kFALSE, Float_t vzmax=10., TString CentEstim = "V0M", TString AnalysisType = "ESD", Int_t year = 2015, Int_t harmonic = 2){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskAllPtcv2", "No analysis manager found.");
    return 0;
  }
 
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskAllPtcv2 *task = new AliAnalysisTaskAllPtcv2(name.Data());
  task->SetIsPrimCut(dcacut); 
  task->SetParticle(ptctype);
  task->SetVzMax(vzmax);
  task->SetCentralityEstimator(CentEstim);
  task->SetAnalysisType(AnalysisType);
  task->SetYear(year);
  task->SetHarmonic(harmonic);
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

