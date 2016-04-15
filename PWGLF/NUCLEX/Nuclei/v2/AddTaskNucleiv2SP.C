class AliAnalysisDataContainer;

AliAnalysisTask *AddTaskNucleiv2SP(TString name="name",Int_t ptctype =1,Bool_t dcacut = kFALSE, Float_t vzmax=10., TString CentEstim = "V0M", TString AnalysisType = "ESD", Bool_t isFlat = kTRUE, Int_t year = 2011, Int_t harmonic = 2){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNucleiv2", "No analysis manager found.");
    return 0;
  }
 
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskNucleiv2SP *task = new AliAnalysisTaskNucleiv2SP(name.Data());
  task->SetIsPrimCut(dcacut); 
  task->SetParticle(ptctype);
  task->SetVzMax(vzmax);
  task->SetCentralityEstimator(CentEstim);
  task->SetAnalysisType(AnalysisType);
  task->SetApplyFlatten(isFlat);
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

