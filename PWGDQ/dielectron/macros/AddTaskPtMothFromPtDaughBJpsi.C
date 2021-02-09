AliAnalysisTask *AddTaskPtMothFromPtDaughBJpsi(Bool_t IsNtuplaCreated=kFALSE, Double_t minPt=0., Double_t maxPt=30., Int_t nbins=150, Double_t alfaExp=1.){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_ffionda_BJPsi", "No analysis manager found.");
    return 0;
  }
  ///
  AliPtMothFromPtDaugh *ptExtr = new AliPtMothFromPtDaugh();
  ptExtr->SetDefaultAnalysis(AliPtMothFromPtDaugh::kBtoJPSI);
  ptExtr->SetBinsPtMoth(minPt,maxPt,nbins,alfaExp);
  ptExtr->SetBinsPtMinMoth(minPt,maxPt,nbins,alfaExp);
  ///
  AliAnalysisTaskPtMothFromPtDaugh *task = new AliAnalysisTaskPtMothFromPtDaugh(IsNtuplaCreated);
  task->SetPtMothFromPtDaugh(ptExtr); // set AliPtMothFromPtDaugh object to the task
  task->SetNtuplaFileName("DecayKine.root");
 
  AliAnalysisDataContainer *cOutput = mgr->CreateContainer("Mothers", TList::Class(), AliAnalysisManager::kOutputContainer,"Mothers.root");
  mgr->ConnectOutput(task, 1, cOutput);
  // optional output for TNtupla
  AliAnalysisDataContainer *cOutput1 = 0x0;
  cOutput1 = mgr->CreateContainer("DecayKine", TNtuple::Class(), AliAnalysisManager::kOutputContainer,"DecayKine.root");
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,2,cOutput1);
  
  return task;
}
