AliAnalysisTask *AddTaskReadAODBranch(Double_t ptLegCut = 1., Bool_t spdFirstRequired=kFALSE, Int_t numClsTPC=90, Int_t pairType = 1, Double_t ptJpsi = 1.3){
  
  Bool_t hasMC = kFALSE;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskReadAODBranch", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskReadAODBranch", "This task requires an input event handler");
    return NULL;
  }

 AliAnalysisTaskDielectronReadAODBranch *readAODBranch = new AliAnalysisTaskDielectronReadAODBranch("ReadAODBranch");
 readAODBranch->SetHasMC(hasMC);
 readAODBranch->SetPtLeg(ptLegCut);
 readAODBranch->SetSpdFirstRequired(spdFirstRequired);
 readAODBranch->SetNclsTPC(numClsTPC); 
 readAODBranch->SetInvMassSignalRegion(2.3,4.);
 readAODBranch->SetInvMassSidebandRegion(2.9,3.2);
 readAODBranch->SetPairType(pairType);
 readAODBranch->SetPtJpsi(ptJpsi);

 mgr->AddTask(readAODBranch);

 AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("resultAOD",
                         TList::Class(), AliAnalysisManager::kOutputContainer,"result.root");

 mgr->ConnectInput(readAODBranch,  0, mgr->GetCommonInputContainer());
 mgr->ConnectOutput(readAODBranch, 1, cOutputHist);

 return readAODBranch;
}
