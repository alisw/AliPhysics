AliAnalysisTaskMuonHadronCorrelations *AddAnalysisTaskMuonHadronCorrelations(const char *centMethod = "V0M") {
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    printf("Error in adding AnalysisTaskMuonHadronCorrelations: no Analysis Manager found!\n");
    return NULL;
  }

  AliAnalysisTaskMuonHadronCorrelations *task = new AliAnalysisTaskMuonHadronCorrelations(Form("AliAnalysisTaskMuonHadronCorrelations_%s",centMethod));

  // Set analysis cuts   
  task->SetFilterBitCentralBarrel(7);  // -> 128
  task->SetMaxEtaCentralBarrel(1.0);
  task->SetMinEtaCentralBarrel(0.0);
  task->SetTriggerMatchLevelMuon(1);

  const Int_t nBinCent = 4;
  Double_t centLimits[nBinCent+1] = {0., 20., 40, 60., 100.};
  task->SetCentBinning(nBinCent, centLimits);

  task->SetCentMethod(centMethod);

  const Int_t nBinPt = 3;
  Double_t ptLimits[nBinPt+1] = {0., 1., 2., 4.};
  task->SetPtBinning(nBinPt, ptLimits);

  const Int_t nBinEta = 3;
  Double_t etaLimits[nBinEta+1] = {-4., -3.6, -3.2, -2.5};
  task->SetEtaBinning(nBinEta, etaLimits);

  mgr->AddTask(task);

  // create output container
  AliAnalysisDataContainer *output = mgr->CreateContainer(Form("MuonHadronCorrHistos_%s",centMethod), TList::Class(), AliAnalysisManager::kOutputContainer,
							  Form("%s:MuonHadronCorrelations_%s", AliAnalysisManager::GetCommonFileName(), centMethod));
  
  // finaly connect input and output
  mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);
    
  return task;
}

