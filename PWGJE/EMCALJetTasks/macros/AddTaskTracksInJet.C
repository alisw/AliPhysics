PWGJE::EMCALJetTasks::AliAnalysisTaskTracksInJet * AddTaskTracksInJet(Bool_t isMC){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  PWGJE::EMCALJetTasks::AliAnalysisTaskTracksInJet *trackInJet = new PWGJE::EMCALJetTasks::AliAnalysisTaskTracksInJet("trackInJet");
  if(isMC){
    trackInJet->SetMC(kTRUE);
    trackInJet->SetOutlierCut(1.2);
  } else
    trackInJet->SetMC(kFALSE);

  TString commonoutput = mgr->GetCommonFileName();
  commonoutput += ":EventHistsMC";
  mgr->ConnectInput(trackInJet, 1, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(trackInJet, 1, mgr->CreateContainer("HistosMC", TList::Class(), AliAnalysisManager::kOutputContainer, commonfile.Data()));
  mgr->ConnectOutput(trackInJet, 2, mgr->CreateContainer("JetTree", TTree::Class(), AliAnalysisManager::kOutputContainer, "JetTree.root"));

  return trackInJet;
}
