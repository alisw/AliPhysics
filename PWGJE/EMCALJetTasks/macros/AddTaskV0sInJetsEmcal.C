AliAnalysisTaskV0sInJetsEmcal* AddTaskV0sInJetsEmcal(
  TString jetBranchName = "",
  Double_t dRadius = 0.2,
  TString jetBranchBgName = "",
  Double_t dRadiusBg = 0.2,
  TString outputFile = "output.root",
  Bool_t bIsMC = kFALSE,
  TString tracksName = "PicoTracks",
  TString clustersCorrName = "CaloClustersCorr",
  TString rhoName = "Rho",
  TString sType = "TPC",
  TString label = ""
)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr)
  {
    Error("AddTaskV0sInJetsEmcal", "No analysis manager found.");
    return 0;
  }

  TString taskName = "V0";
  TString containerName = "V0histo";
  if(jetBranchName.Length())
  {
    taskName += Form("_%s", jetBranchName.Data());
    containerName += Form("_%s", jetBranchName.Data());
  }
  if(label.Length())
  {
    taskName += Form("_%s", label.Data());
    containerName += Form("_%s", label.Data());
  }
  AliAnalysisTaskV0sInJetsEmcal* mytask = new AliAnalysisTaskV0sInJetsEmcal(taskName.Data());
  // Configure task
  mytask->SetMCAnalysis(bIsMC);

  AliParticleContainer* trackCont = mytask->AddTrackContainer(tracksName);
  AliClusterContainer* clusterCont = mytask->AddClusterContainer(clustersCorrName);

  AliJetContainer* jetCont = mytask->AddJetContainer(jetBranchName, sType, dRadius);
  if(jetCont)
  {
    jetCont->SetRhoName(rhoName);
    jetCont->SetLeadingHadronType(0);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
  }
  AliJetContainer* jetContBg = mytask->AddJetContainer(jetBranchBgName, sType, dRadiusBg);
  if(jetContBg)
  {
    jetContBg->SetJetAreaCut(0.01);
    jetContBg->SetAreaEmcCut(0);
    jetContBg->SetJetPtCut(0);
    jetContBg->ConnectParticleContainer(trackCont);
    jetContBg->ConnectClusterContainer(clusterCont);
  }

  // Add task
  mgr->AddTask(mytask);

  // Create containers for input/output
  AliAnalysisDataContainer* cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "Std"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput2 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "QA"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput3 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "Cuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput4 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "MC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));

  // Connect input/output
  mgr->ConnectInput(mytask, 0, cinput0);
  mgr->ConnectOutput(mytask, 1, coutput1);
  mgr->ConnectOutput(mytask, 2, coutput2);
  mgr->ConnectOutput(mytask, 3, coutput3);
  mgr->ConnectOutput(mytask, 4, coutput4);

  return mytask;
}
