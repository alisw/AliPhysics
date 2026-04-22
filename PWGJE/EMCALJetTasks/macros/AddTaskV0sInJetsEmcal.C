AliAnalysisTaskV0sInJetsEmcal* AddTaskV0sInJetsEmcal(
  TString jetBranchName = "",
  Double_t dRadius = 0.2,
  TString jetBranchBgName = "",
  Double_t dRadiusBg = 0.2,
  TString outputFile = "output.root",
  Bool_t bIsMC = kFALSE,
  TString tracksName = "tracks",
  TString clustersCorrName = "CaloClustersCorr",
  TString rhoName = "Rho",
  TString sType = "TPC",
  TString label = "",
  Bool_t bEmb = kFALSE,
  TString jetMCGenBranchName = "",
  TString tracksMCGenName = "mcparticles",
  TString jetMCRecBranchName = "",
  TString tracksMCRecName = "tracks"  
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

  //embedding containers
  AliParticleContainer *trackContMCGen = 0x0;
  AliParticleContainer *trackContMCRec = 0x0;
  AliJetContainer* jetContMCGen = 0x0;
  AliJetContainer* jetContMCRec = 0x0;
  
  if(bEmb) {
    taskName += "_Emb";
    trackContMCGen = mytask->AddMCParticleContainer(tracksMCGenName);
    if(trackContMCGen) trackContMCGen->SetIsEmbedding(kTRUE);
    trackContMCRec = mytask->AddTrackContainer(tracksMCRecName);
    if(trackContMCRec) trackContMCRec->SetIsEmbedding(kTRUE);
    
    if(jetCont) jetCont->SetName("hybridLevelJets");

    jetContMCGen = mytask->AddJetContainer(jetMCGenBranchName, sType, dRadius);
    if(jetContMCGen)
    {
      jetContMCGen->SetRhoName(rhoName);
      jetContMCGen->SetLeadingHadronType(0);
      jetContMCGen->ConnectParticleContainer(trackContMCGen);
      jetContMCGen->ConnectClusterContainer(clusterCont);
      jetContMCGen->SetName("partLevelJets");
    }
    jetContMCRec = mytask->AddJetContainer(jetMCRecBranchName, sType, dRadius);
    if(jetContMCRec)
    {
      jetContMCRec->SetRhoName(rhoName);
      jetContMCRec->SetLeadingHadronType(0);
      jetContMCRec->ConnectParticleContainer(trackContMCRec);
      jetContMCRec->ConnectClusterContainer(clusterCont);
      jetContMCRec->SetName("detLevelJets");
    }
  }
  
  // Add task
  mgr->AddTask(mytask);

  // Create containers for input/output
  AliAnalysisDataContainer* cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "Std"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput2 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "Jet"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput3 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "QA"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput4 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "Cuts"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput5 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "MC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput6 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "Std_Cascade"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput7 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "Jet_Cascade"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));
  AliAnalysisDataContainer* coutput8 = mgr->CreateContainer(Form("%s_%s", containerName.Data(), "QA_Cascade"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFile.Data(), taskName.Data()));

  // Connect input/output
  mgr->ConnectInput(mytask, 0, cinput0);
  mgr->ConnectOutput(mytask, 1, coutput1);
  mgr->ConnectOutput(mytask, 2, coutput2);
  mgr->ConnectOutput(mytask, 3, coutput3);
  mgr->ConnectOutput(mytask, 4, coutput4);
  mgr->ConnectOutput(mytask, 5, coutput5);
  mgr->ConnectOutput(mytask, 6, coutput6);
  mgr->ConnectOutput(mytask, 7, coutput7);
  mgr->ConnectOutput(mytask, 8, coutput8);


  return mytask;
}
