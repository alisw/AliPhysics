AliAnalysisTaskEmcalTmpSparseMaker* AddTaskEmcalTmpSparseMaker(
  const TString sUsedTrks,
  const TString sUsedClus,
  const TString sUsedJets,
  const TString sUsedRho,
  const TString sCutType,
  const Double_t dJetR,
  const Double_t dTrkPtMin   = 0.15,
  const Double_t dCluEnMin   = 0.30,
  const Double_t dJetPtMin   = 0.1,
  const Double_t dJetAreaMin = 0.,
  const Int_t kLeadingType   = 0)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTaskEmcalTmpSparseMaker", "No analysis manager to connect to.");
    return NULL;
  }
//=============================================================================

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEmcalTmpSparseMaker", "This task requires an input event handler");
    return NULL;
  }
//=============================================================================

  const TString sTag = Form("%s_%s_%s", sUsedTrks.Data(), sUsedClus.Data(), sCutType.Data());
//=============================================================================

  AliAnalysisTaskEmcalTmpSparseMaker *taskEmcalTmpSM = new AliAnalysisTaskEmcalTmpSparseMaker(Form("AliAnalysisTaskEmcalTmpSM_%s",sTag.Data()));
//taskEmcalTmpSM->SetForceBeamType(0);
//taskEmcalTmpSM->SetIsPythia(kTRUE);

/*taskEmcalTmpSM->SetCaloTriggerPatchInfoName("EmcalTriggers");
  taskEmcalTmpSM->SetTriggerTypeSel(AliAnalysisTaskEmcal::kJ2);
  taskEmcalTmpSM->SetMainPatchType(AliAnalysisTaskEmcal::kTriggerLevel1Jet);*/
//=============================================================================

  AliParticleContainer *pContTrks = 0;

  if (!sUsedTrks.IsNull()) {
    pContTrks = taskEmcalTmpSM->AddParticleContainer(sUsedTrks.Data());
    pContTrks->SetParticlePtCut(dTrkPtMin);
  }
//=============================================================================

  AliClusterContainer  *pContClus = 0;

  if (!sUsedClus.IsNull()) {
    pContClus = taskEmcalTmpSM->AddClusterContainer(sUsedClus.Data());
    pContClus->SetClusPtCut(dCluEnMin);
  }
//=============================================================================

  AliJetContainer *pContJets = 0;

  if (!sUsedJets.IsNull()) {
    pContJets = taskEmcalTmpSM->AddJetContainer(sUsedJets.Data(), sCutType.Data(), dJetR);
    pContJets->SetNameTitle(taskEmcalJetSM->GetNameJet().Data(),taskEmcalJetSM->GetNameJet().Data());
    pContJets->SetPercAreaCut(dJetAreaMin);
    pContJets->SetJetPtCut(dJetPtMin);

//  pContJets->SetLocalRhoName();
    pContJets->SetRhoName(sUsedRho.Data());

    pContJets->SetLeadingHadronType(kLeadingType);
    if (pContTrks) pContJets->ConnectParticleContainer(pContTrks);
    if (pContClus) pContJets->ConnectClusterContainer(pContClus);
  }
//=============================================================================

  mgr->AddTask(taskEmcalTmpSM);
  mgr->ConnectInput(taskEmcalTmpSM,  0, mgr->GetCommonInputContainer());

  mgr->ConnectOutput(taskEmcalTmpSM, 1, mgr->CreateContainer(Form("listEmcalTmpSM_%s_General",sTag.Data()),
                                                             TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectOutput(taskEmcalTmpSM, 2, mgr->CreateContainer(Form("listEmcalTmpSM_%s_EvH",sTag.Data()),
                                                             TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             "AnalysisResults_EvH.root"));

  if (pContTrks) mgr->ConnectOutput(taskEmcalTmpSM, 3, mgr->CreateContainer(Form("listEmcalTmpSM_%s_Trk",sTag.Data()),
                                                                            TList::Class(),
                                                                            AliAnalysisManager::kOutputContainer,
                                                                            "AnalysisResults_Trk.root"));

  if (pContClus) mgr->ConnectOutput(taskEmcalTmpSM, 4, mgr->CreateContainer(Form("listEmcalTmpSM_%s_Clu",sTag.Data()),
                                                                            TList::Class(),
                                                                            AliAnalysisManager::kOutputContainer,
                                                                            "AnalysisResults_Clu.root"));

  if (pContJets) mgr->ConnectOutput(taskEmcalTmpSM, 5, mgr->CreateContainer(Form("listEmcalTmpSM_%s_Jet",sTag.Data()),
                                                                            TList::Class(),
                                                                            AliAnalysisManager::kOutputContainer,
                                                                            "AnalysisResults_Jet.root"));
//=============================================================================

  return taskEmcalTmpSM;
}
