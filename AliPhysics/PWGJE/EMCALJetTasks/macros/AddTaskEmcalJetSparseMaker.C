AliAnalysisTaskEmcalJetSparseMaker* AddTaskEmcalJetSparseMaker(
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
    ::Error("AddTaskEmcalJetSparseMaker", "No analysis manager to connect to.");
    return NULL;
  }
//=============================================================================

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEmcalJetSparseMaker", "This task requires an input event handler");
    return NULL;
  }
//=============================================================================

  const TString sTag = Form("%s_%s_%s", sUsedJets.Data(), sUsedRho.Data(), sCutType.Data());
//=============================================================================

  AliAnalysisTaskEmcalJetSparseMaker *taskEmcalJetSM = new AliAnalysisTaskEmcalJetSparseMaker(Form("AliAnalysisTaskEmcalJetSM_%s",sTag.Data()));
//taskEmcalJetSM->SetForceBeamType(0);
//taskEmcalJetSM->SetIsPythia(kTRUE);

/*taskEmcalJetSM->SetCaloTriggerPatchInfoName("EmcalTriggers");
  taskEmcalJetSM->SetTriggerTypeSel(AliAnalysisTaskEmcal::kJ2);
  taskEmcalJetSM->SetMainPatchType(AliAnalysisTaskEmcal::kTriggerLevel1Jet);*/
//=============================================================================

  AliParticleContainer *pContTrks = 0;

  if (!sUsedTrks.IsNull()) {
    pContTrks = taskEmcalJetSM->AddParticleContainer(sUsedTrks.Data());
    pContTrks->SetParticlePtCut(dTrkPtMin);
  }
//=============================================================================

  AliClusterContainer  *pContClus = 0;

  if (!sUsedClus.IsNull()) {
    pContClus = taskEmcalJetSM->AddClusterContainer(sUsedClus.Data());
    pContClus->SetClusPtCut(dCluEnMin);
  }
//=============================================================================

  AliJetContainer *pContJets = 0;

  if (!sUsedJets.IsNull()) {
    pContJets = taskEmcalJetSM->AddJetContainer(sUsedJets.Data(), sCutType.Data(), dJetR);
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

  mgr->AddTask(taskEmcalJetSM);
  mgr->ConnectInput(taskEmcalJetSM,  0, mgr->GetCommonInputContainer());

  mgr->ConnectOutput(taskEmcalJetSM, 1, mgr->CreateContainer(Form("listEmcalJetSM_%s_General",sTag.Data()),
                                                             TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectOutput(taskEmcalJetSM, 2, mgr->CreateContainer(Form("listEmcalJetSM_%s_EvH",sTag.Data()),
                                                             TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             "AnalysisResults_EvH.root"));

  mgr->ConnectOutput(taskEmcalJetSM, 3, mgr->CreateContainer(Form("listEmcalJetSM_%s_Jet",sTag.Data()),
                                                             TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             "AnalysisResults_Jet.root"));
//=============================================================================

  return taskEmcalJetSM;
}
