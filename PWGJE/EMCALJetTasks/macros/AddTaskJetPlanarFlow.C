AliAnalysisTaskJetPlanarFlow* AddTaskJetPlanarFlow(
                    const char * njets1,
                    const char * ntracks1,
								    const char * njets2,
                    const char * ntracks2,
								    const char * njets3,
                    const char * ntracks3, 
								    const Double_t R,
								    const char * nrhoBase, 
								    const char *type,				      
								    const char *CentEst,
								    Int_t       pSel,
								    AliAnalysisTaskJetPlanarFlow::JetAnalysisType jetAnalysisType = AliAnalysisTaskJetPlanarFlow::kData,
								    AliAnalysisTaskJetPlanarFlow::JetSubType jetSubType = AliAnalysisTaskJetPlanarFlow::kNoSub) {
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskJetPlanarFlow", "No analysis manager found.");
    return 0;
  }
  Bool_t ismc = kFALSE;
  ismc = (mgr->GetMCtruthEventHandler()) ? kTRUE : kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AliAnalysisTaskJetPlanarFlow", "This task requires an input event handler");
    return NULL;
  }

  TString wagonName1, wagonName2, wagonName3;

  wagonName1 = Form("AliAnalysisTaskJetPlanarFlow_%s_Histograms", njets1);
  wagonName2 = Form("AliAnalysisTaskJetPlanarFlow_%s_JetTree", njets1);
  wagonName3 = Form("AliAnalysisTaskJetPlanarFlow_%s_ParticlTree", njets1);

  // Configure jet tagger task
  AliAnalysisTaskJetPlanarFlow *task = new AliAnalysisTaskJetPlanarFlow(wagonName1.Data());

  //  task->SetNCentBins(4);
  task->SetJetAnalysisType(jetAnalysisType);
  task->SetJetSubType(jetSubType);
  task->SetJetRadius(R);

  // TString thename(njetsBase);
  // if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
  // task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont1 = 0x0;
  AliParticleContainer *trackCont2 = 0x0;
  AliParticleContainer *trackCont3 = 0x0;

  if (jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kData || jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kDet) {
    trackCont1 = task->AddTrackContainer(ntracks1);
  }
  if (jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kTrue) {
    trackCont1 = task->AddMCParticleContainer(ntracks1);
  }
  if (jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kTrueDet) {
    trackCont1 = task->AddTrackContainer(ntracks1);
    trackCont2 = task->AddMCParticleContainer(ntracks2);
  }
  if (jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kEmb) {
    trackCont1 = task->AddParticleContainer(ntracks1);
    trackCont2 = task->AddParticleContainer(ntracks2);
    trackCont3 = task->AddMCParticleContainer(ntracks3);
  }

  AliJetContainer *JetCont1 = 0x0;
  AliJetContainer *JetCont2 = 0x0;
  AliJetContainer *JetCont3 = 0x0;

  AliClusterContainer *clusterCont = task->AddClusterContainer("");
  TString strType(type);

  JetCont1 = task->AddJetContainer(njets1, strType, R);
  if (JetCont1) {
    JetCont1->SetRhoName(nrhoBase);
    JetCont1->ConnectParticleContainer(trackCont1);
    JetCont1->ConnectClusterContainer(clusterCont);
    JetCont1->SetPercAreaCut(0.6);
    JetCont1->SetJetRadius(R);
    JetCont1->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
  }

  if (jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kTrueDet || jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kEmb) {
    JetCont2 = task->AddJetContainer(njets2, strType, R);
    if (JetCont2) {
      JetCont2->SetRhoName(nrhoBase);
      JetCont2->ConnectParticleContainer(trackCont2);
      JetCont2->ConnectClusterContainer(clusterCont);
      JetCont2->SetPercAreaCut(0.6);
      JetCont2->SetJetRadius(R);
      JetCont2->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
    }
  }

  if (jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kEmb) {
    JetCont3 = task->AddJetContainer(njets3, strType, R);
    if (JetCont3) {
      JetCont3->SetRhoName(nrhoBase);
      JetCont3->ConnectParticleContainer(trackCont3);
      JetCont3->ConnectClusterContainer(clusterCont);
      JetCont3->SetPercAreaCut(0.6);
      JetCont3->SetJetRadius(R);
      JetCont3->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
    }
  }

  task->SetCaloTriggerPatchInfoName("");
  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);
  task->SetUseAliAnaUtils(kFALSE);

  mgr->AddTask(task);

  // Connnect input
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // Connect output
  TString contName1(wagonName1);
  TString contName2(wagonName2);
  TString contName3(wagonName3);

  if (jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kData) {
    contName1 += "_Data";
    contName2 += "_Data";
    contName3 += "_Data";
  }
  if (jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kDet) {
    contName1 += "_Det";
    contName2 += "_Det";
    contName3 += "_Det";
  }
  if (jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kTrue) {
    contName1 += "_True";
    contName2 += "_True";
    contName3 += "_True";
  }
  if (jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kTrueDet) {
    contName1 += "_TrueDet";
    contName2 += "_TrueDet";
    contName3 += "_TrueDet";
  }
  if (jetAnalysisType == AliAnalysisTaskJetPlanarFlow::kEmb) {
    contName1 += "_Emb";
    contName2 += "_Emb";
    contName3 += "_Emb";
  }

  if (jetSubType == AliAnalysisTaskJetPlanarFlow::kNoSub) {
    contName1 += "_NoSub";
    contName2 += "_NoSub";
    contName3 += "_NoSub";
  }
  if (jetSubType == AliAnalysisTaskJetPlanarFlow::kAreaSub) {
    contName1 += "_AreaSub";
    contName2 += "_AreaSub";
    contName3 += "_AreaSub";
  }

  TString outputfile = Form("%s", AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 1, coutput1);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 2, coutput2);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(contName3.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 3, coutput3);

  return task;
}
