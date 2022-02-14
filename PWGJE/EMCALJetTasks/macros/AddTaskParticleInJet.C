AliEmcalTrackSelection *TrackSelectionFactory(Bool_t isAOD);

AliAnalysisTaskParticleInJet *AddTaskParticleInJet(
    Bool_t isMC,
    TString mcparticleContainer,
    TString recparticleContainer,
    TString mcjetcontainer,
    TString recjetcontainer
    )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisTaskParticleInJet *jettask = new AliAnalysisTaskParticleInJet(jettask);
  mgr->AddTask(jettask);
  jettask->SetTrackSelection(TrackSelectionFactory(mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class()));

  // Handle particle and jet containers
  AliParticleContainer *recpcont = jettask->AddParticleContainer(recparticleContainer.Data());
  recpcont->SetParticleEtaLimits(-0.5, 0.5);
  recpcont->SetParticlePhiLimits(1.4, 3.1);
  jettask->SetParticleContainerNameRec(recparticleContainer.Data());

  AliJetContainer *recjcont = jettask->AddJetContainer(recjetcontainer.Data(), AliEmcalJet::kEMCAL, 0.4);
  jettask->SetJetContainerNameRec(recjetcontainer.Data());

  if(isMC){
    AliMCParticleContainer *mcpcont = jettask->AddMCParticleContainer(mcparticleContainer.Data());
    mcpcont->SetEtaLimits(-0.5, 0.5);
    mcpcont->SetParticlePhiLimits(1.4, 3.1);
    jettask->SetParticleContainerNameMC(mcparticleContainer.Data());

    AliJetContainer *mcjcont = jettask->AddJetContainer(mcjetcontainer.Data(), AliEmcalJet::kEMCAL, 0.4);
    jettask->SetJetContainerNameMC(mcjetcontainer.Data());
  }

  TString outname = mgr->GetCommonFileName() + ":TracksInJet";
  mgr->ConnectInput(jettask, 1, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(jettask, 1, mgr->CreateContainer("TrackInJetHistos", TList::Class(), AliAnalysisManager::kOutputContainer, outname.Data()));

  return jettask;
}

AliEmcalTrackSelection *TrackSelectionFactory(Bool_t isAOD){
  AliEmcalTrackSelection *result = 0;
  if(isAOD){
    AliEmcalTrackSelectionAOD *aodcuts = new AliEmcalTrackSelectionAOD();
    aodcuts->AddFilterBit(AliAODTrack::kTrkGlobal);
    PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts *extracuts = new PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts;
    extracuts->SetMinTPCCrossedRows(120);
    result = aodcuts;
  } else {
    result = new AliEmcalTrackSelectionESD();
    AliESDtrackCuts *cuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 0);
	  cuts->SetName("Standard Track cuts");
	  cuts->SetMinNCrossedRowsTPC(120);
	  cuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    result->AddTrackCuts(cuts);
  }
  return result;
}
