/**
 * \brief Creating reduced jet tree creator task
 *
 * Creates a reduced jet tree creator task, initialises track cuts and adds it to the
 * analysis manager.
 *
 * \return The tree creator task
 */
AliAnalysisTask *AddTaskHighPtReconstructionEfficiency(){
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	bool isAOD = mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();

	HighPtTracks::AliHighPtReconstructionEfficiency *efficiencyTask = new HighPtTracks::AliHighPtReconstructionEfficiency("MarkusParticlesInJets");
	efficiencyTask->SetStandardTrackCuts(CreateStandardCuts(isAOD));
	efficiencyTask->SetHybridTrackCuts(CreateHybridTrackCuts(isAOD));
	mgr->AddTask(efficiencyTask);

	TString outputcontname = mgr->GetCommonFileName();
	outputcontname += ":ParticleTree";

	efficiencyTask->ConnectInput(0, mgr->GetCommonInputContainer());
	efficiencyTask->ConnectOutput(1, mgr->CreateContainer("ParticleTree", TNtuple::Class(), AliAnalysisManager::kOutputContainer, outputcontname));


	return efficiencyTask;
}

/**
 * \brief Creating standard track cuts
 *
 * Creating standard track cuts as used in the \f$ R_{AA} \f$ analysis
 *
 * \param isAOD
 * \return Standard track cuts
 */
AliESDtrackCuts *CreateStandardCuts(Bool_t isAOD){
	// Create std track cuts
	AliESDtrackCuts *trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
	trackCuts->SetName("Standard Track cuts");
	trackCuts->SetMinNCrossedRowsTPC(120);
	trackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
	if(isAOD){
		// Switch off the golden cut
		trackCuts->SetMaxChi2TPCConstrainedGlobal(1e10);
	}
	return trackCuts;
}

/**
 * \brief Creating hybrid track cuts
 *
 * Creating hybrid track cuts
 *
 * \param isAOD
 * \return Hybrid track cuts
 */
AliESDtrackCuts *CreateHybridTrackCuts(Bool_t isAOD){
  AliESDtrackCuts* hybridTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
  hybridTrackCuts->SetName("Global Hybrid tracks, loose DCA");
  hybridTrackCuts->SetMaxDCAToVertexXY(2.4);
  hybridTrackCuts->SetMaxDCAToVertexZ(3.2);
  hybridTrackCuts->SetDCAToVertex2D(kTRUE);
  hybridTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  hybridTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  if(isAOD){
    // Switch off the golden cut
    hybridTrackCuts->SetMaxChi2TPCConstrainedGlobal(1e10);
  }

  return hybridTrackCuts;
}


