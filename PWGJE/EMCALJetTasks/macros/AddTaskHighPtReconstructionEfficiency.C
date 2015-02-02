AliAnalysisTask *AddTaskHighPtReconstructionEfficiency(){
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	bool isAOD = mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();

	HighPtTracks::AliHighPtReconstructionEfficiency *efficiencyTask = new HighPtTracks::AliHighPtReconstructionEfficiency("MarkusParticlesInJets");
	efficiencyTask->SetTrackCuts(CreateStandardCuts(isAOD));
	mgr->AddTask(efficiencyTask);

	TString outputcontname = mgr->GetCommonFileName();
	outputcontname += ":ParticleTree";

	efficiencyTask->ConnectInput(0, mgr->GetCommonInputContainer());
	efficiencyTask->ConnectOutput(1, mgr->CreateContainer("ParticleTree", TNtuple::Class(), AliAnalysisManager::kOutputContainer, outputcontname));


	return efficiencyTask;
}

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
