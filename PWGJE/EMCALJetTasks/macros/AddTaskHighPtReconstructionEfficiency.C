AliAnalysisTask *AddTaskHighPtReconstructionEfficiency(){
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	HighPtTracks::AliHighPtReconstructionEfficiency *efficiencyTask = new HighPtTracks::AliHighPtReconstructionEfficiency("MarkusParticlesInJets");
	mgr->AddTask(efficiencyTask);

	TString outputcontname = mgr->GetCommonFileName();
	outputcontname += ":ParticleTree";

	efficiencyTask->ConnectInput(0, mgr->GetCommonInputContainer());
	efficiencyTask->ConnectOutput(1, mgr->CreateContainer("ParticleTree", TNtuple::Class(), AliAnalysisManager::kOutputContainer, outputcontname));

	return efficiencyTask;
}
