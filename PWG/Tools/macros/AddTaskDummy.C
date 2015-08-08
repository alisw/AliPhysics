void AddTaskDummy() {
  /// Add dummy task for read check
  /// Switch off FMD branch for a test

	AliESDInputHandler *inh = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	inh->SetInactiveBranches("AliESDFMD"); // Disable FMD branch
  
}
