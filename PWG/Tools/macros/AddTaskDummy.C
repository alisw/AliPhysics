void AddTaskDummy() {
  /// Add dummy task for read check
  /// Switch off FMD branch for a test

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskDummy", "No analysis manager found.");
		return ;
	}

	
	AliESDInputHandler *inh = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	inh->SetInactiveBranches("AliESDFMD"); // Disable FMD branch
  
// 	AliAnalysisTaskDummy *task=NULL;
// 	task= new AliAnalysisTaskDummy();
// 	
// 	mgr->AddTask(task);
}
