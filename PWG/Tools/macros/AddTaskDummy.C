AliAnalysisTaskDummy *AddTaskDummy(const char *name, Int_t timeout) {
  /// Add dummy task for read check
  /// Switch off FMD branch for a test

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskDummy", "No analysis manager found.");
		return nullptr;
	}

	
	// AliESDInputHandler *inh = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	// inh->SetInactiveBranches("AliESDFMD"); // Disable FMD branch
  
// 	AliAnalysisTaskDummy *task=NULL;
 	AliAnalysisTaskDummy *task = new AliAnalysisTaskDummy(name, timeout);
 	mgr->AddTask(task);
 	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
// 	
// 	mgr->AddTask(task);
 	return task;
}
