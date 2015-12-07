AliAnalysisTaskCDTree *AddTaskCDTree() {

	 //--- get the current analysis manager ---//
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTask_CDTree", "No analysis manager found.");
		return 0;
	}
	// Check the analysis type using the event handlers connected to the analysis manager.
	//==============================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTask_CDTree", "This task requires an input event handler");
		return 0;
	}
	TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	Bool_t isMC;
	if(mgr->GetMCtruthEventHandler()) isMC = kTRUE;

	// Create tasks
	AliAnalysisTaskCDTree *task = new AliAnalysisTaskCDTree(inputDataType.Data());
	task->SelectCollisionCandidates(AliVEvent::kMB);
	mgr->AddTask(task);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer("output", TTree::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("output2", TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");

	// Connect input/output
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput);
	mgr->ConnectOutput(task, 2, coutput2);

	return task;
}
