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
	//ESD handler
	AliInputEventHandler *hdl = (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
	if(hdl) hdl->SetNeedField(kTRUE);

	TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	Bool_t isMC = kFALSE;
	if(mgr->GetMCtruthEventHandler()) isMC = kTRUE;

	// Create tasks
	AliAnalysisTaskCDTree *task = new AliAnalysisTaskCDTree("test");
	//task->SetRunTree(runTree);
	//task->SetRunHist(runHist);
	//task->SetIsMC(isMC);
	//task->SetRunSyst(runSyst);
	//task->SelectCollisionCandidates(AliVEvent::kMB);
	//mgr->AddTask(task);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer("output", TTree::Class(), AliAnalysisManager::kOutputContainer,
			Form("%s:PWATask_7TeV",AliAnalysisManager::GetCommonFileName()));
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("output2", TList::Class(), AliAnalysisManager::kOutputContainer,
			Form("%s:PWATask_7TeV",AliAnalysisManager::GetCommonFileName()));

	mgr->AddTask(task);
	mgr->SetDebugLevel(0);

	// Connect input/output
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput);
	mgr->ConnectOutput(task, 2, coutput2);

	return task;
}
