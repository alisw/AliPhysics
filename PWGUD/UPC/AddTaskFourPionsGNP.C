AliAnalysisTaskUpcFourPions* AddTaskFourPions(Bool_t runTree = kFALSE, Bool_t runHist = kTRUE, Bool_t runSyst = kFALSE, Int_t tracking = 0) {

	//--- get the current analysis manager ---//
	AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTask_UpcFourPionsGNP", "No analysis manager found.");
		return 0;
	}

	// Check the analysis type using the event handlers connected to the analysis manager.
	//==============================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTask_UpcFourPionsGNP", "This task requires an input event handler");
		return 0;
	}

	TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

	// Create tasks
	AliAnalysisTaskUpcEtaCAWP* task = new AliAnalysisTaskUpcEtaCAWP(inputDataType.Data());
	task->SetRunTree(runTree);
	task->SetRunHist(runHist);
	task->SetRunSyst(runSyst);
	task->SetTracking(tracking);
	mgr->AddTask(task);

	// Create containers for input/output
	AliAnalysisDataContainer* cinput1 = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer* coutput1 = mgr->CreateContainer("TriggerList", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));
	AliAnalysisDataContainer* coutput2 = mgr->CreateContainer("HistogramList", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));
	AliAnalysisDataContainer* coutput3 = mgr->CreateContainer("KstarList", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));
	AliAnalysisDataContainer* coutput4 = mgr->CreateContainer("FourPionList", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));
	AliAnalysisDataContainer* coutput5 = mgr->CreateContainer("K0s3PiPi4K2K4PiList", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));

	// Connect input/output
	mgr->ConnectInput(task, 0, cinput1);
	mgr->ConnectOutput(task, 1, coutput1);
	mgr->ConnectOutput(task, 2, coutput2);
	mgr->ConnectOutput(task, 3, coutput3);
	mgr->ConnectOutput(task, 4, coutput4);
	mgr->ConnectOutput(task, 5, coutput5);

	return task;
}