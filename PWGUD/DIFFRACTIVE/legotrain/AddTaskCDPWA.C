AliAnalysisTaskCDPWA *AddTaskCDPWA() {

	 //--- get the current analysis manager ---//
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTask_CDPWA", "No analysis manager found.");
		return 0;
	}
	// Check the analysis type using the event handlers connected to the analysis manager.
	//==============================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTask_CDPWA", "This task requires an input event handler");
		return 0;
	}
	//ESD handler
	AliInputEventHandler* hdl = (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
	if (hdl) hdl->SetNeedField(kTRUE);

	TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	Bool_t isMC = kFALSE;
	if(mgr->GetMCtruthEventHandler()) isMC = kTRUE;

	// Create tasks
	AliAnalysisTaskCDPWA *task = new AliAnalysisTaskCDPWA("test");
	//task->SetRunTree(runTree);
	//task->SetRunHist(runHist);
	//task->SetIsMC(isMC);
	//task->SetRunSyst(runSyst);
	//task->SelectCollisionCandidates(AliVEvent::kMB);

	/*
	// Load other task
	gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");

	AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(isMC);

	//OADB Physics selection
	AliOADBPhysicsSelection *oadb = new AliOADBPhysicsSelection("oadb_custom");
	oadb->AddCollisionTriggerClass(AliVEvent::kUserDefined,"+CINT10-B-NOPF-ALLNOTRD","B",0);
	oadb->SetHardwareTrigger(0,"V0A || V0C");
	oadb->SetOfflineTrigger(0,"(V0A || V0C || ADA || ADC) && !V0ABG && !V0CBG && !ADABG && !ADCBG && !TPCLaserWarmUp");

	oadb->AddCollisionTriggerClass(AliVEvent::kUserDefined,"+C0SMB-B-NOPF-ALLNOTRD","B",1);
	oadb->SetHardwareTrigger(1,"SPDGFO >= 1");
	oadb->SetOfflineTrigger(1,"SPDGFO >= 1 && !V0ABG && !V0CBG && !ADABG && !ADCBG && !TPCLaserWarmUp");
	physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(oadb,0);
	*/

	// Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer("output", TTree::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("output2", TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");

	// Add task
	mgr->AddTask(task);

	// Connect input/output
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput);
	mgr->ConnectOutput(task, 2, coutput2);

	return task;
}
