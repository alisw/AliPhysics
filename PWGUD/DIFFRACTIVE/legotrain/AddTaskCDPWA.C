AliAnalysisTaskCDPWA *AddTaskCDPWA(
		Bool_t IsRun2 = kTRUE,//Select Run1/Run2
		Bool_t IsSaveGap = kFALSE,//Store Gap events
		Bool_t IsComb = kTRUE,//Store combinatorics
		Bool_t IsSaveGen = kFALSE,//For DIME/DRgen and PWA
		Bool_t IsPythia8 = kFALSE//For Pythia8
		)
{

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
	// Event Handler
	AliInputEventHandler* hdl = (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
	if (hdl) hdl->SetNeedField(kTRUE);

	TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	Bool_t isMC = kFALSE;
	if(mgr->GetMCtruthEventHandler()) isMC = kTRUE;

	// Create tasks
	AliAnalysisTaskCDPWA *task = new AliAnalysisTaskCDPWA("test");
	task->SetIsRun2(IsRun2);
	task->SetIsMC(isMC);
	task->SetSaveGapEvents(IsSaveGap);
	task->SetCombinatoricsMode(IsComb);
	task->SetSaveGenParticle(IsSaveGen);
	task->SetIsPythia8(IsPythia8);

	/* For the local-test
	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
	AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(isMC);

	// Load other task
	gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");

	AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(isMC);

	//OADB Physics selection
	if (IsRun2) {
		AliOADBPhysicsSelection *oadb = new AliOADBPhysicsSelection("oadb_custom");
		oadb->AddCollisionTriggerClass(AliVEvent::kUserDefined,"+CINT10-B-NOPF-ALLNOTRD","B",0);
		oadb->SetHardwareTrigger(0,"V0A || V0C || ADA || ADC");
		oadb->SetOfflineTrigger(0,"(V0A || V0C || ADA || ADC) && !V0ABG && !V0CBG && !ADABG && !ADCBG && !TPCLaserWarmUp");

		oadb->AddCollisionTriggerClass(AliVEvent::kUserDefined,"+C0SMB-B-NOPF-ALLNOTRD","B",1);
		oadb->SetHardwareTrigger(1,"SPDGFO >= 1");
		oadb->SetOfflineTrigger(1,"SPDGFO >= 1 && !V0ABG && !V0CBG && !ADABG && !ADCBG && !TPCLaserWarmUp");
		physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(oadb,0);
	}
	else {
		task->SelectCollisionCandidates(AliVEvent::kMB);
	}
	*/

	// Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer("output", TTree::Class(), AliAnalysisManager::kOutputContainer,
			Form("%s:PWATask", AliAnalysisManager::GetCommonFileName()));
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("output2", TTree::Class(), AliAnalysisManager::kOutputContainer,
			Form("%s:PWATask", AliAnalysisManager::GetCommonFileName()));
	AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("output3", TList::Class(), AliAnalysisManager::kOutputContainer,
			Form("%s:PWATask", AliAnalysisManager::GetCommonFileName()));

	// Add task
	mgr->AddTask(task);

	// Connect input/output
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput);
	mgr->ConnectOutput(task, 2, coutput2);
	mgr->ConnectOutput(task, 3, coutput3);

	return task;
}
