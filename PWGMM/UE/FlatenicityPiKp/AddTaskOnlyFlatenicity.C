///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskOnlyFlatenicity macro                       //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskOnlyFlatenicity *
AddTaskOnlyFlatenicity(const Char_t *taskname = "V0_Calibrated", TString detForFlat = "V0",
		bool woTrivialscaling = kFALSE, const char* suffix = "")

{
	// detForFlat: "V0", "TPC", "V0_TPC"
	// get the manager via the static access member. since it's static, you don't
	// need an instance of the class to call the function

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		return 0x0;
	}
	// get the input event handler this handler is part of the managing system and
	// feeds events to your task
	if (!mgr->GetInputEventHandler()) {
		return 0x0;
	}

	// now you create an instance of your task
	AliAnalysisTaskOnlyFlatenicity *taskFlat = new AliAnalysisTaskOnlyFlatenicity("taskFlat");
	if (!taskFlat) { return 0x0; }

	taskFlat->SetUseMC(false);
	taskFlat->SetMCclosureTest(kFALSE);
	taskFlat->SetDetectorForFlatenicity(detForFlat);
	taskFlat->SetDeltaV0(kTRUE); 				//@ Set DeltaV0 scaling
	taskFlat->SetRemoveTrivialScaling(woTrivialscaling);	//@ Trivial Nch scaling
	taskFlat->IsdEdxCalibrated(kTRUE);
	taskFlat->SetDataPeriod("18f");
	taskFlat->IsSystVarTrkCuts(kTRUE);
	taskFlat->SetSystVarTrkCuts(1);
	taskFlat->SetPtMin(0.15);
	taskFlat->SetNcl(70);
	mgr->AddTask(taskFlat);

	const char *complement;
	if (woTrivialscaling) {
		complement = "wotrivialscal";
	} else {
		complement = "wtrivialscal";
	}
	const char *complement2;
	if (detForFlat == "V0") {
		complement2 = "V0";
	}
	if (detForFlat == "TPC") {
		complement2 = "TPC";
	}
	if (detForFlat == "V0_TPC") {
		complement2 = "V0_TPC";
	}

	mgr->ConnectInput(taskFlat, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(
			taskFlat, 1,
			mgr->CreateContainer(
				Form("cList%s_%s_%s_%s", taskname, complement, complement2, suffix),
				TList::Class(), AliAnalysisManager::kOutputContainer,
				Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname)));

	return taskFlat;
}
