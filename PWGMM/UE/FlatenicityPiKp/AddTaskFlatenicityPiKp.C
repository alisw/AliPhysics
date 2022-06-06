///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskFlatenicityPiKp macro                           //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskFlatenicityPiKp *
AddTaskFlatenicityPiKp(const Char_t *taskname = "Flat", TString detForFlat = "V0",
		Bool_t woTrivialscaling = kFALSE, Bool_t useMC = kTRUE,
		Bool_t performMCclosuretest = kFALSE, TString V0Mbin = "0_1")

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
	AliAnalysisTaskFlatenicityPiKp *taskFlat =
		new AliAnalysisTaskFlatenicityPiKp("taskFlat");
	if (!taskFlat)
		return 0x0;
	taskFlat->SetUseMC(useMC);
	taskFlat->SetDetectorForFlatenicity(detForFlat);
	taskFlat->IsdEdxCalibrated(kTRUE);
	taskFlat->SetDataPeriod("16k");
	taskFlat->SetV0MBin(V0Mbin);
	taskFlat->SetMCclosureTest(performMCclosuretest);
	taskFlat->SetPtMin(0.15);
	taskFlat->SetNcl(70);
	taskFlat->SetRemoveTrivialScaling(woTrivialscaling);
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

	// complement += detForFlat;

	mgr->ConnectInput(taskFlat, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(
			taskFlat, 1,
			mgr->CreateContainer(
				Form("cList%s_%s_%s_V0MPerc_%s", taskname, complement, complement2,V0Mbin.Data()),
				TList::Class(), AliAnalysisManager::kOutputContainer,
				Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname)));

	return taskFlat;
}
