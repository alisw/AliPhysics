///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskFlatenicityPiKp macro                       //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskMCCorrections *
AddTaskMCCorrections(const Char_t *taskname = "V0_Calibrated", TString detForFlat = "V0",
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
	AliAnalysisTaskMCCorrections *taskFlat =
		new AliAnalysisTaskMCCorrections("taskFlat");
	if (!taskFlat)
		return 0x0;

	taskFlat->SetUseMC(true);
	taskFlat->SetDeltaV0(kTRUE); 				//@ Set DeltaV0 scaling
	taskFlat->SetRemoveTrivialScaling(woTrivialscaling);	//@ Trivial Nch scaling
	taskFlat->SetDataPeriod("16k");
	taskFlat->SetSystVarTrkCuts(1);
	taskFlat->SetPtMin(0.15);
	taskFlat->SetNcl(70);
	mgr->AddTask(taskFlat);

	mgr->ConnectInput(taskFlat, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(
			taskFlat, 1,
			mgr->CreateContainer(
				Form("cList%s_%s", taskname,  suffix),
				TList::Class(), AliAnalysisManager::kOutputContainer,
				Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname)));

	return taskFlat;
}
