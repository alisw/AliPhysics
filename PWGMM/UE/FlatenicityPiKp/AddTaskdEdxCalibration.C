///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskdEdxCalibration macro                       //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskdEdxCalibration *
AddTaskdEdxCalibration(const char *taskname = "dEdxCalibration", const bool isdEdxCal = false, const char* suffix = "")

{
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
	AliAnalysisTaskdEdxCalibration *taskFlat =
		new AliAnalysisTaskdEdxCalibration("taskFlat");
	if (!taskFlat)
		return 0x0;

	taskFlat->SetUseMC(true);
	taskFlat->IsdEdxCalibrated(isdEdxCal);
	taskFlat->SetDataPeriod("16l");
	taskFlat->SetPtMin(0.15);
	taskFlat->SetNcl(70);
	mgr->AddTask(taskFlat);

	const char* complement = "";
	if (isdEdxCal) { complement = "dEdxCal"; }
	else { complement = "dEdxNotCal"; }

	mgr->ConnectInput(taskFlat, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(
			taskFlat, 1,
			mgr->CreateContainer(
				Form("cList%s_%s_%s",taskname,complement,suffix),
				TList::Class(), AliAnalysisManager::kOutputContainer,
				Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname)));

	return taskFlat;
}
