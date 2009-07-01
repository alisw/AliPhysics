//=============================================================================
//
// *** AddTaskTOFCalib.C ***
//
// This macro initialize a complete AnalysisTask object for TOF Calibration.
//
//=============================================================================

AliTOFCalibTask *AddTaskTOFCalib()
{
	// A. Get the pointer to the existing analysis manager via the static access method.
	//==============================================================================

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskTOFCalib", "No analysis manager to connect to.");
		return NULL;
	}  

	// B. Check the analysis type using the event handlers connected to the analysis
	//    manager. The availability of MC handler cann also be checked here.
	//==============================================================================

	if (!mgr->GetInputEventHandler()) {
		::Error("AddTask", "This task requires an input event handler");
		return NULL;
	}  
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

	// C. Create the task, add it to manager.
	//===========================================================================

	AliTOFCalibTask *taskTOF = new AliTOFCalibTask("TOFCalibTask");
	mgr->AddTask(taskTOF);

	// D. Configure the analysis task. Extra parameters can be used via optional
	// arguments of the AddTaskXXX() function.
	//===========================================================================
	
	// E. Create ONLY the output containers for the data produced by the task.
	// Get and connect other common input/output containers via the manager as below
	//==============================================================================
	AliAnalysisDataContainer *coutput1  = mgr->CreateContainer("histoList",  TList::Class(),
							       AliAnalysisManager::kOutputContainer,"outputHistos.root");
	AliAnalysisDataContainer *coutput2  = mgr->CreateContainer("tofArray",  TList::Class(),
							       AliAnalysisManager::kOutputContainer);

	mgr->ConnectInput(taskTOF, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(taskTOF, 1, coutput1);
	mgr->ConnectOutput(taskTOF, 2, coutput2);

	// Return task pointer at the end
	return taskTOF;
}
