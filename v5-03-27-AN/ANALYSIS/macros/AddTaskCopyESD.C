//=============================================================================
//
// *** AddTaskCopyESD.C ***
//
// This macro initialize a complete AnalysisTask object for Copying ESD.
//
//=============================================================================

AliAnalysisTaskCopyESD *AddTaskCopyESD()
{

	// pointer to the analysis manager
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskTOFCalib", "No analysis manager to connect to.");
		return NULL;
	}  

	// check the input handler
	if (!mgr->GetInputEventHandler()) {
		::Error("AddTask", "This task requires an input event handler");
		return NULL;
	}  
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

	// create the task
	AliAnalysisTaskCopyESD *copy = new AliAnalysisTaskCopyESD("ESD copying task");
	mgr->AddTask(copy);        

	// connecting the input/output containers
	AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput0 = mgr->GetCommonOutputContainer();

	mgr->ConnectInput (copy, 0, cinput0 );
	mgr->ConnectOutput(copy, 0, coutput0);

	return copy;
}
