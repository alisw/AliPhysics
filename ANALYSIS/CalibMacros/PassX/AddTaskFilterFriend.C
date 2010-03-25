//=============================================================================
//
// *** AddTaskFilterFriend.C ***
//
// This macro initialize a complete AnalysisTask object for filtering ESD with AliAnalysisTaskFilterFriend.
//
//=============================================================================

AliAnalysisTaskFilterFriend *AddTaskFilterFriend()
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
	AliAnalysisTaskFilterFriend* filter = new AliAnalysisTaskFilterFriend("filter_1");
	mgr->AddTask(filter);

	// connecting the input/output containers
	AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput0 = mgr->GetCommonOutputContainer();

	mgr->ConnectInput (filter, 0, cinput0 );
	//mgr->ConnectOutput(filter, 0, coutput0);

	return filter;
}
