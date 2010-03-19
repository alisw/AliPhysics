//=============================================================================
//
// *** AddTaskAddObject.C ***
//
// This macro initialize a complete AnalysisTask object for filtering ESD with AliAnalysisTaskFilterFriendSecond.
//
//=============================================================================

AliAnalysisTaskAddObject *AddTaskAddObject()
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
	AliAnalysisTaskAddObject* add = new AliAnalysisTaskAddObject("addObj");
	mgr->AddTask(add);

	// connecting the input/output containers
	AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("histo",  TH1::Class(), AliAnalysisManager::kOutputContainer, "AliESDfriends_v1.root");

	mgr->ConnectInput (add, 0, cinput0 );
	mgr->ConnectOutput(add, 0, coutput0);

	return add;
}
