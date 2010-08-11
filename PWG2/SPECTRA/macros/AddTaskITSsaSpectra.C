AliAnalysisTaskSEITSsaSpectra *AddTaskITSsaSpectra(Bool_t optNtuple=kFALSE){
	// Creates, configures and attaches to the train the task for pi, K , spectra
	// with ITS standalone tracks
	// Get the pointer to the existing analysis manager via the static access method.
	//==============================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		::Error("AddTaskITSsaSpectra", "No analysis manager to connect to.");
		return NULL;
	}   

	// Check the analysis type using the event handlers connected to the analysis manager.
	//==============================================================================
	if (!mgr->GetInputEventHandler()) {
		::Error("AddTaskITSsaSpectra", "This task requires an input event handler");
		return NULL;
	}   

	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	if(type.Contains("AOD")){
		::Error("AddTaskITSsaSpectra", "This task requires to run on ESD");
		return NULL;
	}

	Bool_t isMC=kFALSE;
	if (mgr->GetMCtruthEventHandler()) isMC=kTRUE;

	// Create and configure the task
	AliAnalysisTaskSEITSsaSpectra *taskits = new AliAnalysisTaskSEITSsaSpectra();
	taskits->SelectCollisionCandidates();
	taskits->SetReadMC(isMC);
	taskits->SetFillNtuple(optNtuple);
	mgr->AddTask(taskits);

	// Create ONLY the output containers for the data produced by the task.
	// Get and connect other common input/output containers via the manager as below
	//==============================================================================
	TString outputFileName = AliAnalysisManager::GetCommonFileName();
	outputFileName += ":PWG2SpectraITSsa";

	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistITSsa",
			TList::Class(),
			AliAnalysisManager::kOutputContainer,
			outputFileName );

	mgr->ConnectInput(taskits, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(taskits, 1, coutput1);
	return taskits;
}   
