AliAnalysisTaskCentral* AddTaskCentral(Bool_t *simulation=kFALSE){

// Get the pointer to the existing analysis manager via the static access method.
//==============================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
    	::Error("AddTaskCentral", "No analysis manager to connect to!");
    	return NULL;
	}

// Check the analysis type using the event handlers connected to the analysis manager.
//==============================================================================
	if (!mgr->GetInputEventHandler()) {
    	::Error("AddTaskCentral", "This task requires an input event handler!");
    	return NULL;
  	}
  
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
	if (type=="AOD"){
		::Error("AddTaskCentral", "This task is not tested for AOD analysis!");
    	return NULL;
	}

	// Create and configure the task
	AliAnalysisTaskCentral *taskcentral = new AliAnalysisTaskCentral("TaskCentral");
	taskcentral->SetSimulation(kTRUE); //kTRUE if we are running on simulated data 
	mgr->AddTask(taskcentral);

	// Create ONLY the output containers for the data produced by the task.
	// Get and connect other common input/output containers via the manager as below
	//==============================================================================
	TString outputFileName = AliAnalysisManager::GetCommonFileName();
	outputFileName += ":PWG2Central";

	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cobCentral",
								  TList::Class(),
								  AliAnalysisManager::kOutputContainer,
								  outputFileName );
	
	mgr->ConnectInput(taskcentral, 0, mgr->GetCommonInputContainer()); 
	mgr->ConnectOutput(taskcentral, 0, coutput1);
	return taskcentral;
}
