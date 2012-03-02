AliAnalysisTaskDiHadronPID *AddTaskDiHadronPID() {
    
	// AddTask Macro (v 8.00).
	// Updated: Mar 2nd. 2012.
	
	// Get the current analysis manager.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskDiHadronPID.C", "No analysis manager found.");
        return 0;
    }
    
    // Add the task to the manager.
    AliAnalysisTaskDiHadronPID *task = new AliAnalysisTaskDiHadronPID("DiHadronPID");
    
    task->SetVerbose(kFALSE);
    task->SetCalculateMixedEvents();
	
	mgr->AddTask(task);
    
	// Data containers.
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
	mgr->ConnectInput(task, 0, cinput); 
	
	AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("DiHadronPID", TList::Class(),
                         AliAnalysisManager::kOutputContainer,"DiHadronPID.root");
	
	mgr->ConnectOutput (task,  1, coutput1);
	
	return task;
	
}


