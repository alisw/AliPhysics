//_____________________________________________________________________
AliAnalysisTask *AddTaskf0f2(TString taskName, Bool_t isAOD, UInt_t filterbit){
	// Load Custom Configuration and parameters
	// override values with parameters

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    if (!isAOD){
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
        AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(0);
        if(!physSelTask) { Printf("no physSelTask"); return; }
        //==== CENTRALITY
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
        AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
        //taskCentrality->SetPass(2);
    }

	AliAnalysisTaskf0f2 *task = new AliAnalysisTaskf0f2(taskName.Data(),taskName.Data());
    task->SetFilterBit(filterbit);
	mgr->AddTask((AliAnalysisTask*)task);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

	// Connect input/output
	mgr->ConnectInput(task, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",task->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), task->GetName()));
	mgr->ConnectOutput(task, 1, jHist );

	return task;
}

