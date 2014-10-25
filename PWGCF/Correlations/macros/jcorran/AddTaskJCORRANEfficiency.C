//_____________________________________________________________________
AliAnalysisTask *AddTaskJCORRANEfficiency(TString taskName, int fTriggerMask){
    // Load Custom Configuration and parameters
    // override values with parameters

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    //==== JCORRAN Efficiency TASK
    AliJEfficiencyTask *jefftask = new AliJEfficiencyTask(taskName.Data(),"JOD");
    jefftask->SetDebugLevel(0);
    jefftask->SetFilterTaskName("PWGCFJCORRANTask");

    AliJEfficiencyScanner *fEffScanner;
    fEffScanner = new AliJEfficiencyScanner("EfficiencyScanner");
    fEffScanner->SetMBTriggMask( fTriggerMask );
    jefftask->SetJEfficiencyScanner( fEffScanner );

    mgr->AddTask((AliAnalysisTask*) jefftask);


    // Create containers for input/output
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

    // Connect input/output
	mgr->ConnectInput(jefftask, 0, cinput);
	// Connect input/output
	AliAnalysisDataContainer *effHist = mgr->CreateContainer(Form("%scontainer",jefftask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), jefftask->GetName()));
	mgr->ConnectOutput(jefftask, 1, effHist );

	return jefftask;
}

