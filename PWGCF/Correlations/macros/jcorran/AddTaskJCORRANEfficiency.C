//_____________________________________________________________________
AliAnalysisTask *AddTaskJCORRANEfficiency(){
    // Load Custom Configuration and parameters
    // override values with parameters

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    //==== JCORRAN Efficiency TASK
    AliJEfficiencyTask *jefftask = new AliJEfficiencyTask("JCORRANEfficiencyTask","JOD");
    jefftask->SetDebugLevel(0);
    jefftask->SetFilterTaskName("PWGCFJCORRANTask");

    AliJEfficiencyScanner *fEffScanner;
    fEffScanner = new AliJEfficiencyScanner("EfficiencyScanner");
    jefftask->SetJEfficiencyScanner( fEffScanner );

    mgr->AddTask((AliAnalysisTask*) jefftask);


    // Create containers for input/output
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

    // Connect input/output
	mgr->ConnectInput(jefftask, 0, cinput);
	// Connect input/output
	AliAnalysisDataContainer *effHist = mgr->CreateContainer("JEffHist",  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:JEffHist",AliAnalysisManager::GetCommonFileName()));
	mgr->ConnectOutput(jefftask, 1, effHist );

	return jefftask;
}

