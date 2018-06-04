//_____________________________________________________________________
AliAnalysisTask *AddTaskJCDijetTask(TString taskName,Bool_t isMC){
	// Load Custom Configuration and parameters
	// override values with parameters

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	//==== Set up the dijet task ====
	AliJCDijetTask *dijetTask = new AliJCDijetTask(taskName.Data(),"AOD");
	dijetTask->SetDebugLevel(5);
  	dijetTask->SetJCatalystTaskName("JCatalystTaskEP");  // AliJCatalystTask has this name hard coded
	dijetTask->SetIsMC(isMC);
	cout << dijetTask->GetName() << endl;


	mgr->AddTask((AliAnalysisTask*) dijetTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


	// Connect input/output
	mgr->ConnectInput(dijetTask, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",dijetTask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), dijetTask->GetName()));
	mgr->ConnectOutput(dijetTask, 1, jHist );

	return dijetTask;
}

