//_____________________________________________________________________
AliAnalysisTask *AddTaskJSound(TString taskName,Bool_t isMC){
	// Load Custom Configuration and parameters
	// override values with parameters

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	//==== Set up di-hadron correlation jT task ====
	AliJSoundsTask *sTask = new AliJSoundsTask(taskName.Data(),"AOD");
	sTask->SetDebugLevel(5);
  	sTask->SetJCatalystTaskName("JCatalystTask");  // AliJCatalystTask has this name hard coded
	sTask->SetIsMC(isMC);
	cout << sTask->GetName() << endl;


	mgr->AddTask((AliAnalysisTask*) sTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


	// Connect input/output
	mgr->ConnectInput(sTask, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",sTask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), sTask->GetName()));
	mgr->ConnectOutput(sTask, 1, jHist );

	return sTask;
}

