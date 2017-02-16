//______________________________________________________________________________

AliAnalysisTaskCEPAnalysis* AddTaskCEPAnalysis()
{

	// get the manager and task
	AliAnalysisManager *aam = AliAnalysisManager::GetAnalysisManager();

  TString name = TString("CEPAnalysis");
	AliAnalysisTaskCEP *task = new AliAnalysisTaskCEP(
    name.Data(),taskConfig, numTracksMax);

	// get input and output managers
	AliAnalysisDataContainer *aadci = aam->GetCommonInputContainer();
	AliAnalysisDataContainer *aadco1 = aam->CreateContainer
	(
		Form("CEPHist"),
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:CEPHist", AliAnalysisManager::GetCommonFileName())
	);

	AliAnalysisDataContainer *aadco2 = aam->CreateContainer
	(
		Form("CEP"),
		TTree::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:CEP", AliAnalysisManager::GetCommonFileName())
	);
	
	// add task and connect input and output managers
	aam->AddTask(task);
	aam->ConnectInput (task,0,aadci);
	aam->ConnectOutput(task,1,aadco1);
	aam->ConnectOutput(task,2,aadco2);

	// return pointer to Task
	return task;

}

//______________________________________________________________________________

