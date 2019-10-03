//______________________________________________________________________________

AliAnalysisTaskCEPQA* AddTaskCEPQA()
{

	// get the manager and task
	AliAnalysisManager *aam = AliAnalysisManager::GetAnalysisManager();

  TString name = TString("CEPQA");
	AliAnalysisTaskCEPQA *task = new AliAnalysisTaskCEPQA(name.Data());

	// get input and output managers
	AliAnalysisDataContainer *aadci = aam->GetCommonInputContainer();
	AliAnalysisDataContainer *aadco1 = aam->CreateContainer
	(
		Form("CEPQAtree"),
		TTree::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:CEPQA", AliAnalysisManager::GetCommonFileName())
	);
	
	AliAnalysisDataContainer *aadco2 = aam->CreateContainer
	(
		Form("CEPQAhisto"),
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:CEPQA", AliAnalysisManager::GetCommonFileName())
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

