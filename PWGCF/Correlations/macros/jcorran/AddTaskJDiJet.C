AliJDiJetTask * AddTaskJDiJet (
		Int_t    trigger       = AliVEvent::kEMCEJE,
		TString  taskName      = "DiJetTask"   ,
		TString  jetTaskName   = "AliJJetTask" ,
		TString  cardName      = "cardAlice_pp.input",
		TString  cardSetting   = "",
		Int_t     debug         = 1
		){
	// Get the pointer to the existing analysis manager via the static access method.
	//==============================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr)
	{
		::Error("AddTaskJJet", "No analysis manager to connect to.");
		return NULL;
	}  

	// Check the analysis type using the event handlers connected to the analysis manager.
	//==============================================================================
	if (!mgr->GetInputEventHandler())
	{
		::Error("AddTaskJJet", "This task requires an input event handler");
		return NULL;
	}


	AliJCard *card = new AliJCard(cardName.Data());
	card->PrintOut();
	card->ReadLine( cardSetting.Data() );
	card->ReCompile();
	card->PrintOut();

	//-------------------------------------------------------
	// Init the task and do settings
	//-------------------------------------------------------
	AliJDiJetTask * dijetTask = new AliJDiJetTask(taskName,"AOD");
	dijetTask->SetJetTaskName( jetTaskName );
	dijetTask->SetCard( card );
	dijetTask->SelectCollisionCandidates(trigger);
	mgr->AddTask(dijetTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

	// Connect input/output
	mgr->ConnectInput(dijetTask, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",dijetTask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), dijetTask->GetName()));
	mgr->ConnectOutput(dijetTask, 1, jHist );

	return dijetTask;


}
