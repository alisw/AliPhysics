AliJJetCORRTask * AddTaskJJetCORR (
		Int_t       trigger            = AliVEvent::kEMCEJE,
		TString  taskName      = "JJetCORRTask"   ,
		TString  jetTaskName   = "AliJJetTask" ,
		TString  cardName      = "cardAlice_pp.input",
		TString  cardSetting   = "",
		int      TargetJetIndex  = 1,
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
	AliJJetCORRTask * jjetcorrTask = new AliJJetCORRTask(taskName,"AOD");
	jjetcorrTask->SetJetTaskName( jetTaskName );
	jjetcorrTask->SetCard( card );
	jjetcorrTask->SetTargetJetIndex( TargetJetIndex );
	jjetcorrTask->SelectCollisionCandidates(trigger);
	jjetcorrTask->SetDebugMode( debug );
	mgr->AddTask(jjetcorrTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

	// Connect input/output
	mgr->ConnectInput(jjetcorrTask, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",jjetcorrTask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), jjetcorrTask->GetName()));
	mgr->ConnectOutput(jjetcorrTask, 1, jHist );

	return jjetcorrTask;


}
