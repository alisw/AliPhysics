AliJJetQATask * AddTaskJJetQA (
		Int_t       trigger            = AliVEvent::kAny,
		TString  taskName      = "JJetQATask"   ,
		TString  jetTaskName   = "AliJJetTask" ,
		UInt_t flags = 0,
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
		::Error("AddTaskJJetQA", "No analysis manager to connect to.");
		return NULL;
	}  

	// Check the analysis type using the event handlers connected to the analysis manager.
	//==============================================================================
	if (!mgr->GetInputEventHandler())
	{
		::Error("AddTaskJJetQA", "This task requires an input event handler");
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
	AliJJetQATask * jjetqaTask = new AliJJetQATask(taskName,"AOD");
	jjetqaTask->SetJetTaskName( jetTaskName );
	jjetqaTask->AddFlags(flags);
	jjetqaTask->SetCard( card );
	jjetqaTask->SetTargetJetIndex( TargetJetIndex );
	jjetqaTask->SelectCollisionCandidates(trigger);
	jjetqaTask->SetDebugMode( debug );
	mgr->AddTask(jjetqaTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

	// Connect input/output
	mgr->ConnectInput(jjetqaTask, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",jjetqaTask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), jjetqaTask->GetName()));
	mgr->ConnectOutput(jjetqaTask, 1, jHist );

	return jjetqaTask;


}
