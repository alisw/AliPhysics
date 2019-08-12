// $Id$

AliJJetJtTask* AddTaskJJetJt(
		Int_t       trigger            = AliVEvent::kEMCEJE,
		TString  taskName      = "AliJJetJtTask"   ,
		TString  jetTaskName   = "AliJJetTask" ,
		TString  cardName      = "card.input",
		TString  cardSetting   = "",
		Int_t	    Nrandom    =  1,
		Int_t	    moveJet =  1,
		Int_t	    doMC       = 0,
    Int_t     leading = 0,
		Int_t       debug 		 = 1,
        double     maxDeltaR = 0.5
		)
{  

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

	//-------------------------------------------------------
	// Init the task and do settingTasks
	//-------------------------------------------------------
	cout<<"card_name input"<<endl;
	AliJCard *card = new AliJCard(cardName.Data());
	card->PrintOut();
	card->ReadLine( cardSetting.Data() );
	card->ReCompile();
	card->PrintOut();

	AliJJetJtTask * jtTask = new AliJJetJtTask(taskName,"AOD");
	jtTask->SetJetTaskName(jetTaskName);
	jtTask->SetMC(doMC);
	jtTask->SetCard( card );
	jtTask->SelectCollisionCandidates(trigger);
	jtTask->SetNrandom(Nrandom);
	jtTask->SetMoveJet(moveJet);
    jtTask->SetLeadingJets(leading);
    jtTask->SetMaxDeltaRCorr(maxDeltaR);
	mgr->AddTask(jtTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


	// Connect input/output
	mgr->ConnectInput(jtTask, 0, cinput);
	AliAnalysisDataContainer *jjtHist = mgr->CreateContainer(Form("%scontainer",jtTask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), jtTask->GetName()));
	mgr->ConnectOutput(jtTask, 1, jjtHist );


	return jtTask;
}
