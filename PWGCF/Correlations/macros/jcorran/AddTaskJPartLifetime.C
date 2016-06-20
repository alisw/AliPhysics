//_____________________________________________________________________
AliAnalysisTask *AddTaskJPartLifetime(TString cardName, TString cardSetting){
	// Load Custom Configuration and parameters
	// override values with parameters
	// surfix in last arguments are added for subwagons.
	//
	//==== Get the pointer to the Analyis mgr
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	AliJCard *pcard = new AliJCard(cardName.Data());
	pcard->PrintOut();
	pcard->ReadLine( cardSetting.Data() );
	pcard->ReCompile();
	pcard->PrintOut();

	//==== JCORRAN TASK
	AliJPartLifetime *lttask = new AliJPartLifetime("JPartLifetime","ESD");
	lttask->SetCard(pcard);

	//==== Add task
	mgr->AddTask((AliAnalysisTask*) lttask);

	//==== Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *hist = mgr->CreateContainer(Form("%scontainer",lttask->GetName()),  TList::Class(), AliAnalysisManager::kOutputContainer,  Form("%s:%s",  AliAnalysisManager::GetCommonFileName(), lttask->GetName()));

	//==== Connect input/output
	mgr->ConnectInput(lttask, 0, cinput);
	mgr->ConnectOutput(lttask, 1, hist);
	return lttask;
}
