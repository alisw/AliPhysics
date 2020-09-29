//_____________________________________________________________________
AliAnalysisTask *AddTaskJHSInterplay(TString taskName, Bool_t ismc, Bool_t iskinematiconly, UInt_t triggSel, TString cardName, TString cardSetting,TString suffix = ""){
	// Load Custom Configuration and parameters
	// override values with parameters
    	TString combinedName = taskName;
        if(suffix.Length() > 0)
                combinedName += "_"+suffix;

	cout<<"### DEGUG Input is "<< combinedName << "\t"<< cardName <<"\t"<< cardSetting <<"\t"<<"#########"<<endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	//==== JCORRAN Efficiency TASK
	AliJHSInterplayTask *jhstask = new AliJHSInterplayTask(combinedName.Data());
	jhstask->SetDebugLevel(0);
	jhstask->SetDebugMode(0);
	jhstask->SetJCatalystTaskName("JCatalystTask");
	jhstask->SetJFJTaskName("JFJTask");
	cout << jhstask->GetName() << endl;


	// === Create AliJCORRAN ====
	AliJCard *card = new AliJCard(cardName.Data());
	card->PrintOut();
	card->ReadLine( cardSetting.Data() );
	card->ReCompile();
	card->PrintOut();
	jhstask->SetCard( card );
	jhstask->SelectCollisionCandidates( triggSel );  //Apply offline trigger selection by AliPhysicsSelectionTask

	mgr->AddTask((AliAnalysisTask*) jhstask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

	// Connect input/output
	mgr->ConnectInput(jhstask, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",jhstask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), jhstask->GetName()));
	mgr->ConnectOutput(jhstask, 1, jHist );

	return jhstask;
}

