//_____________________________________________________________________
AliAnalysisTask *AddTaskJEbECORR(TString taskName, Bool_t ismc, Bool_t iskinematiconly, UInt_t triggSel, TString cardName, TString cardSetting, TString ebeCentFile, Bool_t enableCORR){
	// Load Custom Configuration and parameters
	// override values with parameters

	cout<<"### DEGUG Input is "<< taskName << "\t"<< cardName <<"\t"<< cardSetting <<"\t"<<"#########"<<endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	//==== JCORRAN Efficiency TASK
	AliJEbECORRTask *jebetask = new AliJEbECORRTask(taskName.Data());
	jebetask->SetDebugLevel(0);
	jebetask->SetDebugMode(0);
	jebetask->SetEbePercentileInputFileName(ebeCentFile);
	jebetask->SetIsMC(ismc);
	jebetask->SetEnableCORR(enableCORR);
	jebetask->SetKineOnly(iskinematiconly);
	cout << jebetask->GetName() << endl;


	// === Create AliJCORRAN ====
	AliJCard *card = new AliJCard(cardName.Data());
	card->PrintOut();
	card->ReadLine( cardSetting.Data() );
	card->ReCompile();
	card->PrintOut();
	jebetask->SetCard( card );
	jebetask->SelectCollisionCandidates( triggSel );  //Apply offline trigger selection by AliPhysicsSelectionTask

	mgr->AddTask((AliAnalysisTask*) jebetask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

	// Connect input/output
	mgr->ConnectInput(jebetask, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",jebetask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), jebetask->GetName()));
	mgr->ConnectOutput(jebetask, 1, jHist );

	return jebetask;
}

