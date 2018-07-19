//_____________________________________________________________________
AliAnalysisTask *AddTaskJCIaaEPMulti(TString taskName, int EPdetID, TString cardName, TString jtrigg, TString jassoc, TString cardSetting, TString inclusFileName=""){
	// Load Custom Configuration and parameters
	// override values with parameters

	cout<<"### DEBUG Input is "<< cardName <<"\t"<<jtrigg<<"\t"<<jassoc<<"\t"<<inclusFileName<<"\t"<<"#########"<<endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	const int NEPbins = 6;
	double EPbins[NEPbins+1] = {0,15,30,45,60,75,90};

	//==== Set up di-hadron correlation jT task ====
	AliJCIaaEPTask *myTask[NEPbins];
	for(int i=0;i<NEPbins;i++) {
		TString mytaskName = Form("%s_EP%.0f_%.0f",taskName.Data(),EPbins[i],EPbins[i+1]);
		myTask[i] = new AliJCIaaEPTask(mytaskName.Data(),"JOD");
		myTask[i]->SetDebugLevel(5);
		myTask[i]->SetJFlowBaseTaskName("JFlowBaseTask");  // AliJFlowBaseTask has this name hard coded
		myTask[i]->SetEPDector( EPdetID );
		myTask[i]->SetEPmin(EPbins[i]);
		myTask[i]->SetEPmax(EPbins[i+1]);
		cout << myTask[i]->GetName() << endl;
	}

	// === Set up JCard ====
	AliJCard *card = new AliJCard(cardName.Data());
	card->PrintOut();
	card->ReadLine( cardSetting.Data() );
	card->ReCompile();
	card->PrintOut();

	// === Create analysis object ===

	AliJIaaAna *fAna[NEPbins];
	for(int i=0;i<NEPbins;i++) {
		fAna[i] = new AliJIaaAna( kFALSE );
		fAna[i]->SetCard( card );
		fAna[i]->SetTrigger( jtrigg.Data() );
		fAna[i]->SetAssoc( jassoc.Data() );
		fAna[i]->SetEnableEP( true );
		if( inclusFileName ) fAna[i]->SetInclusiveFile(inclusFileName.Data());
	}

	for(int i=0;i<NEPbins;i++) {
		myTask[i]->SetAnalysis( fAna[i] );
		mgr->AddTask((AliAnalysisTask*) myTask[i]);
	}


	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


	// Connect input/output
	for(int i=0;i<NEPbins;i++) {
		mgr->ConnectInput(myTask[i], 0, cinput);
		AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",myTask[i]->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), myTask[i]->GetName()));
		mgr->ConnectOutput(myTask[i], 1, jHist );
	}

	return myTask[0];
}

