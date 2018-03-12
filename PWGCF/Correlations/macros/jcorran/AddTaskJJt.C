//_____________________________________________________________________
AliAnalysisTask *AddTaskJJt(TString taskName, TString cardName, TString jtrigg, TString jassoc, TString cardSetting, TString inclusFileName=""){
	// Load Custom Configuration and parameters
	// override values with parameters

	cout<<"### DEBUG Input is "<< cardName <<"\t"<<jtrigg<<"\t"<<jassoc<<"\t"<<inclusFileName<<"\t"<<"#########"<<endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	//==== Set up di-hadron correlation jT task ====
	AliJJtTask *jtTask = new AliJJtTask(taskName.Data(),"JOD");
	jtTask->SetDebugLevel(5);
  	jtTask->SetJCatalystTaskName("JCatalystTask");  // AliJCatalystTask has this name hard coded
	cout << jtTask->GetName() << endl;


	// === Set up JCard ====
	AliJCard *card = new AliJCard(cardName.Data());
	card->PrintOut();
	card->ReadLine( cardSetting.Data() );
	card->ReCompile();
	card->PrintOut();

  // === Create analysis object ===
  
	AliJJtAna *fJtAna;
	fJtAna = new AliJJtAna( kFALSE );

	fJtAna->SetCard( card );
	fJtAna->SetTrigger( jtrigg.Data() );
	fJtAna->SetAssoc( jassoc.Data() );
	if( inclusFileName ) fJtAna->SetInclusiveFile(inclusFileName.Data());

	jtTask->SetJtAnalysis( fJtAna );

	mgr->AddTask((AliAnalysisTask*) jtTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


	// Connect input/output
	mgr->ConnectInput(jtTask, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",jtTask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), jtTask->GetName()));
	mgr->ConnectOutput(jtTask, 1, jHist );

	return jtTask;
}

