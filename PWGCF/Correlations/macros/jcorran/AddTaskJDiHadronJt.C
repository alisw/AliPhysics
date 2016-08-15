//_____________________________________________________________________
AliAnalysisTask *AddTaskJDiHadronJt(TString taskName, TString cardName, TString jtrigg, TString jassoc, TString cardSetting, TString inclusFileName=""){
	// Load Custom Configuration and parameters
	// override values with parameters

	cout<<"### DEBUG Input is "<< cardName <<"\t"<<jtrigg<<"\t"<<jassoc<<"\t"<<inclusFileName<<"\t"<<"#########"<<endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	//==== Set up di-hadron correlation jT task ====
	AliJDiHadronJtTask *jtTask = new AliJDiHadronJtTask(taskName.Data(),"JOD");
	jtTask->SetDebugLevel(5);
  jtTask->SetFilterTaskName("PWGCFJCORRANTask");  // JCORRAN filter has this name hard coded
	cout << jtTask->GetName() << endl;


	// === Set up JCard ====
	AliJCard *card = new AliJCard(cardName.Data());
	card->PrintOut();
	card->ReadLine( cardSetting.Data() );
	card->ReCompile();
	card->PrintOut();

  // === Create analysis object ===
  
	AliJJtAnalysis *fJtAnalysis;
	fJtAnalysis = new AliJJtAnalysis( kFALSE );

	fJtAnalysis->SetCard( card );
	fJtAnalysis->SetTrigger( jtrigg.Data() );
	fJtAnalysis->SetAssoc( jassoc.Data() );
	if( inclusFileName ) fJtAnalysis->SetInclusiveFile(inclusFileName.Data());

	jtTask->SetJtAnalysis( fJtAnalysis );

	mgr->AddTask((AliAnalysisTask*) jtTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


	// Connect input/output
	mgr->ConnectInput(jtTask, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",jtTask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), jtTask->GetName()));
	mgr->ConnectOutput(jtTask, 1, jHist );

	return jtTask;
}

