//_____________________________________________________________________
AliAnalysisTask *AddTaskJDiHadronCorr(TString taskName, TString cardName, TString jtrigg, TString jassoc, TString cardSetting, TString inclusFileName=""){
	// Load Custom Configuration and parameters
	// override values with parameters

	cout<<"### DEGUG Input is "<< cardName <<"\t"<<jtrigg<<"\t"<<jassoc<<"\t"<<inclusFileName<<"\t"<<"#########"<<endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	//==== JCORRAN Efficiency TASK
	AliJDiHadronCorrTask *jdihadtask = new AliJDiHadronCorrTask(taskName.Data(),"JOD");
	jdihadtask->SetDebugLevel(5);
	jdihadtask->SetFilterTaskName("PWGCFJCORRANTask");
	cout << jdihadtask->GetName() << endl;


	// === Create AliJCORRAN ====
	AliJCard *card = new AliJCard(cardName.Data());
	card->PrintOut();
	card->ReadLine( cardSetting.Data() );
	card->ReCompile();
	card->PrintOut();

	AliJCORRAN *fJCORRAN;
	fJCORRAN = new AliJCORRAN( kFALSE );

	fJCORRAN->SetCard( card );
	fJCORRAN->SetTrigger( jtrigg.Data() );
	fJCORRAN->SetAssoc( jassoc.Data() );
	if( inclusFileName ) fJCORRAN->SetInclusiveFile(inclusFileName.Data());

	jdihadtask->SetJCORRAN( fJCORRAN );

	mgr->AddTask((AliAnalysisTask*) jdihadtask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


	// Connect input/output
	mgr->ConnectInput(jdihadtask, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",jdihadtask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), jdihadtask->GetName()));
	mgr->ConnectOutput(jdihadtask, 1, jHist );

	return jdihadtask;
}

