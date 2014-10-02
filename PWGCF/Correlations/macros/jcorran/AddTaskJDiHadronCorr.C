//_____________________________________________________________________
AliAnalysisTask *AddTaskJDiHadronCorr(TString cardName, TString jtrigg, TString jassoc, TString inclusFileName=""){
    // Load Custom Configuration and parameters
    // override values with parameters

    cout<<"### DEGUG Input is "<< cardName <<"\t"<<jtrigg<<"\t"<<jassoc<<"\t"<<inclusFileName<<"\t"<<"#########"<<endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	//==== JCORRAN Efficiency TASK
	AliJDiHadronCorrTask *jdihadtask = new AliJDiHadronCorrTask("JDiHadronCorrTask","JOD");
	jdihadtask->SetDebugLevel(5);
	jdihadtask->SetFilterTaskName("PWGCFJCORRANTask");


	// === Create AliJCORRAN ====
	AliJCard *card = new AliJCard(cardName.Data());
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
	AliAnalysisDataContainer *jHist = mgr->CreateContainer("JDiHadronCorr",  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:JDiHadronCorr",AliAnalysisManager::GetCommonFileName()));
	mgr->ConnectOutput(jdihadtask, 1, jHist );

	return jdihadtask;
}

