//_____________________________________________________________________
AliAnalysisTask *AddTaskJFFlucMapMCMaster(TString taskName="JFFlucMaster", double ptmin = 0.5){
	// Load Custom Configuration and parameters
	// period < 0 MC
	cout << "AddTaskJFFlucMapMaster:: period=" << period <<"\t ptmin="<< ptmin << endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	//-------- Correction Maps ----------
    //-------- JFlucWagons -------
    const int Nsets  = 7; // number of configurations
 	AliJFFlucTask *myTask[Nsets];
	TString configNames[Nsets] = {
		"hybrid", // 0
		"global",
		"nqq",    // 2
		"subA",
		"SPD",    // 4
		"zvtx",
		"pqq"  // 6
	};
   	Int_t hybridCut = 768;
    Int_t globalCut = 96;
    UInt_t selEvt;
	
	selEvt = AliVEvent::kAny;
	
	// --- track related common configurations -----
	for(int i=0;i<Nsets;i++) {
		myTask[i] = new AliJFFlucTask(Form("%s_s_%s",taskName.Data(), configNames[i].Data()));
	    myTask[i]->AddFlags(AliJFFlucTask::FLUC_MC|AliJFFlucTask::FLUC_EXCLUDEWDECAY);
		myTask[i]->SelectCollisionCandidates( selEvt );
		myTask[i]->SetCentDetName("V0M");
		myTask[i]->SelectSubevents(AliJFFlucTask::SUBEVENT_A|AliJFFlucTask::SUBEVENT_B);
		myTask[i]->SetTestFilterBit(hybridCut);
		myTask[i]->SetEtaRange(0.4, 0.8);
		myTask[i]->SetPtRange(ptmin, 5.0);
		myTask[i]->SetEffConfig(0,hybridCut);
	}
	// s_global
	int iS = 1;
	myTask[iS]->SetTestFilterBit(globalCut);
	myTask[iS]->SetEffConfig(0,globalCut);
	// s_nqq
	iS = 2;
	myTask[iS]->SetParticleCharge(-1);
    // s_subA
    iS = 3;
	myTask[iS]->SelectSubevents(AliJFFlucTask::SUBEVENT_A); // subA
	//
	//----------- Event related Check -------------
	// s_SPD
	iS = 4;
	myTask[iS]->SetCentDetName("CL1");
	// s_zvtx
	iS = 5;
	myTask[iS]->SetZVertexCut(8);
	// s_pqq
	myTask[iS]->SetParticleCharge(1);
	
	// Must add the tasks
	for(int i=0;i<Nsets;i++) mgr->AddTask((AliAnalysisTask*) myTask[i]);
	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

	// Connect input/output
	for(int i=0;i<Nsets;i++) {
		mgr->ConnectInput(myTask[i], 0, cinput);
		AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",myTask[i]->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), myTask[i]->GetName()));
		mgr->ConnectOutput(myTask[i], 1, jHist );
	}

	return myTask[0];
}

