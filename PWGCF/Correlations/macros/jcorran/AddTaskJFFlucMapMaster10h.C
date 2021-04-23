//_____________________________________________________________________
AliAnalysisTask *AddTaskJFFlucMapMaster10h(TString taskName="JFFlucMaster", double ptmin = 0.2){
	// Load Custom Configuration and parameters
	cout << "AddTaskJFFlucMapMaster  ptmin="<< ptmin << endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	//-------- Correction Maps ----------
	
    //-------- JFlucWagons -------
    const int Nsets  = 7; // number of configurations
 	AliJFFlucTask *myTask[Nsets];
	TString configNames[Nsets] = {
		"tpconly", 
		"hybrid", // 1
		"nqq",    // 
		"subA",
		"SPD",    // 4
		"zvtx",
		"pileup"  // 6
	};
	Int_t tpconlyCut = 128;
   	Int_t hybridCut = 768;

    UInt_t selEvt;
	selEvt = AliVEvent::kMB;
	// --- track related common configurations -----
	for(int i=0;i<Nsets;i++) {
		myTask[i] = new AliJFFlucTask(Form("%s_s_%s",taskName.Data(), configNames[i].Data()));
		if(i!=6) {
			myTask[i]->AddFlags(AliJFFlucTask::FLUC_CUT_OUTLIERS|AliJFFlucTask::FLUC_EBE_WEIGHTING);
		}
		myTask[i]->SelectCollisionCandidates( selEvt );
		myTask[i]->SetCentDetName("V0M");
		myTask[i]->SelectSubevents(AliJFFlucTask::SUBEVENT_A|AliJFFlucTask::SUBEVENT_B);
		myTask[i]->SetTestFilterBit(tpconlyCut);
		myTask[i]->SetEtaRange(0.4, 0.8);
		myTask[i]->SetPtRange(ptmin, 5.0);
		myTask[i]->SetEffConfig(0,tpconlyCut);
	}
	// s_hybrid
	int iS = 1;
	myTask[iS]->SetTestFilterBit(hybridCut);
	myTask[iS]->SetEffConfig(0,hybridCut);
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
	// s_pileup
	iS = 6;
	myTask[iS]->AddFlags(AliJFFlucTask::FLUC_EBE_WEIGHTING); 	
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

