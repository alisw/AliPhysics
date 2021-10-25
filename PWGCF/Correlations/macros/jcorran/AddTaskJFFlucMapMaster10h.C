//_____________________________________________________________________
AliAnalysisTask *AddTaskJFFlucMapMaster10h(TString taskName="JFFlucMaster", double ptmin = 0.2){
	// Load Custom Configuration and parameters
	cout << "AddTaskJFFlucMapMaster  ptmin="<< ptmin << endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	//-------- Correction Maps ----------
	/// What about the correction map for the efficiency in LHC10h?

  //-------- JFlucWagons -------
  const int Nsets  = 15; // number of configurations
 	AliJFFlucTask *myTask[Nsets];
	TString configNames[Nsets] = {
		"tpconly", 
		"hybrid", // 1
		"nqq",    // 
		"subA",
		"V0M",    // 4
		"zvtx6",
		"zvtx12",	// 6
		"pileup",	// 7
		"ntpc100",
		"ntpc90",
		"ntpc80",	// 10
		"dcaxy1",
		"dcaz2",	// 12
		"chi03",
		"chi35"	// 14
	};

	Int_t tpconlyCut = 128;
  Int_t hybridCut = 768;
  
  UInt_t selEvt;
	selEvt = AliVEvent::kMB;

	// --- track related common configurations -----
	for(int i=0;i<Nsets;i++) {
		myTask[i] = new AliJFFlucTask(Form("%s_s_%s",taskName.Data(), configNames[i].Data()));

		if (i!=7) {myTask[i]->AddFlags(AliJFFlucTask::FLUC_CUT_OUTLIERS);}	// Apply the outlier cut for all selections except the pileup systematics.

		myTask[i]->SelectCollisionCandidates( selEvt );
		myTask[i]->SetCentDetName("CL1");
		myTask[i]->SelectSubevents(AliJFFlucTask::SUBEVENT_A|AliJFFlucTask::SUBEVENT_B);
		myTask[i]->SetTestFilterBit(tpconlyCut);
		myTask[i]->SetEtaRange(0.4, 0.8); // nothing to do with map
		myTask[i]->SetPtRange(ptmin, 5.0);
		myTask[i]->SetEffConfig(0, tpconlyCut);
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
	// s_VOM
	iS = 4;
	myTask[iS]->SetCentDetName("V0M");
	// s_zvtx 6
	iS = 5;
	myTask[iS]->SetZVertexCut(6);
  // s_zvtx 12
	iS = 6;
	myTask[iS]->SetZVertexCut(12);
	// s_ntpc 100
	iS = 8;
	myTask[iS]->SetNumTPCClusters(100);
	// s_ntpc 90
	iS = 9;
	myTask[iS]->SetNumTPCClusters(90);
	// s_ntpc 80
	iS = 10;
	myTask[iS]->SetNumTPCClusters(80);
	// s_dcaxy 1
	iS = 11;
	myTask[iS]->SetDCAxyCut(1);
	// s_dcaz 2
	iS = 12;
	myTask[iS]->SetDCAzCut(2);
	// s_chi2/ndf in [0.3, 4.0]
	iS = 13;
	myTask[iS]->SetChi2Cuts(0.3, 4.0);
	// s_chi2/ndf in [0.1, 3.5]
	iS = 14;
	myTask[iS]->SetChi2Cuts(0.1, 3.5);

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

