//_____________________________________________________________________
AliAnalysisTask *AddTaskJFFlucMapMaster(TString taskName="JFFlucMaster", UInt_t period = 0, double ptmin = 0.5){
	// Load Custom Configuration and parameters
	enum { lhc15o=0, lhc18q=1, lhc18r=2, lhc10h=3};
	cout << "AddTaskJFFlucMapMaster:: period=" << period <<"\t ptmin="<< ptmin << endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	//-------- Correction Maps ----------
	

	if(period == lhc18q) {
		AliJCorrectionMapTask *cmaptask = new AliJCorrectionMapTask("JCorrectionMapTask"); 
    	cmaptask->EnableCentFlattening("alien:///alice/cern.ch/user/j/jparkkil/legotrain/Cent/CentWeights_LHC18q_pass13.root");
    	cmaptask->EnableEffCorrection("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC18q-LHC18l8-0-Lists.root"); // needed but not used!
    	mgr->AddTask((AliAnalysisTask*) cmaptask);
		
    } else if(period == lhc18r) {
    	AliJCorrectionMapTask *cmaptask = new AliJCorrectionMapTask("JCorrectionMapTask"); 
    	cmaptask->EnableCentFlattening("alien:///alice/cern.ch/user/j/jparkkil/legotrain/Cent/CentWeights_LHC18r_pass13.root");
    	cmaptask->EnableEffCorrection("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC18q-LHC18l8-0-Lists.root"); // needed but not used!
    	mgr->AddTask((AliAnalysisTask*) cmaptask);
    } else if(period == lhc10h) {
    	AliJCorrectionMapTask *cmaptask = new AliJCorrectionMapTask("JCorrectionMapTask"); 
    	cmaptask->EnableEffCorrection("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC10h-LHC11a10a_bis-0-Lists.root"); // needed but not used!
    	mgr->AddTask((AliAnalysisTask*) cmaptask);
    }
   
    //-------- JFlucWagons -------
    const int Nsets  = 14; // number of configurations
 	AliJFFlucTask *myTask[Nsets];
	TString configNames[Nsets] = {
		"hybrid", // 0
		"global",
		"nqq",    // 2
		"subA",
		"SPD",    // 4
		"zvtx",
		"pileup",  // 6
		"ntpc100",
		"ntpc90",	// 8
		"ntpc80",
		"dcaxy1",	// 10
		"dcaz2",
		"chi03",	// 12
		"chi35"
	};
   	Int_t hybridCut = 768;
    Int_t globalCut = 96;
    UInt_t selEvt;
	if(period == lhc15o)  {
		selEvt = AliVEvent::kINT7;
	} else if(period == lhc18q || period == lhc18r) {
		selEvt = AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral;
	} else if(period == lhc10h) {
		selEvt = AliVEvent::kMB;
	}
	// --- track related common configurations -----
	for(int i=0;i<Nsets;i++) {
		myTask[i] = new AliJFFlucTask(Form("%s_s_%s",taskName.Data(), configNames[i].Data()));
		if(i!=6) {
			if(period != lhc10h) myTask[i]->AddFlags(AliJFFlucTask::FLUC_EBE_WEIGHTING|AliJFFlucTask::FLUC_CUT_OUTLIERS);
			if(period == lhc18q || period == lhc18r) myTask[i]->AddFlags(AliJFFlucTask::FLUC_CENT_FLATTENING);
		}
		myTask[i]->SelectCollisionCandidates( selEvt );
		myTask[i]->SetCentDetName("V0M");
		myTask[i]->SelectSubevents(AliJFFlucTask::SUBEVENT_A|AliJFFlucTask::SUBEVENT_B);
		myTask[i]->SetTestFilterBit(hybridCut);
		myTask[i]->SetEtaRange(0.4, 0.8);
		myTask[i]->SetPtRange(ptmin, 5.0);
		myTask[i]->SetEffConfig(0,hybridCut);
		myTask[i]->SetZVertexCut(8.);
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
	myTask[iS]->SetZVertexCut(10.);
	// s_pileup
	iS = 6;
	if(period==lhc18q || period==lhc18r) {
		myTask[iS]->AddFlags(AliJFFlucTask::FLUC_EBE_WEIGHTING|AliJFFlucTask::FLUC_CENT_FLATTENING); // CAN|T overwrite Flags!!!
	} else if(period==lhc15o) {
		myTask[iS]->AddFlags(AliJFFlucTask::FLUC_EBE_WEIGHTING); // CAN|T overwrite Flags!!! no weightening for 15o
	}
	//
		//----------- Track related Check -------------
	// s_ntpc 100
	iS = 7;
	myTask[iS]->SetNumTPCClusters(100);
	// s_ntpc 90
	iS = 8;
	myTask[iS]->SetNumTPCClusters(90);
	// s_ntpc 80
	iS = 9;
	myTask[iS]->SetNumTPCClusters(80);
	// s_dcaxy 1 cm
	iS = 10;
	myTask[iS]->SetDCAxyCut(1);
		// s_dcaz 2cm
	iS = 11;
	myTask[iS]->SetDCAzCut(2);
	// s_chi2/ndf in [0.3, 4.0]
	iS = 12;
	myTask[iS]->SetChi2Cuts(0.3, 4.0);
	// s_chi2/ndf in [0.1, 3.5]
	iS = 13;
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

