//_____________________________________________________________________
AliAnalysisTask *AddTaskJFFlucMasterPubMaps(TString taskName="JFFlucMasterPubMaps", UInt_t period = 0, double ptmin = 0.5, Bool_t removebadarea = kFALSE, Bool_t mapsmooth = kFALSE){
	// Load Custom Configuration and parameters
	enum { lhc15o=0, lhc18q=1, lhc18r=2 };
	const TString speriod[3]= {"15o","18q","18r"}; //needed string to load correct map config based on string
	cout << "AddTaskJFFlucMaster:: period=" << period <<"\t ptmin="<< ptmin << endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	//-------- Loading Correction Maps ----------
	
	//-------- JFlucWagons -------
	const int Nsets  = 1; // number of configurations
 	AliJFFlucTask *myTask[Nsets];
	const TString configNames[Nsets] = {
		"hybrid", // 0
		//"global",
		//"nqq",    // 2
		//"subA",
		//"SPD",    // 4
		//"zvtx",
		//"pileup"  // 6
	};
	//loading correction map
	TString s_smooth = "";
	if(mapsmooth) s_smooth = "Smooth_"; //enable when needed again
	//const TString MAPdirname="alien:///alice/cern.ch/user/a/aonnerst/legotrain/NUAError/";
	const TString MAPdirname="alien:///alice/cern.ch/user/j/jparkkil/legotrain/NUA/";
	AliJCorrectionMapTask *cmaptask = new AliJCorrectionMapTask("JCorrectionMapTask");
	if(period == lhc18q || period == lhc18r) {
		cmaptask->EnableCentFlattening(Form("alien:///alice/cern.ch/user/j/jparkkil/legotrain/Cent/CentWeights_LHC%s_pass13.root",speriod[period].Data()));//centrality flattening
		cmaptask->EnableEffCorrection(Form("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC%s-LHC18l8-0-Lists.root",speriod[period].Data()));//efficiency cirrection
	}
	if(period == lhc15o) {
		cmaptask->EnableEffCorrection(Form("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC%s-LHC16g-0-Lists.root",speriod[period].Data()));//efficiency cirrection
	}
	const TString mapFilesLHC15o[Nsets] = { //correction maps used for 2020 publication
		"PhiWeights_LHC15o_hybrid_3d-rebin2-ROOT5.root",
		//"PhiWeights_LHC15o_global_3d-rebin-ROOT5.root",
		//"PhiWeights_LHC15o_s_charge_nqq_3d-rebin2-ROOT5.root",
		//"PhiWeights_LHC15o_hybrid_3d-rebin2-ROOT5.root",
		//"PhiWeights_LHC15o_s_SPD_3d-rebin2-ROOT5.root",
		//"PhiWeights_LHC15o_s_zvtx_3d-rebin2-ROOT5.root",
		//"PhiWeights_LHC15o_s_pileup_3d-rebin2-ROOT5.root"
	};
	//cout<<"Using LHC15o correction maps.\n";
	//Load phi correction maps for a limited set of configurations. The rest of the configs will use the map for the default configuration.
	const UInt_t mapIndices[] = {0};//,1,2,4,5};
	for(UInt_t i = 0; i < sizeof(mapIndices)/sizeof(mapIndices[0]); ++i){
		UInt_t ci = mapIndices[i];
		//TString MAPfilename = Form("%sPhiWeights_LHC%s_Error_%spt%02d_s_%s.root",MAPdirname.Data(), speriod[period].Data(),s_smooth.Data(), Int_t(ptmin*10), configNames[ci].Data()); //azimuthal correction
		TString MAPfilename = MAPdirname+mapFilesLHC15o[ci];
		cout << MAPfilename.Data() << endl;
		cmaptask->EnablePhiCorrection(i,MAPfilename);
	}
	mgr->AddTask((AliAnalysisTask*) cmaptask);
    
   
   	Int_t hybridCut = 768;
    Int_t globalCut = 96;
    UInt_t selEvt;
	if(period == lhc15o)  selEvt = AliVEvent::kINT7;
	else if(period == lhc18q || period == lhc18r)  selEvt = AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral;
	// --- track related common configurations -----
	for(int i=0;i<Nsets;i++) {
		myTask[i] = new AliJFFlucTask(Form("%s_s_%s",taskName.Data(), configNames[i].Data()));
		if(i!=6) {
			myTask[i]->AddFlags(AliJFFlucTask::FLUC_EBE_WEIGHTING|AliJFFlucTask::FLUC_CUT_OUTLIERS|AliJFFlucTask::FLUC_PHI_CORRECTION);// AliJFFlucTask::FLUC_PHI_CORRECTION must be enabled
			if(period == lhc18q || period == lhc18r) myTask[i]->AddFlags(AliJFFlucTask::FLUC_CENT_FLATTENING);
		}
		myTask[i]->SelectCollisionCandidates( selEvt );
		myTask[i]->SetCentDetName("V0M");
		myTask[i]->SelectSubevents(AliJFFlucTask::SUBEVENT_A|AliJFFlucTask::SUBEVENT_B);
		myTask[i]->SetTestFilterBit(hybridCut);
		myTask[i]->SetEtaRange(0.4, 0.8);
		myTask[i]->SetPtRange(ptmin, 5.0);
		myTask[i]->SetEffConfig(1,hybridCut);
		//myTask[i]->SetPhiCorrectionIndex(i);//cmaptask->EnablePhiCorrection(i,MAPfilenames[i]);
		//Find if this configuration uses its own phi correction map, otherwise, use default.
		UInt_t mapIndex = 0;
		for(UInt_t j = 0; j < sizeof(mapIndices)/sizeof(mapIndices[0]); ++j)
			if(mapIndices[j] == i){
				mapIndex = j;
				break;
			}
		myTask[i]->SetPhiCorrectionIndex(mapIndex);//cmaptask->EnablePhiCorrection(i,MAPfilenames[i]);
		myTask[i]->SetRemoveBadArea(removebadarea);
		myTask[i]->SetZVertexCut(8.);
	}
	// s_global
	/*int iS = 1;
	myTask[iS]->SetTestFilterBit(globalCut);
	myTask[iS]->SetEffConfig(1,globalCut);
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
	myTask[iS]->AddFlags(AliJFFlucTask::FLUC_EBE_WEIGHTING|AliJFFlucTask::FLUC_PHI_CORRECTION); // CAN|T overwrite Flags!!! no weightening for 15o
	if(period == lhc18q || period == lhc18r) myTask[iS]->AddFlags(AliJFFlucTask::FLUC_CENT_FLATTENING);*/
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

