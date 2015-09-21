void AddTask_GammaCaloMerged_pp( 	Int_t 		trainConfig 				= 1,  								// change different set of cuts
									Int_t 		isMC   						= 0, 								// run MC
									Int_t 		enableQAMesonTask 			= 0, 								// enable QA in AliAnalysisTaskGammaCalo
									Int_t 		enableQAClusterTask 		= 0, 								// enable additional QA task
									TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
									TString 	cutnumberAODBranch 			= "000000006008400001001500000",
									TString 	periodname 					= "LHC12f1x", 						// period name
									Bool_t 		doWeighting 				= kFALSE,							// enables weighting
									Int_t 		enableExtQA					= 0,								// enable QA(3), disabled (0)
									Bool_t 		enableTriggerMimicking		= kFALSE,							// enable trigger mimicking
									Bool_t 		enableTriggerOverlapRej		= kFALSE,							// enable trigger overlap rejection
									Float_t		maxFacPtHard				= 3.,								// maximum factor between hardest jet and ptHard generated
									TString		periodNameV0Reader			= "",								// period Name for respective period selected in V0Reader
									Int_t 		selectedMeson				=1 
) {
	
	Int_t isHeavyIon = 0;
	
	// ================== GetAnalysisManager ===============================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error(Form("AddTask_GammaCaloMerged_pp_%i",trainConfig), "No analysis manager found.");
		return ;
	}

	// ================== GetInputEventHandler =============================
	AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

	Bool_t isMCForOtherSettings = 0;
	if (isMC > 0) isMCForOtherSettings = 1;
	//========= Add PID Reponse to ANALYSIS manager ====
	if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
		gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
		AddTaskPIDResponse(isMCForOtherSettings);
	}
	
	Printf("here \n");
	
	//=========  Set Cutnumber for V0Reader ================================
	TString cutnumberPhoton = "00000008400100001500000000";
	TString cutnumberEvent = "00000003";
	Bool_t doEtaShift = kFALSE;
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

	//========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
	if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
		AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
		if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
		fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
		fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
		fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

		if (!mgr) {
			Error("AddTask_V0ReaderV1", "No analysis manager found.");
			return;
		}
		AliConvEventCuts *fEventCuts=NULL;
		if(cutnumberEvent!=""){
			fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
			fEventCuts->SetPreSelectionCutFlag(kTRUE);
			if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
				fEventCuts->DoEtaShift(doEtaShift);
				fV0ReaderV1->SetEventCuts(fEventCuts);
				fEventCuts->SetFillCutHistograms("",kTRUE);
			}
		}
		// Set AnalysisCut Number
		AliConversionPhotonCuts *fCuts=NULL;
		if(cutnumberPhoton!=""){
			fCuts= new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
			fCuts->SetPreSelectionCutFlag(kTRUE);
			fCuts->SetIsHeavyIon(isHeavyIon);
			if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
				fV0ReaderV1->SetConversionCuts(fCuts);
				fCuts->SetFillCutHistograms("",kTRUE);
			}
		}
		if(inputHandler->IsA()==AliAODInputHandler::Class()){
		// AOD mode
			fV0ReaderV1->SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
		}
		fV0ReaderV1->Init();

		AliLog::SetGlobalLogLevel(AliLog::kFatal);

		//connect input V0Reader
		mgr->AddTask(fV0ReaderV1);
		mgr->ConnectInput(fV0ReaderV1,0,cinput);
	}

	//================================================
	//========= Add task to the ANALYSIS manager =====
	//================================================
	AliAnalysisTaskGammaCaloMerged *task=NULL;
	task= new AliAnalysisTaskGammaCaloMerged(Form("GammaCaloMerged_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 2;
	if (trainConfig == 111 	|| trainConfig == 112 	|| trainConfig == 113	|| trainConfig == 114  )
		numberOfCuts = 3;
	if (trainConfig == 1 	|| trainConfig == 2 	||
		trainConfig == 11	|| trainConfig == 12 	|| trainConfig == 13	|| trainConfig == 14 )
		numberOfCuts = 4;

	TString *eventCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *clusterMergedCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// cluster cuts
	// 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
	// 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
	
	// ************************************* EMCAL cuts ****************************************************
	// LHC11a
	if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111100050032200000"; clusterMergedCutArray[0] = "1111100050022110001"; mesonCutArray[0] = "0163101100000000"; // MB NLM 1
		eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111100050032200000"; clusterMergedCutArray[1] = "1111100050022110001"; mesonCutArray[1] = "0163101100000000"; // EMC1 NLM 2
		eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111100050032200000"; clusterMergedCutArray[2] = "1111100050022110002"; mesonCutArray[2] = "0163102200000000"; // MB NLM 2
		eventCutArray[ 3] = "00051113"; clusterCutArray[3] = "1111100050032200000"; clusterMergedCutArray[3] = "1111100050022110002"; mesonCutArray[3] = "0163102200000000"; // EMC1 NLM 2
	} else if (trainConfig == 2){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050032200000"; clusterMergedCutArray[0] = "1111121050022110001"; mesonCutArray[0] = "0163101100000000"; // MB NLM 1
		eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121050032200000"; clusterMergedCutArray[1] = "1111121050022110001"; mesonCutArray[1] = "0163101100000000"; // EMC1 NLM 2
		eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111121050032200000"; clusterMergedCutArray[2] = "1111121050022110002"; mesonCutArray[2] = "0163102200000000"; // MB NLM 2
		eventCutArray[ 3] = "00051113"; clusterCutArray[3] = "1111121050032200000"; clusterMergedCutArray[3] = "1111121050022110002"; mesonCutArray[3] = "0163102200000000"; // EMC1 NLM 2

	// LHC13g
	} else if (trainConfig == 10){  // EMCAL clusters, EMCEGA triggers
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111100050032200000"; clusterMergedCutArray[0] = "1111100050022110001"; mesonCutArray[0] = "0163100000000000"; // INT7
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111100050032200000"; clusterMergedCutArray[1] = "1111100050022110002"; mesonCutArray[1] = "0163102200000000"; // INT7
	} else if (trainConfig == 11){  // EMCAL clusters, EMCEGA trigger, 1 NLM
		eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111100050032200000"; clusterMergedCutArray[0] = "1111100050022110001"; mesonCutArray[0] = "0163100000000000"; // INT7
		eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111100050032200000"; clusterMergedCutArray[1] = "1111100050022110001"; mesonCutArray[1] = "0163100000000000"; // INT7
		eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111100050032200000"; clusterMergedCutArray[2] = "1111100050022110001"; mesonCutArray[2] = "0163100000000000"; // INT7
		eventCutArray[ 3] = "00052113"; clusterCutArray[3] = "1111100050032200000"; clusterMergedCutArray[3] = "1111100050022110001"; mesonCutArray[3] = "0163100000000000"; // INT7
	} else if (trainConfig == 12){  // EMCAL clusters, EMCEGA trigger, 2 NLM
		eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111100050032200000"; clusterMergedCutArray[0] = "1111100050022110002"; mesonCutArray[0] = "0163102200000000"; // INT7
		eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111100050032200000"; clusterMergedCutArray[1] = "1111100050022110002"; mesonCutArray[1] = "0163102200000000"; // INT7
		eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111100050032200000"; clusterMergedCutArray[2] = "1111100050022110002"; mesonCutArray[2] = "0163102200000000"; // INT7
		eventCutArray[ 3] = "00052113"; clusterCutArray[3] = "1111100050032200000"; clusterMergedCutArray[3] = "1111100050022110002"; mesonCutArray[3] = "0163102200000000"; // INT7
	} else if (trainConfig == 13){  // EMCAL clusters, EMCEGA trigger, 1 NLM
		eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111121050032200000"; clusterMergedCutArray[0] = "1111121050022110001"; mesonCutArray[0] = "0163100000000000"; // INT7
		eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111121050032200000"; clusterMergedCutArray[1] = "1111121050022110001"; mesonCutArray[1] = "0163100000000000"; // INT7
		eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111121050032200000"; clusterMergedCutArray[2] = "1111121050022110001"; mesonCutArray[2] = "0163100000000000"; // INT7
		eventCutArray[ 3] = "00052113"; clusterCutArray[3] = "1111121050032200000"; clusterMergedCutArray[3] = "1111121050022110001"; mesonCutArray[3] = "0163100000000000"; // INT7
	} else if (trainConfig == 14){  // EMCAL clusters, EMCEGA trigger, 2 NLM
		eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111121050032200000"; clusterMergedCutArray[0] = "1111121050022110002"; mesonCutArray[0] = "0163102200000000"; // INT7
		eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111121050032200000"; clusterMergedCutArray[1] = "1111121050022110002"; mesonCutArray[1] = "0163102200000000"; // INT7
		eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111121050032200000"; clusterMergedCutArray[2] = "1111121050022110002"; mesonCutArray[2] = "0163102200000000"; // INT7
		eventCutArray[ 3] = "00052113"; clusterCutArray[3] = "1111121050032200000"; clusterMergedCutArray[3] = "1111121050022110002"; mesonCutArray[3] = "0163102200000000"; // INT7

	// LHC12
    // default with three cuts
	} else if (trainConfig == 111){  // EMCAL clusters, different triggers no NonLinearity
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111100060032200000"; clusterMergedCutArray[0] = "1111100060022110001"; mesonCutArray[0] = "0163100000000000"; // INT7
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111100060032200000"; clusterMergedCutArray[1] = "1111100060022110001"; mesonCutArray[1] = "0163100000000000"; // EMC7
		eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111100060032200000"; clusterMergedCutArray[2] = "1111100060022110001"; mesonCutArray[2] = "0163100000000000"; // EMCEG1,
	} else if (trainConfig == 112){  // EMCAL clusters, different triggers no NonLinearity
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111100060032200000"; clusterMergedCutArray[0] = "1111100060022110002"; mesonCutArray[0] = "0163102200000000"; // INT7
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111100060032200000"; clusterMergedCutArray[1] = "1111100060022110002"; mesonCutArray[1] = "0163102200000000"; // EMC7
		eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111100060032200000"; clusterMergedCutArray[2] = "1111100060022110002"; mesonCutArray[2] = "0163102200000000"; // EMCEG1,
	} else if (trainConfig == 113){  // EMCAL clusters, different triggers with NonLinearity
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060032200000"; clusterMergedCutArray[0] = "1111111060022110001"; mesonCutArray[0] = "0163100000000000"; // INT7
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111111060032200000"; clusterMergedCutArray[1] = "1111111060022110001"; mesonCutArray[1] = "0163100000000000"; // EMC7
		eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111111060032200000"; clusterMergedCutArray[2] = "1111111060022110001"; mesonCutArray[2] = "0163100000000000"; // EMCEG1,
	} else if (trainConfig == 114){  // EMCAL clusters, different triggers with NonLinearity
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060032200000"; clusterMergedCutArray[0] = "1111111060022110002"; mesonCutArray[0] = "0163102200000000"; // INT7
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111111060032200000"; clusterMergedCutArray[1] = "1111111060022110002"; mesonCutArray[1] = "0163102200000000"; // EMC7
		eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111111060032200000"; clusterMergedCutArray[2] = "1111111060022110002"; mesonCutArray[2] = "0163102200000000"; // EMCEG1,
    
    // all the cut variations
	} else {
		Error(Form("GammaCaloMerged_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
	TList *ClusterCutList = new TList();
	TList *ClusterMergedCutList = new TList();
	TList *MesonCutList = new TList();

	TList *HeaderList = new TList();
	if (periodname.Contains("LHC12i3")){	
		TObjString *Header2 = new TObjString("BOX");
		HeaderList->Add(Header2);
	} else if (periodname.CompareTo("LHC14e2b")==0){
		TObjString *Header2 = new TObjString("pi0_1");
		HeaderList->Add(Header2);
		TObjString *Header3 = new TObjString("eta_2");
		HeaderList->Add(Header3);
	}	
	
	TString energy = "";
	TString mcName = "";
	TString mcNameAdd = "";
	if (periodname.Contains("WOSDD")){
		mcNameAdd = "_WOSDD";
	} else if (periodname.Contains("WSDD")){
		mcNameAdd = "_WSDD";
	} 	
	if (periodname.Contains("LHC12i3")){
		energy = "2760GeV";
		mcName = "Pythia8_LHC12i3";
	} else if (periodname.Contains("LHC12f1a")){	
		energy = "2760GeV";
		mcName = "Pythia8_LHC12f1a";	
	} else if (periodname.Contains("LHC12f1b")){	
		energy = "2760GeV";
		mcName = "Phojet_LHC12f1b";			
	} else if (periodname.Contains("LHC14e2a")){	
		energy = "8TeV";
		mcName = "Pythia8_LHC14e2a";			
	} else if (periodname.Contains("LHC14e2b")){	
		energy = "8TeV";
		mcName = "Pythia8_LHC14e2b";				
	} else if (periodname.Contains("LHC14e2c")){		
		energy = "8TeV";
		mcName = "Phojet_LHC14e2c";					
	}	
	
	EventCutList->SetOwner(kTRUE);
	AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
	ClusterCutList->SetOwner(kTRUE);
	AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
	ClusterMergedCutList->SetOwner(kTRUE);
	AliCaloPhotonCuts **analysisClusterMergedCuts = new AliCaloPhotonCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

	for(Int_t i = 0; i<numberOfCuts; i++){
		analysisEventCuts[i] = new AliConvEventCuts();   
		
		// definition of weighting input
		TString fitNamePi0 = Form("Pi0_Fit_Data_%s",energy.Data());
		TString fitNameEta = Form("Eta_Fit_Data_%s",energy.Data());

		TString mcInputNamePi0 = "";
		TString mcInputNameEta = "";
		mcInputNamePi0 = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
		mcInputNameEta = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
		
		if (doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

		analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
		analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
		analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
		analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
		EventCutList->Add(analysisEventCuts[i]);
		analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
		
		analysisClusterCuts[i] = new AliCaloPhotonCuts();
		analysisClusterCuts[i]->InitializeCutsFromCutString(clusterCutArray[i].Data());
		ClusterCutList->Add(analysisClusterCuts[i]);
		analysisClusterCuts[i]->SetExtendedQA(enableExtQA);
		analysisClusterCuts[i]->SetFillCutHistograms("");

		analysisClusterMergedCuts[i] = new AliCaloPhotonCuts();
		analysisClusterMergedCuts[i]->SetIsMergedClusterCut(kTRUE);
		analysisClusterMergedCuts[i]->InitializeCutsFromCutString(clusterMergedCutArray[i].Data());
		ClusterMergedCutList->Add(analysisClusterMergedCuts[i]);
		analysisClusterMergedCuts[i]->SetExtendedQA(enableExtQA);
		analysisClusterMergedCuts[i]->SetFillCutHistograms("");

		
		analysisMesonCuts[i] = new AliConversionMesonCuts();
		analysisMesonCuts[i]->SetEnableOpeningAngleCut(kFALSE);
		analysisMesonCuts[i]->SetIsMergedClusterCut(kTRUE);
		analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
		MesonCutList->Add(analysisMesonCuts[i]);
		analysisMesonCuts[i]->SetFillCutHistograms("");
		analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
	}
	task->SetSelectedMesonID(selectedMeson);
	task->SetEventCutList(numberOfCuts,EventCutList);
	task->SetCaloCutList(numberOfCuts,ClusterCutList);
	task->SetCaloMergedCutList(numberOfCuts,ClusterMergedCutList);
	task->SetMesonCutList(numberOfCuts,MesonCutList);
	task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
	task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
	if(enableExtQA == 3){ task->SetPlotHistsExtQA(kTRUE);}
	
	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer(Form("GammaCaloMerged_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaCaloMerged_%i.root",trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
