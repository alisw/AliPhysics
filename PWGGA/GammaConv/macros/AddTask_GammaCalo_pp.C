void AddTask_GammaCalo_pp(  Int_t 		trainConfig 				= 1,  								// change different set of cuts
							Int_t 		isMC   						= 0, 							// run MC
							Int_t 		enableQAMesonTask 			= 0, 								// enable QA in AliAnalysisTaskGammaCalo
							Int_t 		enableQAClusterTask 		= 0, 								// enable additional QA task
							TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
                            TString 	cutnumberAODBranch 			= "000000006008400001001500000",
							TString 	periodname 					= "LHC12f1x", 						// period name
							Bool_t 		doWeighting 				= kFALSE,							// enables weighting
							Bool_t 		isUsingTHnSparse 			= kTRUE,							// enable or disable usage of THnSparses for background estimation
							Int_t 		enableExtQA					= 0,								// enable QA(3), disabled (0)
							Bool_t 		enableTriggerMimicking		= kFALSE,							// enable trigger mimicking
							Bool_t 		enableTriggerOverlapRej		= kFALSE,							// enable trigger overlap rejection
							Float_t		maxFacPtHard				= 3.									// maximum factor between hardest jet and ptHard generated
) {

	// ================= Load Librariers =================================
	gSystem->Load("libCore");  
	gSystem->Load("libTree");
	gSystem->Load("libGeom");
	gSystem->Load("libVMC");
	gSystem->Load("libPhysics");
	gSystem->Load("libMinuit");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");  
	gSystem->Load("libCDB");
	gSystem->Load("libSTEER");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libTender");
	gSystem->Load("libTenderSupplies");
	gSystem->Load("libPWGflowBase");
	gSystem->Load("libPWGflowTasks");
	gSystem->Load("libPWGGAGammaConv");
	
	Int_t isHeavyIon = 0;
	
	// ================== GetAnalysisManager ===============================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error(Form("AddTask_GammaCalo_pp_%i",trainConfig), "No analysis manager found.");
		return ;
	}

	// ================== GetInputEventHandler =============================
	AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
	//**********************************************************************************************
	//**********************************************************************************************
		AliInputEventHandler *inputEventHandler=mgr->GetInputEventHandler();
		inputEventHandler->SetInactiveBranches("AliESDFMD"); // Disable FMD branch, see ALIROOT-6222
	//**********************************************************************************************
	//**********************************************************************************************
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
	AliAnalysisTaskGammaCalo *task=NULL;
	task= new AliAnalysisTaskGammaCalo(Form("GammaCalo_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 2;
	if (trainConfig == 32 	|| trainConfig == 101 || trainConfig == 201 || trainConfig == 63)
		numberOfCuts = 1;
	if (trainConfig == 111 	|| trainConfig == 114 	|| trainConfig == 117 	|| trainConfig == 120 	|| trainConfig == 121 || 
		trainConfig == 72 	|| trainConfig == 122 	|| trainConfig == 123 	|| trainConfig == 124 	|| trainConfig == 75  || 
		trainConfig == 78 	|| trainConfig == 81	 )
		numberOfCuts = 3;
	if (trainConfig == 31 	|| trainConfig == 4 	|| trainConfig == 6 	|| trainConfig == 53  || trainConfig == 56 	||
		trainConfig == 59 	|| trainConfig == 62 	|| trainConfig == 64 	|| trainConfig == 102 || trainConfig == 110 )
		numberOfCuts = 4;
	if (trainConfig == 2 	|| trainConfig == 3		|| trainConfig == 103	|| trainConfig == 104)
		numberOfCuts = 5;
	if (trainConfig == 5) 
		numberOfCuts = 6;


  
	TString *eventCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// cluster cuts
	// 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
	// 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
	
	// ************************************* EMCAL cuts ****************************************************
	// LHC11a
	if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050032230000"; mesonCutArray[0] = "0163103100000000"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121050032230000"; mesonCutArray[1] = "0163103100000000"; // 400 MeV cluster min energy
	} else if (trainConfig == 2){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050012230000"; mesonCutArray[0] = "0163103100000000"; //0.2 GeV/c
		eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050022230000"; mesonCutArray[1] = "0163103100000000"; //0.3 GeV/c
		eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111121050032230000"; mesonCutArray[2] = "0163103100000000"; //0.4 GeV/c default
		eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1111121050042230000"; mesonCutArray[3] = "0163103100000000"; //0.5 GeV/c
		eventCutArray[ 4] = "00003113"; clusterCutArray[4] = "1111121050052230000"; mesonCutArray[4] = "0163103100000000"; //0.6 GeV/c
	} else if (trainConfig == 3){ //EMCAL minNCells variation
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050031230000"; mesonCutArray[0] = "0163103100000000"; //n cells >= 1
		eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050033230000"; mesonCutArray[1] = "0163103100000000"; //n cells >= 3
		eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111121050032000000"; mesonCutArray[2] = "0163103100000000"; //no M02 cut
		eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1113121050032230000"; mesonCutArray[3] = "0163103100000000"; //only modules with TRD infront
		eventCutArray[ 4] = "00003113"; clusterCutArray[4] = "1111221050032230000"; mesonCutArray[4] = "0163103100000000"; //no modules with TRD infront
	}else if (trainConfig == 4){ // EMCAL clusters 2.76 TeV NonLinearity
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111101050032230000"; mesonCutArray[0] = "0163103100000000"; // NonLinearity kSDMv5
		eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050032230000"; mesonCutArray[1] = "0163103100000000"; // NonLinearity LHC11a ConvCalo
		eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111122050032230000"; mesonCutArray[2] = "0163103100000000"; // NonLinearity LHC11a Calo
		eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1111100050032230000"; mesonCutArray[3] = "0163103100000000"; // NonLinearity none
	} else if (trainConfig == 6){ // EMCAL clusters 2.76 TeV NonLinearity
		eventCutArray[ 0] = "00051113"; clusterCutArray[0] = "1111101050032230000"; mesonCutArray[0] = "0163103100000000"; // NonLinearity kSDMv5
		eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121050032230000"; mesonCutArray[1] = "0163103100000000"; // NonLinearity LHC11a ConvCalo
		eventCutArray[ 2] = "00051113"; clusterCutArray[2] = "1111122050032230000"; mesonCutArray[2] = "0163103100000000"; // NonLinearity LHC11a Calo
		eventCutArray[ 3] = "00051113"; clusterCutArray[3] = "1111100050032230000"; mesonCutArray[3] = "0163103100000000"; // NonLinearity none
	// LHC13g	
	} else if (trainConfig == 5){  // EMCAL clusters, EMCEGA triggers
		eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111100050032220000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1,
		eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111100050032220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG2,
		eventCutArray[ 2] = "00093113"; clusterCutArray[2] = "1111100050032220000"; mesonCutArray[2] = "0163103100000000"; // EMCEJ1,
		eventCutArray[ 3] = "00095113"; clusterCutArray[3] = "1111100050032220000"; mesonCutArray[3] = "0163103100000000"; // EMCEJ2,
		eventCutArray[ 4] = "00000113"; clusterCutArray[4] = "1111100050032220000"; mesonCutArray[4] = "0163103100000000"; // INT7
		eventCutArray[ 5] = "00052113"; clusterCutArray[5] = "1111100050032220000"; mesonCutArray[5] = "0163103100000000"; // EMC7
	} else if (trainConfig == 12){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0) without and with added signals
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111100050032220000"; mesonCutArray[0] = "0163103100000000"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "00003123"; clusterCutArray[1] = "1111100050032220000"; mesonCutArray[1] = "0163103100000000"; // 400 MeV cluster min energy
	} else if (trainConfig == 13){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
		eventCutArray[ 0] = "00003123"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[0] = "0163103100000000"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "00051123"; clusterCutArray[1] = "1111100050032230000"; mesonCutArray[1] = "0163103100000000"; // 400 MeV cluster min energy
	// ************************************* PHOS cuts ****************************************************
	} else if (trainConfig == 31) { //PHOS clusters
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "2444400040033200000"; mesonCutArray[0] = "0163103100000000"; //pp LHC11a with SDD, PHOS
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "2444400040033200000"; mesonCutArray[1] = "0163103100000000"; //pp LHC13g default MB
		eventCutArray[ 2] = "00061113"; clusterCutArray[2] = "2444400040033200000"; mesonCutArray[2] = "0163103100000000"; //pp LHC11a PHI1
		eventCutArray[ 3] = "00062113"; clusterCutArray[3] = "2444400040033200000"; mesonCutArray[3] = "0163103100000000"; //pp LHC11a PHI7
	} else if (trainConfig == 32){ // Validation PHOS
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "2444400040033200000"; mesonCutArray[0] = "0163003100900000";
	} else if (trainConfig == 33){ // PHOS clusters, without and with added signals
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "2444400040033200000"; mesonCutArray[0] = "0163003100900000";
		eventCutArray[ 1] = "00003123"; clusterCutArray[1] = "2444400040033200000"; mesonCutArray[1] = "0163003100900000";

    // LHC13g cut studies
	} else if (trainConfig == 51){  // EMCAL clusters, EMCEG1 trigger
		eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1, 400 MeV min energy, NCells >=2, M02 default cut, 50ns timing
		eventCutArray[ 1] = "00083113"; clusterCutArray[1] = "1111121060052220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG1, 600 MeV min energy
	} else if (trainConfig == 52){  // EMCAL clusters, EMCEG1 trigger
		eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111121060031220000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1,                     NCells >=1
		eventCutArray[ 1] = "00083113"; clusterCutArray[1] = "1111121060033220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG1,                     NCells >=3
	} else if (trainConfig == 53){  // EMCAL clusters, EMCEG1 trigger
		eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111121060032000000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1,                                 no M02 cut
		eventCutArray[ 1] = "00083113"; clusterCutArray[1] = "1111101060032220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG1, 						 	standard kSDMv5
		eventCutArray[ 2] = "00083113"; clusterCutArray[2] = "1111122060032220000"; mesonCutArray[2] = "0163103100000000"; // EMCEG1,							NonLinearity LHC11a Calo
		eventCutArray[ 3] = "00083113"; clusterCutArray[3] = "1111100060032220000"; mesonCutArray[3] = "0163103100000000"; // EMCEG1,							NonLinearity none
		
	} else if (trainConfig == 54){  // EMCAL clusters, EMCEG2 trigger
		eventCutArray[ 0] = "00085113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000000"; // EMCEG2, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111121060052220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG2, 600 MeV min energy
	} else if (trainConfig == 55){  // EMCAL clusters, EMCEG2 trigger
		eventCutArray[ 0] = "00085113"; clusterCutArray[0] = "1111121060031220000"; mesonCutArray[0] = "0163103100000000"; // EMCEG2,                     NCells >=1
		eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111121060033220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG2,                     NCells >=3
	} else if (trainConfig == 56){  // EMCAL clusters, EMCEG2 trigger
		eventCutArray[ 0] = "00085113"; clusterCutArray[0] = "1111121060032000000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1,                                 no M02 cut
		eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111101060032220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG1, 						 	standard kSDMv5
		eventCutArray[ 2] = "00085113"; clusterCutArray[2] = "1111122060032220000"; mesonCutArray[2] = "0163103100000000"; // EMCEG1,							NonLinearity LHC11a Calo
		eventCutArray[ 3] = "00085113"; clusterCutArray[3] = "1111100060032220000"; mesonCutArray[3] = "0163103100000000"; // EMCEG1,							NonLinearity none
	
	} else if (trainConfig == 57){  // EMCAL clusters, INT7 trigger
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000000"; // INT7, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111121060052220000"; mesonCutArray[1] = "0163103100000000"; // INT7, 600 MeV min energy
	} else if (trainConfig == 58){  // EMCAL clusters, INT7 trigger
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111121060031220000"; mesonCutArray[0] = "0163103100000000"; // INT7,                       NCells >=1
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111121060033220000"; mesonCutArray[1] = "0163103100000000"; // INT7,                       NCells >=3
	} else if (trainConfig == 59){  // EMCAL clusters, INT7 trigger
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111121060032000000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1,                                 no M02 cut
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111101060032220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG1, 						 	standard kSDMv5
		eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111122060032220000"; mesonCutArray[2] = "0163103100000000"; // EMCEG1,							NonLinearity LHC11a Calo
		eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111100060032220000"; mesonCutArray[3] = "0163103100000000"; // EMCEG1,							NonLinearity none
		
	} else if (trainConfig == 60){  // EMCAL clusters, EMC7 trigger
		eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000000"; // EMC7, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111121060052220000"; mesonCutArray[1] = "0163103100000000"; // EMC7, 600 MeV min energy
	} else if (trainConfig == 61){  // EMCAL clusters, EMC7 trigger
		eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111121060031220000"; mesonCutArray[0] = "0163103100000000"; // EMC7,                     NCells >=1
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111121060033220000"; mesonCutArray[1] = "0163103100000000"; // EMC7,                     NCells >=3
	} else if (trainConfig == 62){  // EMCAL clusters, EMC7 trigger
		eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111121060032000000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1,                                 no M02 cut
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111101060032220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG1, 						 	standard kSDMv5
		eventCutArray[ 2] = "00052113"; clusterCutArray[2] = "1111122060032220000"; mesonCutArray[2] = "0163103100000000"; // EMCEG1,							NonLinearity LHC11a Calo
		eventCutArray[ 3] = "00052113"; clusterCutArray[3] = "1111100060032220000"; mesonCutArray[3] = "0163103100000000"; // EMCEG1,							NonLinearity none

		
	} else if (trainConfig == 63){  // EMCAL clusters, INT7 trigger
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000000"; // INT7, 400 MeV min energy, NCells >=2, M02 default cut
	} else if (trainConfig == 64){  // no NonLinearity
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111100060032220000"; mesonCutArray[0] = "0163103100000000"; // INT7
		eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111100060032220000"; mesonCutArray[0] = "0163103100000000"; // EMC7
		eventCutArray[ 0] = "00085113"; clusterCutArray[0] = "1111100060032220000"; mesonCutArray[0] = "0163103100000000"; // EG2
		eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111100060032220000"; mesonCutArray[0] = "0163103100000000"; // EG1
    
    // LHC11a cut studies
	} else if (trainConfig == 70){  // EMCAL clusters, MB (INT1) trigger
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000000"; // MB, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050052220000"; mesonCutArray[1] = "0163103100000000"; // MB, 600 MeV min energy
	} else if (trainConfig == 71){  // EMCAL clusters, MB (INT1) trigger
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050031220000"; mesonCutArray[0] = "0163103100000000"; // MB,                     NCells >=1
		eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050033220000"; mesonCutArray[1] = "0163103100000000"; // MB,                     NCells >=3
	} else if (trainConfig == 72){  // EMCAL clusters, MB (INT1) trigger
		eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050032000000"; mesonCutArray[0] = "0163103100000000"; // MB,                                 no M02 cut
		eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121060032220000"; mesonCutArray[1] = "0163103100000000"; // MB,                                                 30ns timing
		eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111121040032220000"; mesonCutArray[2] = "0163103100000000"; // MB,                                                 100ns timing

	} else if (trainConfig == 73){  // EMCAL clusters, EMC1 trigger
		eventCutArray[ 0] = "00051113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000000"; // EMC1, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121050052220000"; mesonCutArray[1] = "0163103100000000"; // EMC1, 600 MeV min energy
	} else if (trainConfig == 74){  // EMCAL clusters, EMC1 trigger
		eventCutArray[ 0] = "00051113"; clusterCutArray[0] = "1111121050031220000"; mesonCutArray[0] = "0163103100000000"; // EMC1,                     NCells >=1
		eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121050033220000"; mesonCutArray[1] = "0163103100000000"; // EMC1,                     NCells >=3
	} else if (trainConfig == 75){  // EMCAL clusters, EMC1 trigger
		eventCutArray[ 0] = "00051113"; clusterCutArray[0] = "1111121050032000000"; mesonCutArray[0] = "0163103100000000"; // EMC1,                                 no M02 cut
		eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121060032220000"; mesonCutArray[1] = "0163103100000000"; // EMC1,                                                 30ns timing
		eventCutArray[ 2] = "00051113"; clusterCutArray[2] = "1111121040032220000"; mesonCutArray[2] = "0163103100000000"; // EMC1,                                                 100ns timing

	} else if (trainConfig == 76){  // EMCAL clusters, MB (INT1) trigger, for added signals
		eventCutArray[ 0] = "00003123"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000000"; // MB, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "00003123"; clusterCutArray[1] = "1111121050052220000"; mesonCutArray[1] = "0163103100000000"; // MB, 600 MeV min energy
	} else if (trainConfig == 77){  // EMCAL clusters, MB (INT1) trigger, for added signals
		eventCutArray[ 0] = "00003123"; clusterCutArray[0] = "1111121050031220000"; mesonCutArray[0] = "0163103100000000"; // MB,                     NCells >=1
		eventCutArray[ 1] = "00003123"; clusterCutArray[1] = "1111121050033220000"; mesonCutArray[1] = "0163103100000000"; // MB,                     NCells >=3
	} else if (trainConfig == 78){  // EMCAL clusters, MB (INT1) trigger, for added signals
		eventCutArray[ 0] = "00003123"; clusterCutArray[0] = "1111121050032000000"; mesonCutArray[0] = "0163103100000000"; // MB,                                 no M02 cut
		eventCutArray[ 1] = "00003123"; clusterCutArray[1] = "1111121020032220000"; mesonCutArray[1] = "0163103100000000"; // MB,                                                 500ns timing
		eventCutArray[ 2] = "00003123"; clusterCutArray[2] = "1111121040032220000"; mesonCutArray[2] = "0163103100000000"; // MB,                                                 100ns timing

	} else if (trainConfig == 79){  // EMCAL clusters, EMC1 special trigger
		eventCutArray[ 0] = "00051123"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000000"; // EMC1, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "00051123"; clusterCutArray[1] = "1111121050052220000"; mesonCutArray[1] = "0163103100000000"; // EMC1, 600 MeV min energy
	} else if (trainConfig == 80){  // EMCAL clusters, EMC1 special trigger
		eventCutArray[ 0] = "00051123"; clusterCutArray[0] = "1111121050031220000"; mesonCutArray[0] = "0163103100000000"; // EMC1,                     NCells >=1
		eventCutArray[ 1] = "00051123"; clusterCutArray[1] = "1111121050033220000"; mesonCutArray[1] = "0163103100000000"; // EMC1,                     NCells >=3
	} else if (trainConfig == 81){  // EMCAL clusters, EMC1 special trigger
		eventCutArray[ 0] = "00051123"; clusterCutArray[0] = "1111121050032000000"; mesonCutArray[0] = "0163103100000000"; // EMC1,                                 no M02 cut
		eventCutArray[ 1] = "00051123"; clusterCutArray[1] = "1111121020032220000"; mesonCutArray[1] = "0163103100000000"; // EMC1,                                                 500ns timing
		eventCutArray[ 2] = "00051123"; clusterCutArray[2] = "1111121040032220000"; mesonCutArray[2] = "0163103100000000"; // EMC1,                                                 100ns timing

  
    // 8 TeV configs
		//standard cut
	} else if (trainConfig == 101){ // EMCAL clusters pp 8 TeV 
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060032230000"; mesonCutArray[0] = "0163103100000000"; // 400 MeV cluster min energy
		// 8 TeV variations
	} else if (trainConfig == 102){ // EMCAL clusters pp 8 TeV, timing variation
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111050032230000"; mesonCutArray[0] = "0163103100000000"; // time -50ns_50ns
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111060032230000"; mesonCutArray[1] = "0163103100000000"; // time -30ns_35ns - standard
		eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111111070032230000"; mesonCutArray[2] = "0163103100000000"; // time -30ns_30ns
		eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111111080032230000"; mesonCutArray[3] = "0163103100000000"; // time -20ns_30ns
	} else if (trainConfig == 103){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060012230000"; mesonCutArray[0] = "0163103100000000"; //0.2 GeV/c
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111060022230000"; mesonCutArray[1] = "0163103100000000"; //0.3 GeV/c
		eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111111060032230000"; mesonCutArray[2] = "0163103100000000"; //0.4 GeV/c default
		eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111111060042230000"; mesonCutArray[3] = "0163103100000000"; //0.5 GeV/c
		eventCutArray[ 4] = "00000113"; clusterCutArray[4] = "1111111060052230000"; mesonCutArray[4] = "0163103100000000"; //0.6 GeV/c
	} else if (trainConfig == 104){ //EMCAL minNCells, M02, with/without TRD variation
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060031230000"; mesonCutArray[0] = "0163103100000000"; //n cells >= 1
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111060033230000"; mesonCutArray[1] = "0163103100000000"; //n cells >= 3
		eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111111060032000000"; mesonCutArray[2] = "0163103100000000"; //no M02 cut
		eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1113111060032230000"; mesonCutArray[3] = "0163103100000000"; //only modules with TRD infront
		eventCutArray[ 4] = "00000113"; clusterCutArray[4] = "1111211060032230000"; mesonCutArray[4] = "0163103100000000"; //no modules with TRD infront

	} else if (trainConfig == 110){ // EMCAL clusters pp 8 TeV, Different NonLinearities
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111101063032230000"; mesonCutArray[0] = "0163103100000000"; // NonLinearity kSDMv5
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000000"; // NonLinearity LHC12 ConvCalo
		eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111112063032230000"; mesonCutArray[2] = "0163103100000000"; // NonLinearity LHC12 Calo
		eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111100063032230000"; mesonCutArray[3] = "0163103100000000"; // NonLinearity none

   // LHC12fa-i and MC
    // default with three cuts
	} else if (trainConfig == 111){  // EMCAL clusters, different triggers, no NonLinearity
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111100060032230000"; mesonCutArray[0] = "0163103100000000"; // INT7
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111100060032230000"; mesonCutArray[1] = "0163103100000000"; // EMC7
		eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111100060032230000"; mesonCutArray[2] = "0163103100000000"; // EMCEG1,
    
    // all the cut variations
	} else if (trainConfig == 112){  // EMCAL clusters, EMCEG1 trigger
		eventCutArray[ 0] = "00081113"; clusterCutArray[0] = "1111111060032230000"; mesonCutArray[0] = "0163103100000000"; // EMCEGA, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "00081113"; clusterCutArray[1] = "1111111060052230000"; mesonCutArray[1] = "0163103100000000"; // EMCEGA, 600 MeV min energy
	} else if (trainConfig == 113){  // EMCAL clusters, EMCEG1 trigger
		eventCutArray[ 0] = "00081113"; clusterCutArray[0] = "1111111060031230000"; mesonCutArray[0] = "0163103100000000"; // EMCEGA,                     NCells >=1
		eventCutArray[ 1] = "00081113"; clusterCutArray[1] = "1111111060033230000"; mesonCutArray[1] = "0163103100000000"; // EMCEGA,                     NCells >=3
	} else if (trainConfig == 114){  // EMCAL clusters, EMCEG1 trigger
		eventCutArray[ 0] = "00081113"; clusterCutArray[0] = "1111111060032000000"; mesonCutArray[0] = "0163103100000000"; // EMCEGA,                                 no M02 cut
		eventCutArray[ 1] = "00081113"; clusterCutArray[1] = "1111111020032230000"; mesonCutArray[1] = "0163103100000000"; // EMCEGA,                                                 500ns timing
		eventCutArray[ 2] = "00085113"; clusterCutArray[2] = "1111111040032230000"; mesonCutArray[2] = "0163103100000000"; // EMCEGA,                                                 100ns timing
	} else if (trainConfig == 115){  // EMCAL clusters, INT7 trigger
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060032230000"; mesonCutArray[0] = "0163103100000000"; // INT7, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111060052230000"; mesonCutArray[1] = "0163103100000000"; // INT7, 600 MeV min energy
	} else if (trainConfig == 116){  // EMCAL clusters, INT7 trigger
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060031230000"; mesonCutArray[0] = "0163103100000000"; // INT7,                       NCells >=1
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111060033230000"; mesonCutArray[1] = "0163103100000000"; // INT7,                       NCells >=3
	} else if (trainConfig == 117){  // EMCAL clusters, INT7 trigger
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060032000000"; mesonCutArray[0] = "0163103100000000"; // INT7,                                   no M02 cut
		eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111020032230000"; mesonCutArray[1] = "0163103100000000"; // INT7,                                                   500ns timing
		eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111111040032230000"; mesonCutArray[2] = "0163103100000000"; // INT7,                                                   100ns timing
	} else if (trainConfig == 118){  // EMCAL clusters, EMC7 trigger
		eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111111060032230000"; mesonCutArray[0] = "0163103100000000"; // EMC7, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111111060052230000"; mesonCutArray[1] = "0163103100000000"; // EMC7, 600 MeV min energy
	} else if (trainConfig == 119){  // EMCAL clusters, EMC7 trigger
		eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111111060031230000"; mesonCutArray[0] = "0163103100000000"; // EMC7,                     NCells >=1
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111111060033230000"; mesonCutArray[1] = "0163103100000000"; // EMC7,                     NCells >=3
	} else if (trainConfig == 120){  // EMCAL clusters, EMC7 trigger
		eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111111060032000000"; mesonCutArray[0] = "0163103100000000"; // EMC7,                                 no M02 cut
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111111020032230000"; mesonCutArray[1] = "0163103100000000"; // EMC7,                                                 500ns timing
		eventCutArray[ 2] = "00052113"; clusterCutArray[2] = "1111111040032230000"; mesonCutArray[2] = "0163103100000000"; // EMC7,                                                 100ns timing

	}else if (trainConfig == 121){ // EMCAL clusters, different special triggers, NonLinearity LHC12 ConvCalo
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060032230000"; mesonCutArray[0] = "0163103100000000"; // INT7
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111111060032230000"; mesonCutArray[1] = "0163103100000000"; // EMC7
		eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111111060032230000"; mesonCutArray[2] = "0163103100000000"; // EMCEG1,
	}else if (trainConfig == 122){ // EMCAL clusters, different special triggers, NonLinearity kSDMv5
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111101060032230000"; mesonCutArray[0] = "0163103100000000"; // INT7
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111101060032230000"; mesonCutArray[1] = "0163103100000000"; // EMC7
		eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111101060032230000"; mesonCutArray[2] = "0163103100000000"; // EMCEG1,
	}else if (trainConfig == 123){ // EMCAL clusters, different special triggers, NonLinearity LHC12 Calo
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111112060032230000"; mesonCutArray[0] = "0163103100000000"; // INT7
		eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111112060032230000"; mesonCutArray[1] = "0163103100000000"; // EMC7
		eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111112060032230000"; mesonCutArray[2] = "0163103100000000"; // EMCEG1,

		// 7 TeV
	} else if (trainConfig == 201){ // EMCAL clusters pp 7 TeV
		eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111100010032230000"; mesonCutArray[0] = "0163103100000000"; // 1000ns timing cut

	} else {
		Error(Form("GammaCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
	TList *ClusterCutList = new TList();
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
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

	for(Int_t i = 0; i<numberOfCuts; i++){
		analysisEventCuts[i] = new AliConvEventCuts();   
		
		// definition of weighting input
		TString fitNamePi0 = Form("Pi0_Fit_Data_%s",energy.Data());
		TString fitNameEta = Form("Eta_Fit_Data_%s",energy.Data());
		TString fAddedSignalString = eventCutArray[i];
		fAddedSignalString = fAddedSignalString(6,1);
		Bool_t fAddedSignal = kFALSE;
		if (fAddedSignalString.CompareTo("2") == 0) fAddedSignal = kTRUE;

		TString mcInputNamePi0 = "";
		TString mcInputNameEta = "";
		if (fAddedSignal && (periodname.Contains("LHC12i3") || periodname.CompareTo("LHC14e2b")==0)){
			mcInputNamePi0 = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
			mcInputNameEta = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
		} else {
			mcInputNamePi0 = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
			mcInputNameEta = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
		}	
		
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
		
		analysisMesonCuts[i] = new AliConversionMesonCuts();
		analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
		MesonCutList->Add(analysisMesonCuts[i]);
		analysisMesonCuts[i]->SetFillCutHistograms("");
		analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
	}
	task->SetEventCutList(numberOfCuts,EventCutList);
	task->SetCaloCutList(numberOfCuts,ClusterCutList);
	task->SetMesonCutList(numberOfCuts,MesonCutList);
	task->SetDoMesonAnalysis(kTRUE);
	task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
	task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
    task->SetDoTHnSparse(isUsingTHnSparse);
	if(enableExtQA == 3){ task->SetPlotHistsExtQA(kTRUE);}
	
	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer(Form("GammaCalo_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaCalo_%i.root",trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
