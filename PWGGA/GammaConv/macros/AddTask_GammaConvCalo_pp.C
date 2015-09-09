void AddTask_GammaConvCalo_pp(  	Int_t 		trainConfig 				= 1,  								//change different set of cuts
									Int_t	 	isMC 						= 0, 								//run MC
									Int_t 		enableQAMesonTask 			= 1, 								//enable QA in AliAnalysisTaskGammaConvV1
									Int_t 		enableQAPhotonTask 			= 1, 								// enable additional QA task
									TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
									TString 	cutnumberAODBranch 			= "000000006008400001001500000",
									Int_t 		enableExtMatchAndQA 		= 0,								// enable matching histograms (1) and extended QA (2), only QA(3), all disabled (0)
									TString 	periodname 					= "LHC12f1x", 						// period name
									Bool_t 		doWeighting 				= kFALSE,							// enables weighting
									Bool_t 		enableV0findingEffi 		= kFALSE,							// enables V0finding efficiency histograms
									Bool_t 		isUsingTHnSparse 			= kTRUE, 							// enable or disable usage of THnSparses for background estimation
								    Bool_t 		enableTriggerMimicking		= kFALSE,							// enable trigger mimicking
									Bool_t 		enableTriggerOverlapRej		= kFALSE,							// enable trigger overlap rejection
									Float_t		maxFacPtHard				= 3.,								// maximum factor between hardest jet and ptHard generated
									TString		periodNameV0Reader			= ""
							) {
	
	Int_t isHeavyIon = 0;
	
	// ================== GetAnalysisManager ===============================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error(Form("AddTask_GammaConvCalo_pp_%i",trainConfig), "No analysis manager found.");
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
		fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);

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
	AliAnalysisTaskGammaConvCalo *task=NULL;
	task= new AliAnalysisTaskGammaConvCalo(Form("GammaConvCalo_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 2;
	if (trainConfig == 101 || 	trainConfig == 131 || 	trainConfig == 120 || 	trainConfig == 121||
		trainConfig == 122 || 	trainConfig == 123) {
		numberOfCuts = 1;
	}
	if (trainConfig == 113 || 	trainConfig == 114 || 	trainConfig == 115 || 	trainConfig == 116 ||	trainConfig == 129) {
		numberOfCuts = 3;
	}
	if (trainConfig == 8 || 	trainConfig == 10 || 	trainConfig == 12 ||	trainConfig == 13 || 	trainConfig == 14 ||
		trainConfig == 21 ||	trainConfig == 23 || 	trainConfig == 24 || 	trainConfig == 26 || 	trainConfig == 27 ||
		trainConfig == 40 ||  	trainConfig == 41 ||  	trainConfig == 48 ||  	trainConfig == 50 || 	trainConfig == 52 ||  	
		trainConfig == 53 ||  	trainConfig == 54 ||  	trainConfig == 61 ||  	trainConfig == 63 || 	trainConfig == 65 || 
		trainConfig == 67 ||  	trainConfig == 68 || 	trainConfig == 75 || 	trainConfig == 77 || 	trainConfig == 79 || 
		trainConfig == 80 ||  	trainConfig == 81 || 	trainConfig == 88 || 	trainConfig == 90 || 	trainConfig == 92 || 
		trainConfig == 93 ||  	trainConfig == 94 || 	trainConfig == 97 ||	trainConfig == 30 ||
		trainConfig == 108 || 	trainConfig == 126 ||	trainConfig == 128 ||
		trainConfig == 111 || 	trainConfig == 117 || 	trainConfig == 118 || 	trainConfig == 119||	trainConfig == 23) {
		numberOfCuts = 4;		
	}
	if (trainConfig == 2 || 	trainConfig == 3 || 	trainConfig == 5 || 	trainConfig == 6 || 	trainConfig == 7 || 
		trainConfig == 11 ||	trainConfig == 15 ||	trainConfig == 16 ||	trainConfig == 18 ||	trainConfig == 19 ||
		trainConfig == 20 ||	trainConfig == 25 ||	trainConfig == 42 ||	trainConfig == 43 ||	trainConfig == 45 ||	
		trainConfig == 46 ||	trainConfig == 47 ||	trainConfig == 51 ||	trainConfig == 55 ||	trainConfig == 56 ||	
		trainConfig == 58 ||	trainConfig == 59 ||	trainConfig == 60 ||	trainConfig == 64 ||	trainConfig == 69 ||	
		trainConfig == 70 ||	trainConfig == 72 ||	trainConfig == 73 ||	trainConfig == 74 ||	trainConfig == 78 ||	
		trainConfig == 82 ||	trainConfig == 83 ||	trainConfig == 85 ||	trainConfig == 86 ||	trainConfig == 87 ||	
		trainConfig == 91 ||
		trainConfig == 102 ||	trainConfig == 103 || 	trainConfig == 105 || 	trainConfig == 106 ||
		trainConfig == 107 ||	trainConfig == 127) {
		numberOfCuts = 5;	
	}
	if (trainConfig == 4  || 	trainConfig == 17 || 	trainConfig == 44 || 	trainConfig == 57 || 	 trainConfig == 71 ||
		trainConfig == 84 ||	trainConfig == 301 ||	trainConfig == 302 ||
		trainConfig == 104 ||	trainConfig == 110) {
		numberOfCuts = 6;
	}

	TString *eventCutArray = new TString[numberOfCuts];
	TString *photonCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// cluster cuts
	// 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
	// 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"

	// ************************************* EMCAL cuts ****************************************************
	// LHC11a with new non linearities
	if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // 400 MeV cluster min energy
		
	// minimum bias variations	
	} else if (trainConfig == 2){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121053012230000"; mesonCutArray[0] = "0163103100000010"; //0.2 GeV/c
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121053022230000"; mesonCutArray[1] = "0163103100000010"; //0.3 GeV/c
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; //0.4 GeV/c default
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121053042230000"; mesonCutArray[3] = "0163103100000010"; //0.5 GeV/c
		eventCutArray[ 4] = "00003113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121053052230000"; mesonCutArray[4] = "0163103100000010"; //0.6 GeV/c
	} else if (trainConfig == 3){ //EMCAL minNCells variation
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121053031230000"; mesonCutArray[0] = "0163103100000010"; //n cells >= 1
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121053033230000"; mesonCutArray[1] = "0163103100000010"; //n cells >= 3
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121053032200000"; mesonCutArray[2] = "0163103100000010"; //no max M02 cut
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1113121053032230000"; mesonCutArray[3] = "0163103100000010"; //only modules with TRD infront
		eventCutArray[ 4] = "00003113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111221053032230000"; mesonCutArray[4] = "0163103100000010"; //no modules with TRD infront
	} else if (trainConfig == 4){ // EMCAL track matching variations 
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121051032230000"; mesonCutArray[0] = "0163103100000010"; //
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121052032230000"; mesonCutArray[1] = "0163103100000010"; //
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; //
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121054032230000"; mesonCutArray[3] = "0163103100000010"; //
		eventCutArray[ 4] = "00003113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121055032230000"; mesonCutArray[4] = "0163103100000010"; //
		eventCutArray[ 5] = "00003113"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "1111121056032230000"; mesonCutArray[5] = "0163103100000010"; //
	} else if (trainConfig == 5){ // PCM variations
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009227000008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx e -3, 5
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009127000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx e -5, 5
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009357000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi 2
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009317000008250400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi 0
		eventCutArray[ 4] = "00003113"; photonCutArray[ 4] = "00200009387300008250400000"; clusterCutArray[4] = "1111121053032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi 2 high 1 (> 3.5 GeV)
	} else if (trainConfig == 6){ // PCM variations
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000009250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // qt 2D 0.03
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327000003250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // qt 1D 0.05
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009327000002250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // qt 1D 0.07
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200049327000008250400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // single pt > 0.075
		eventCutArray[ 4] = "00003113"; photonCutArray[ 4] = "00200019327000008250400000"; clusterCutArray[4] = "1111121053032230000"; mesonCutArray[4] = "0163103100000010"; // single pt > 0.1
	} else if (trainConfig == 7){ // PCM variations
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008850400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327000008260400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009327000008860400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009327000008280400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 4] = "00003113"; photonCutArray[ 4] = "00200009327000008880400000"; clusterCutArray[4] = "1111121053032230000"; mesonCutArray[4] = "0163103100000010"; // 2D psi pair chi2 var
	} else if (trainConfig == 8){ // PCM variations
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200006327000008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // min TPC cl > 0.7
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200008327000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // min TPC cl > 0.35
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163106100000010"; // alpha < 0.8
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163105100000010"; // alpha < 0.75
	} else if (trainConfig == 9){ // PCM variations
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00202209327000008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // restrict acceptance to EMCAL loose
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00204409327000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // restrict acceptance to EMCAL tight
	} else if (trainConfig == 10){ // PCM variations pi dEdx
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009317300008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327300008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009325000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 1: 0.3 ->
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009320000008250400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 1: 0.5 ->
	} else if (trainConfig == 11){ // PCM variations pi dEdx	
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327600008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327400008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009315600008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009367400008250400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
		eventCutArray[ 4] = "00003113"; photonCutArray[ 4] = "00200009347400008250400000"; clusterCutArray[4] = "1111121053032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi: 3: 0.4-3, 1: 3. ->
	} else if (trainConfig == 12){ // PCM variations to close V0s	
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327000008250401000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009327000008250402000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009327000008250403000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
	} else if (trainConfig == 13){ // EMCAL clusters 2.76 TeV LHC11a, timing variation
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // time 50ns
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121043032230000"; mesonCutArray[1] = "0163103100000010"; // time 100ns
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121033032230000"; mesonCutArray[2] = "0163103100000010"; // time 200ns
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121023032230000"; mesonCutArray[3] = "0163103100000010"; // time 500ns
	}else if (trainConfig == 14){  //LHC11a NonLinearity variations
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111100053032230000"; mesonCutArray[0] = "0163103100000010"; // NonLinearity none
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111101053032230000"; mesonCutArray[1] = "0163103100000010"; // NonLinearity kSDMv5
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // NonLinearity LHC11a ConvCalo
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111122053032230000"; mesonCutArray[3] = "0163103100000010"; // NonLinearity LHC11a Calo

	// EMC1 variations	
	} else if (trainConfig == 15){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121053012230000"; mesonCutArray[0] = "0163103100000010"; //0.2 GeV/c
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121053022230000"; mesonCutArray[1] = "0163103100000010"; //0.3 GeV/c
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; //0.4 GeV/c default
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121053042230000"; mesonCutArray[3] = "0163103100000010"; //0.5 GeV/c
		eventCutArray[ 4] = "00051113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121053052230000"; mesonCutArray[4] = "0163103100000010"; //0.6 GeV/c
	} else if (trainConfig == 16){ //EMCAL minNCells variation
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121053031230000"; mesonCutArray[0] = "0163103100000010"; //n cells >= 1
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121053033230000"; mesonCutArray[1] = "0163103100000010"; //n cells >= 3
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121053032200000"; mesonCutArray[2] = "0163103100000010"; //no max M02 cut
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1113121053032230000"; mesonCutArray[3] = "0163103100000010"; //only modules with TRD infront
		eventCutArray[ 4] = "00051113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111221053032230000"; mesonCutArray[4] = "0163103100000010"; //no modules with TRD infront
	} else if (trainConfig == 17){ // EMCAL track matching variations 
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121051032230000"; mesonCutArray[0] = "0163103100000010"; //
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121052032230000"; mesonCutArray[1] = "0163103100000010"; //
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; //
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121054032230000"; mesonCutArray[3] = "0163103100000010"; //
		eventCutArray[ 4] = "00051113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121055032230000"; mesonCutArray[4] = "0163103100000010"; //
		eventCutArray[ 5] = "00051113"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "1111121056032230000"; mesonCutArray[5] = "0163103100000010"; //
	} else if (trainConfig == 18){ // PCM variations
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009227000008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx e -3, 5
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009127000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx e -5, 5
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009357000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi 2
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200009317000008250400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi 0
		eventCutArray[ 4] = "00051113"; photonCutArray[ 4] = "00200009387300008250400000"; clusterCutArray[4] = "1111121053032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi 2 high 1 (> 3.5 GeV)
	} else if (trainConfig == 19){ // PCM variations
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009327000009250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // qt 2D 0.03
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327000003250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // qt 1D 0.05
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009327000002250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // qt 1D 0.07
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200049327000008250400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // single pt > 0.075
		eventCutArray[ 4] = "00051113"; photonCutArray[ 4] = "00200019327000008250400000"; clusterCutArray[4] = "1111121053032230000"; mesonCutArray[4] = "0163103100000010"; // single pt > 0.1
	} else if (trainConfig == 20){ // PCM variations
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009327000008850400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327000008260400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009327000008860400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200009327000008280400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 4] = "00051113"; photonCutArray[ 4] = "00200009327000008880400000"; clusterCutArray[4] = "1111121053032230000"; mesonCutArray[4] = "0163103100000010"; // 2D psi pair chi2 var
	} else if (trainConfig == 21){ // PCM variations
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200006327000008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // min TPC cl > 0.7
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200008327000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // min TPC cl > 0.35
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163106100000010"; // alpha < 0.8
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163105100000010"; // alpha < 0.75
	} else if (trainConfig == 22){ // PCM variations
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00202209327000008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // restrict acceptance to EMCAL loose
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00204409327000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // restrict acceptance to EMCAL tight
	} else if (trainConfig == 23){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121053062230000"; mesonCutArray[0] = "0163103100000010"; // min Energy cluster = 4.5 GeV
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121053072230000"; mesonCutArray[1] = "0163103100000010"; // min Energy cluster = 5.0 GeV
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121053082230000"; mesonCutArray[2] = "0163103100000010"; // min Energy cluster = 5.5 GeV
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121053092230000"; mesonCutArray[3] = "0163103100000010"; // min Energy cluster = 6.0 GeV
	} else if (trainConfig == 24){ // PCM variations pi dEdx
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009317300008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327300008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009325000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 1: 0.3 ->
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200009320000008250400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 1: 0.5 ->
	} else if (trainConfig == 25){ // PCM variations pi dEdx	
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009327600008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327400008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009315600008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200009367400008250400000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
		eventCutArray[ 4] = "00051113"; photonCutArray[ 4] = "00200009347400008250400000"; clusterCutArray[4] = "1111121053032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi: 3: 0.4-3, 1: 3. ->
	} else if (trainConfig == 26){ // PCM variations to close V0s	
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327000008250401000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009327000008250402000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200009327000008250403000"; clusterCutArray[3] = "1111121053032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
	} else if (trainConfig == 27){ // LHC11a NonLinearity variations
		eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111100053032230000"; mesonCutArray[0] = "0163103100000010"; // NonLinearity none
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111101053032230000"; mesonCutArray[1] = "0163103100000010"; // NonLinearity kSDMv5
		eventCutArray[ 2] = "00051113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121053032230000"; mesonCutArray[2] = "0163103100000010"; // NonLinearity LHC11a ConvCalo
		eventCutArray[ 3] = "00051113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111122053032230000"; mesonCutArray[3] = "0163103100000010"; // NonLinearity LHC11a Calo
		
	//LHC11a EMCal no non linearity internally	
	} else if (trainConfig == 28){ 
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111100053032230000"; mesonCutArray[0] = "0163103100000010"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111100053032230000"; mesonCutArray[1] = "0163103100000010"; // 400 MeV cluster min energy

	// LHC11a for Gustavo with SetDoTreeConvGammaShowerShape
	} else if (trainConfig == 29){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121053032230000"; mesonCutArray[0] = "0163103100000010"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "00051113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // 400 MeV cluster min energy
	}else if (trainConfig == 30){  //LHC11a additional NonLinearity variations
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111123053032230000"; mesonCutArray[0] = "0163103100000010"; // NonLinearity kTestBeamv2 + ConvCalo
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111124053032230000"; mesonCutArray[1] = "0163103100000010"; // NonLinearity kTestBeamv2 + Calo
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111125053032230000"; mesonCutArray[2] = "0163103100000010"; // NonLinearity LHC11a ConvCalo kSDM
		eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111126053032230000"; mesonCutArray[3] = "0163103100000010"; // NonLinearity LHC11a Calo kSDM

	// LHC13g	
	} else if (trainConfig == 40){  // LHC13g without non linearity
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111100063032230000"; mesonCutArray[0] = "0163103100000010"; // INT7
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111100063032230000"; mesonCutArray[1] = "0163103100000010"; // EMC7
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111100063032230000"; mesonCutArray[2] = "0163103100000010"; // EMCEG1,
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111100063032230000"; mesonCutArray[3] = "0163103100000010"; // EMCEG2,
		
	// LHC13g new conv calo non lienarity
	} else if (trainConfig == 41){  // EMCal, all triggers
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // INT7
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // EMC7
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // EMCEG1,
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // EMCEG2,

	// LHC13g variations MB cut	
	} else if (trainConfig == 42){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063012230000"; mesonCutArray[0] = "0163103100000010"; //0.2 GeV/c
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063022230000"; mesonCutArray[1] = "0163103100000010"; //0.3 GeV/c
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; //0.4 GeV/c default
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121063042230000"; mesonCutArray[3] = "0163103100000010"; //0.5 GeV/c
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121063052230000"; mesonCutArray[4] = "0163103100000010"; //0.6 GeV/c
	} else if (trainConfig == 43){ //EMCAL minNCells variation
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063031230000"; mesonCutArray[0] = "0163103100000010"; //n cells >= 1
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063033230000"; mesonCutArray[1] = "0163103100000010"; //n cells >= 3
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032200000"; mesonCutArray[2] = "0163103100000010"; //no max M02 cut
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1112121063032230000"; mesonCutArray[3] = "0163103100000010"; //only modules with TRD infront
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111321063032230000"; mesonCutArray[4] = "0163103100000010"; //no modules with TRD infront
	} else if (trainConfig == 44){ // EMCAL track matching variations 
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121061032230000"; mesonCutArray[0] = "0163103100000010"; //
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121062032230000"; mesonCutArray[1] = "0163103100000010"; //
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; //
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121064032230000"; mesonCutArray[3] = "0163103100000010"; //
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121065032230000"; mesonCutArray[4] = "0163103100000010"; //
		eventCutArray[ 5] = "00000113"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "1111121066032230000"; mesonCutArray[5] = "0163103100000010"; //
	} else if (trainConfig == 45){ // PCM variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx e -3, 5
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009127000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx e -5, 5
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009357000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi 2
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009317000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi 0
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009387300008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi 2 high 1 (> 3.5 GeV)
	} else if (trainConfig == 46){ // PCM variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000009250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // qt 2D 0.03
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000003250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // qt 1D 0.05
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000002250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // qt 1D 0.07
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200049327000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // single pt > 0.075
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200019327000008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // single pt > 0.1
	} else if (trainConfig == 47){ // PCM variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008850400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008260400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008860400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008280400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009327000008880400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // 2D psi pair chi2 var
	} else if (trainConfig == 48){ // PCM variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200006327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // min TPC cl > 0.7
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200008327000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // min TPC cl > 0.35
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163106100000010"; // alpha < 0.8
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163105100000010"; // alpha < 0.75
	} else if (trainConfig == 49){ // PCM variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00202209327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // restrict acceptance to EMCAL loose
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00204409327000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // restrict acceptance to EMCAL tight
	} else if (trainConfig == 50){ // PCM variations pi dEdx
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009317300008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327300008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009325000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 1: 0.3 ->
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009320000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 1: 0.5 ->
	} else if (trainConfig == 51){ // PCM variations pi dEdx	
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327600008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327400008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009315600008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009367400008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009347400008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi: 3: 0.4-3, 1: 3. ->
	} else if (trainConfig == 52){ // PCM variations to close V0s	
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250401000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250402000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250403000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
	} else if (trainConfig == 53){ // EMCAL clusters 2.76 TeV LHC11a, timing variation
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // time 30ns
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // time 50ns
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121043032230000"; mesonCutArray[2] = "0163103100000010"; // time 100ns
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121033032230000"; mesonCutArray[3] = "0163103100000010"; // time 200ns
	} else if (trainConfig == 54){  //LHC11a NonLinearity variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111100063032230000"; mesonCutArray[0] = "0163103100000010"; // NonLinearity none
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111101063032230000"; mesonCutArray[1] = "0163103100000010"; // NonLinearity kSDMv5
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // NonLinearity LHC11a ConvCalo
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111122063032230000"; mesonCutArray[3] = "0163103100000010"; // NonLinearity LHC11a Calo
		
	// LHC13g variations EMC7 cut	
	} else if (trainConfig == 55){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063012230000"; mesonCutArray[0] = "0163103100000010"; //0.2 GeV/c
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063022230000"; mesonCutArray[1] = "0163103100000010"; //0.3 GeV/c
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; //0.4 GeV/c default
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121063042230000"; mesonCutArray[3] = "0163103100000010"; //0.5 GeV/c
		eventCutArray[ 4] = "00052113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121063052230000"; mesonCutArray[4] = "0163103100000010"; //0.6 GeV/c
	} else if (trainConfig == 56){ //EMCAL minNCells variation
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063031230000"; mesonCutArray[0] = "0163103100000010"; //n cells >= 1
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063033230000"; mesonCutArray[1] = "0163103100000010"; //n cells >= 3
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032200000"; mesonCutArray[2] = "0163103100000010"; //no max M02 cut
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1112121063032230000"; mesonCutArray[3] = "0163103100000010"; //only modules with TRD infront
		eventCutArray[ 4] = "00052113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111321063032230000"; mesonCutArray[4] = "0163103100000010"; //no modules with TRD infront
	} else if (trainConfig == 57){ // EMCAL track matching variations 
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121061032230000"; mesonCutArray[0] = "0163103100000010"; //
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121062032230000"; mesonCutArray[1] = "0163103100000010"; //
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; //
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121064032230000"; mesonCutArray[3] = "0163103100000010"; //
		eventCutArray[ 4] = "00052113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121065032230000"; mesonCutArray[4] = "0163103100000010"; //
		eventCutArray[ 5] = "00052113"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "1111121066032230000"; mesonCutArray[5] = "0163103100000010"; //
	} else if (trainConfig == 58){ // PCM variations
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009227000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx e -3, 5
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009127000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx e -5, 5
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009357000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi 2
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200009317000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi 0
		eventCutArray[ 4] = "00052113"; photonCutArray[ 4] = "00200009387300008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi 2 high 1 (> 3.5 GeV)
	} else if (trainConfig == 59){ // PCM variations
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009327000009250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // qt 2D 0.03
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000003250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // qt 1D 0.05
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009327000002250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // qt 1D 0.07
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200049327000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // single pt > 0.075
		eventCutArray[ 4] = "00052113"; photonCutArray[ 4] = "00200019327000008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // single pt > 0.1
	} else if (trainConfig == 60){ // PCM variations
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009327000008850400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008260400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009327000008860400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200009327000008280400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 4] = "00052113"; photonCutArray[ 4] = "00200009327000008880400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // 2D psi pair chi2 var
	} else if (trainConfig == 61){ // PCM variations
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200006327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // min TPC cl > 0.7
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200008327000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // min TPC cl > 0.35
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163106100000010"; // alpha < 0.8
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163105100000010"; // alpha < 0.75
	} else if (trainConfig == 62){ // PCM variations
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00202209327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // restrict acceptance to EMCAL loose
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00204409327000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // restrict acceptance to EMCAL tight
	} else if (trainConfig == 63){ // PCM variations pi dEdx
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009317300008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327300008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009325000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 1: 0.3 ->
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200009320000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 1: 0.5 ->
	} else if (trainConfig == 64){ // PCM variations pi dEdx	
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009327600008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327400008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009315600008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200009367400008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
		eventCutArray[ 4] = "00052113"; photonCutArray[ 4] = "00200009347400008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi: 3: 0.4-3, 1: 3. ->
	} else if (trainConfig == 65){ // PCM variations to close V0s	
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250401000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009327000008250402000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200009327000008250403000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
	} else if (trainConfig == 67){ // EMCAL clusters 2.76 TeV LHC11a, timing variation
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // time 50ns
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // time 100ns
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121043032230000"; mesonCutArray[2] = "0163103100000010"; // time 200ns
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121033032230000"; mesonCutArray[3] = "0163103100000010"; // time 500ns
	} else if (trainConfig == 68){  //LHC11a NonLinearity variations
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111100063032230000"; mesonCutArray[0] = "0163103100000010"; // NonLinearity none
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111101063032230000"; mesonCutArray[1] = "0163103100000010"; // NonLinearity kSDMv5
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // NonLinearity LHC11a ConvCalo
		eventCutArray[ 3] = "00052113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111122063032230000"; mesonCutArray[3] = "0163103100000010"; // NonLinearity LHC11a Calo
		
	// LHC13g variations EG2 cut	
	} else if (trainConfig == 69){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063012230000"; mesonCutArray[0] = "0163103100000010"; //0.2 GeV/c
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063022230000"; mesonCutArray[1] = "0163103100000010"; //0.3 GeV/c
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; //0.4 GeV/c default
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121063042230000"; mesonCutArray[3] = "0163103100000010"; //0.5 GeV/c
		eventCutArray[ 4] = "00085113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121063052230000"; mesonCutArray[4] = "0163103100000010"; //0.6 GeV/c
	} else if (trainConfig == 70){ //EMCAL minNCells variation
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063031230000"; mesonCutArray[0] = "0163103100000010"; //n cells >= 1
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063033230000"; mesonCutArray[1] = "0163103100000010"; //n cells >= 3
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032200000"; mesonCutArray[2] = "0163103100000010"; //no max M02 cut
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1113121063032230000"; mesonCutArray[3] = "0163103100000010"; //only modules with TRD infront
		eventCutArray[ 4] = "00085113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111321063032230000"; mesonCutArray[4] = "0163103100000010"; //no modules with TRD infront
	} else if (trainConfig == 71){ // EMCAL track matching variations 
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121061032230000"; mesonCutArray[0] = "0163103100000010"; //
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121062032230000"; mesonCutArray[1] = "0163103100000010"; //
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; //
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121064032230000"; mesonCutArray[3] = "0163103100000010"; //
		eventCutArray[ 4] = "00085113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121065032230000"; mesonCutArray[4] = "0163103100000010"; //
		eventCutArray[ 5] = "00085113"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "1111121066032230000"; mesonCutArray[5] = "0163103100000010"; //
	} else if (trainConfig == 72){ // PCM variations
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009227000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx e -3, 5
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009127000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx e -5, 5
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009357000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi 2
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009317000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi 0
		eventCutArray[ 4] = "00085113"; photonCutArray[ 4] = "00200009387300008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi 2 high 1 (> 3.5 GeV)
	} else if (trainConfig == 73){ // PCM variations
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009327000009250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // qt 2D 0.03
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009327000003250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // qt 1D 0.05
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009327000002250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // qt 1D 0.07
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200049327000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // single pt > 0.075
		eventCutArray[ 4] = "00085113"; photonCutArray[ 4] = "00200019327000008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // single pt > 0.1
	} else if (trainConfig == 74){ // PCM variations
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009327000008850400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009327000008260400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009327000008860400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009327000008280400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 4] = "00085113"; photonCutArray[ 4] = "00200009327000008880400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // 2D psi pair chi2 var
	} else if (trainConfig == 75){ // PCM variations
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200006327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // min TPC cl > 0.7
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200008327000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // min TPC cl > 0.35
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163106100000010"; // alpha < 0.8
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163105100000010"; // alpha < 0.75
	} else if (trainConfig == 76){ // PCM variations
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00202209327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // restrict acceptance to EMCAL loose
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00204409327000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // restrict acceptance to EMCAL tight
	} else if (trainConfig == 77){ // PCM variations pi dEdx
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009317300008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009327300008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009325000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 1: 0.3 ->
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009320000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 1: 0.5 ->
	} else if (trainConfig == 78){ // PCM variations pi dEdx	
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009327600008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009327400008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009315600008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009367400008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
		eventCutArray[ 4] = "00085113"; photonCutArray[ 4] = "00200009347400008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi: 3: 0.4-3, 1: 3. ->
	} else if (trainConfig == 79){ // PCM variations to close V0s	
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009327000008250401000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009327000008250402000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009327000008250403000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
	} else if (trainConfig == 80){ // EMCAL clusters 2.76 TeV LHC11a, timing variation
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // time 50ns
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // time 100ns
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121043032230000"; mesonCutArray[2] = "0163103100000010"; // time 200ns
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121033032230000"; mesonCutArray[3] = "0163103100000010"; // time 500ns
	} else if (trainConfig == 81){  //LHC11a NonLinearity variations
		eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111100063032230000"; mesonCutArray[0] = "0163103100000010"; // NonLinearity none
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111101063032230000"; mesonCutArray[1] = "0163103100000010"; // NonLinearity kSDMv5
		eventCutArray[ 2] = "00085113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // NonLinearity LHC11a ConvCalo
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111122063032230000"; mesonCutArray[3] = "0163103100000010"; // NonLinearity LHC11a Calo

	// LHC13g variations EG2 cut	
	} else if (trainConfig == 82){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063012230000"; mesonCutArray[0] = "0163103100000010"; //0.2 GeV/c
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063022230000"; mesonCutArray[1] = "0163103100000010"; //0.3 GeV/c
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; //0.4 GeV/c default
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121063042230000"; mesonCutArray[3] = "0163103100000010"; //0.5 GeV/c
		eventCutArray[ 4] = "00083113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121063052230000"; mesonCutArray[4] = "0163103100000010"; //0.6 GeV/c
	} else if (trainConfig == 83){ //EMCAL minNCells variation
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063031230000"; mesonCutArray[0] = "0163103100000010"; //n cells >= 1
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063033230000"; mesonCutArray[1] = "0163103100000010"; //n cells >= 3
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032200000"; mesonCutArray[2] = "0163103100000010"; //no max M02 cut
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1112121063032230000"; mesonCutArray[3] = "0163103100000010"; //only modules with TRD infront
		eventCutArray[ 4] = "00083113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111321063032230000"; mesonCutArray[4] = "0163103100000010"; //no modules with TRD infront
	} else if (trainConfig == 84){ // EMCAL track matching variations 
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121061032230000"; mesonCutArray[0] = "0163103100000010"; //
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121062032230000"; mesonCutArray[1] = "0163103100000010"; //
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; //
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121064032230000"; mesonCutArray[3] = "0163103100000010"; //
		eventCutArray[ 4] = "00083113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111121065032230000"; mesonCutArray[4] = "0163103100000010"; //
		eventCutArray[ 5] = "00083113"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "1111121066032230000"; mesonCutArray[5] = "0163103100000010"; //
	} else if (trainConfig == 85){ // PCM variations
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009227000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx e -3, 5
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200009127000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx e -5, 5
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009357000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi 2
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200009317000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi 0
		eventCutArray[ 4] = "00083113"; photonCutArray[ 4] = "00200009387300008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi 2 high 1 (> 3.5 GeV)
	} else if (trainConfig == 86){ // PCM variations
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009327000009250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // qt 2D 0.03
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200009327000003250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // qt 1D 0.05
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000002250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // qt 1D 0.07
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200049327000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // single pt > 0.075
		eventCutArray[ 4] = "00083113"; photonCutArray[ 4] = "00200019327000008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // single pt > 0.1
	} else if (trainConfig == 87){ // PCM variations
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009327000008850400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200009327000008260400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000008860400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200009327000008280400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 4] = "00083113"; photonCutArray[ 4] = "00200009327000008880400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // 2D psi pair chi2 var
	} else if (trainConfig == 88){ // PCM variations
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200006327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // min TPC cl > 0.7
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200008327000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // min TPC cl > 0.35
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163106100000010"; // alpha < 0.8
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163105100000010"; // alpha < 0.75
	} else if (trainConfig == 89){ // PCM variations
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00202209327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // restrict acceptance to EMCAL loose
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00204409327000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // restrict acceptance to EMCAL tight
	} else if (trainConfig == 90){ // PCM variations pi dEdx
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009317300008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200009327300008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009325000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 1: 0.3 ->
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200009320000008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 1: 0.5 ->
	} else if (trainConfig == 91){ // PCM variations pi dEdx	
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009327600008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200009327400008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009315600008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200009367400008250400000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
		eventCutArray[ 4] = "00083113"; photonCutArray[ 4] = "00200009347400008250400000"; clusterCutArray[4] = "1111121063032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi: 3: 0.4-3, 1: 3. ->
	} else if (trainConfig == 92){ // PCM variations to close V0s	
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200009327000008250401000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000008250402000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200009327000008250403000"; clusterCutArray[3] = "1111121063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
	} else if (trainConfig == 93){ // EMCAL clusters 2.76 TeV LHC11a, timing variation
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // time 50ns
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121053032230000"; mesonCutArray[1] = "0163103100000010"; // time 100ns
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121043032230000"; mesonCutArray[2] = "0163103100000010"; // time 200ns
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121033032230000"; mesonCutArray[3] = "0163103100000010"; // time 500ns
	} else if (trainConfig == 94){  //LHC11a NonLinearity variations
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111100063032230000"; mesonCutArray[0] = "0163103100000010"; // NonLinearity none
		eventCutArray[ 1] = "00083113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111101063032230000"; mesonCutArray[1] = "0163103100000010"; // NonLinearity kSDMv5
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032230000"; mesonCutArray[2] = "0163103100000010"; // NonLinearity LHC11a ConvCalo
		eventCutArray[ 3] = "00083113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111122063032230000"; mesonCutArray[3] = "0163103100000010"; // NonLinearity LHC11a Calo
		
	//LHC13g
	} else if (trainConfig == 95){  // EMCAL clusters, EMCEGA triggers, track matching 0.035
		eventCutArray[ 0] = "00083113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // EMCEG1,
		eventCutArray[ 1] = "00085113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // EMCEG2,
	} else if (trainConfig == 96){  // EMCAL clusters, kEMC trigger, track matching 0.035
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032230000"; mesonCutArray[0] = "0163103100000010"; // INT7
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063032230000"; mesonCutArray[1] = "0163103100000010"; // EMC7
	} else if (trainConfig == 97){  // EMCAL clusters, EMCEGA triggers, track matching 0.035
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111121063032200000"; mesonCutArray[0] = "0163103100000010"; // INT7, max M02 off
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111121063032200000"; mesonCutArray[1] = "0163103100000010"; // EMC7, max M02 off
		eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111121063032200000"; mesonCutArray[2] = "0163103100000010"; // EMCEG1, max M02 off
		eventCutArray[ 3] = "00085113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111121063032200000"; mesonCutArray[3] = "0163103100000010"; // EMCEG2, max M02 off


	// ************************************* EMCAL cuts ****************************************************
	// LHC12
	} else if (trainConfig == 101){ // EMCAL clusters 8 TeV LHC12
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // 400 MeV cluster min energy
	} else if (trainConfig == 102){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063012230000"; mesonCutArray[0] = "0163103100000010"; //0.2 GeV/c
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111063022230000"; mesonCutArray[1] = "0163103100000010"; //0.3 GeV/c
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000010"; //0.4 GeV/c default
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111111063042230000"; mesonCutArray[3] = "0163103100000010"; //0.5 GeV/c
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111111063052230000"; mesonCutArray[4] = "0163103100000010"; //0.6 GeV/c
	} else if (trainConfig == 103){ //EMCAL minNCells variation
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063031230000"; mesonCutArray[0] = "0163103100000010"; //n cells >= 1
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111063033230000"; mesonCutArray[1] = "0163103100000010"; //n cells >= 3
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111063032200000"; mesonCutArray[2] = "0163103100000010"; //no max M02 cut
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1112111063032230000"; mesonCutArray[3] = "0163103100000010"; //only modules with TRD infront
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111311063032230000"; mesonCutArray[4] = "0163103100000010"; //no modules with TRD infront
	} else if (trainConfig == 104){ // EMCAL track matching variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111061032230000"; mesonCutArray[0] = "0163103100000010"; // track matching variations
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111062032230000"; mesonCutArray[1] = "0163103100000010"; //
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000010"; //
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111111064032230000"; mesonCutArray[3] = "0163103100000010"; //
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111111065032230000"; mesonCutArray[4] = "0163103100000010"; //
		eventCutArray[ 5] = "00000113"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "1111111066032230000"; mesonCutArray[5] = "0163103100000010"; //
	} else if (trainConfig == 105){ // PCM variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx e -3, 5
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009127000008250400000"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx e -5, 5
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009357000008250400000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi 2
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009317000008250400000"; clusterCutArray[3] = "1111111063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi 0
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009387300008250400000"; clusterCutArray[4] = "1111111063032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi 2 high 1 (> 3.5 GeV)
	} else if (trainConfig == 106){ // PCM variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000009250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // qt 2D 0.03
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000003250400000"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000010"; // qt 1D 0.05
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000002250400000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000010"; // qt 1D 0.07
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200049327000008250400000"; clusterCutArray[3] = "1111111063032230000"; mesonCutArray[3] = "0163103100000010"; // single pt > 0.075
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200019327000008250400000"; clusterCutArray[4] = "1111111063032230000"; mesonCutArray[4] = "0163103100000010"; // single pt > 0.1
	} else if (trainConfig == 107){ // PCM variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008850400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008260400000"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008860400000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008280400000"; clusterCutArray[3] = "1111111063032230000"; mesonCutArray[3] = "0163103100000010"; // 2D psi pair chi2 var
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009327000008880400000"; clusterCutArray[4] = "1111111063032230000"; mesonCutArray[4] = "0163103100000010"; // 2D psi pair chi2 var
	} else if (trainConfig == 108){ // PCM variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200006327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // min TPC cl > 0.7
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200008327000008250400000"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000010"; // min TPC cl > 0.35
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163106100000010"; // alpha < 0.8
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111111063032230000"; mesonCutArray[3] = "0163105100000010"; // alpha < 0.75
	} else if (trainConfig == 109){ // PCM variations
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00202209327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // restrict acceptance to EMCAL loose
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00204409327000008250400000"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000010"; // restrict acceptance to EMCAL tight
	} else if (trainConfig == 110){ // Different NonLinearities
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111100063032230000"; mesonCutArray[0] = "0163103100000010"; // NonLinearity none
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111101063032230000"; mesonCutArray[1] = "0163103100000010"; // NonLinearity kSDMv5
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000010"; // NonLinearity LHC12 ConvCalo
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111112063032230000"; mesonCutArray[3] = "0163103100000010"; // NonLinearity LHC12 Calo
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "1111113063032230000"; mesonCutArray[4] = "0163103100000010"; // NonLinearity LHC12 kTestBeamv2+ConvCalo
		eventCutArray[ 5] = "00000113"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "1111114063032230000"; mesonCutArray[5] = "0163103100000010"; // NonLinearity LHC12 kTestBeamv2+Calo
	} else if (trainConfig == 111){  // EMCAL clusters, EMC triggers (EMC7, EMCEGA, EMCEJE)
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // INT7
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000010"; // EMC7
		eventCutArray[ 2] = "00081113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000010"; // EMCEGA
		eventCutArray[ 3] = "00091113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111111063032230000"; mesonCutArray[3] = "0163103100000010"; // EMCEJE
	} else if (trainConfig == 112){ // With/without Added Signals
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; //
		eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000010"; //
	} else if (trainConfig == 113){  // EMCAL clusters, MB with minEnergy variation
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // min Energy cluster = 0.4 GeV
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111063062230000"; mesonCutArray[1] = "0163103100000010"; // min Energy cluster = 4.5 GeV
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111063092230000"; mesonCutArray[2] = "0163103100000010"; // min Energy cluster = 6.0 GeV
	} else if (trainConfig == 114){ // EMCAL clusters, EMC7 with minEnergy variation
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // min Energy cluster = 0.4 GeV
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111063062230000"; mesonCutArray[1] = "0163103100000010"; // min Energy cluster = 4.5 GeV
		eventCutArray[ 2] = "00052113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111063092230000"; mesonCutArray[2] = "0163103100000010"; // min Energy cluster = 6.0 GeV
	} else if (trainConfig == 115){ // EMCAL clusters, EMCEGA with minEnergy variation
		eventCutArray[ 0] = "00081113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // min Energy cluster = 0.4 GeV
		eventCutArray[ 1] = "00081113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111063062230000"; mesonCutArray[1] = "0163103100000010"; // min Energy cluster = 4.5 GeV
		eventCutArray[ 2] = "00081113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111063092230000"; mesonCutArray[2] = "0163103100000010"; // min Energy cluster = 6.0 GeV
	} else if (trainConfig == 116){ // EMCAL clusters, EMCEJE with minEnergy variation
		eventCutArray[ 0] = "00091113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // min Energy cluster = 0.4 GeV
		eventCutArray[ 1] = "00091113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111063062230000"; mesonCutArray[1] = "0163103100000010"; // min Energy cluster = 4.5 GeV
		eventCutArray[ 2] = "00091113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111063092230000"; mesonCutArray[2] = "0163103100000010"; // min Energy cluster = 6.0 GeV
	} else if (trainConfig == 117){  // EMCAL clusters, EMC triggers (EMC7, EMCEGA, EMCEJE)
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // INT7
		eventCutArray[ 1] = "00052113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000010"; // EMC7
		eventCutArray[ 2] = "00081113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000010"; // EMCEGA
		eventCutArray[ 3] = "00091113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111111063032230000"; mesonCutArray[3] = "0163103100000010"; // EMCEJE
	} else if (trainConfig == 118){ // EMCAL clusters, timing variation
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111053032230000"; mesonCutArray[0] = "0163103100000010"; // time 50ns
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111043032230000"; mesonCutArray[1] = "0163103100000010"; // time 100ns
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111033032230000"; mesonCutArray[2] = "0163103100000010"; // time 200ns
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111111023032230000"; mesonCutArray[3] = "0163103100000010"; // time 500ns
	} else if (trainConfig == 119){ // EMCAL clusters, timing variation
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // time
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111111073032230000"; mesonCutArray[1] = "0163103100000010"; // time
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111111083032230000"; mesonCutArray[2] = "0163103100000010"; // time
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "1111111093032230000"; mesonCutArray[3] = "0163103100000010"; // time
	} else if (trainConfig == 120){ // EMCAL clusters, MB for extendedQA
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // time
	} else if (trainConfig == 121){ // EMCAL clusters kEMC for extQA
		eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; //
	} else if (trainConfig == 122){ // EMCAL clusters EMCEGA for extQA
		eventCutArray[ 0] = "00081113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; //
	} else if (trainConfig == 123){ // EMCAL clusters EMCEJE for extQA
		eventCutArray[ 0] = "00091113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; //




	} else if (trainConfig == 126){ // PCM variations pi dEdx
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009317300008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327300008250400000"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009325000008250400000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 1: 0.3 ->
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009320000008250400000"; clusterCutArray[3] = "1111111063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 1: 0.5 ->
	} else if (trainConfig == 127){ // PCM variations pi dEdx
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327600008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327400008250400000"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009315600008250400000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009367400008250400000"; clusterCutArray[3] = "1111111063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
		eventCutArray[ 4] = "00000113"; photonCutArray[ 4] = "00200009347400008250400000"; clusterCutArray[4] = "1111111063032230000"; mesonCutArray[4] = "0163103100000010"; // dEdx pi: 3: 0.4-3, 1: 3. ->
	} else if (trainConfig == 128){ // PCM variations to close V0s
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111111063032230000"; mesonCutArray[0] = "0163103100000010"; // dEdx pi: 1: 0.4-2, -10: 2. ->
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250401000"; clusterCutArray[1] = "1111111063032230000"; mesonCutArray[1] = "0163103100000010"; // dEdx pi: 1: 0.4-3, -10: 3. ->
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250402000"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000010"; // dEdx pi: 0: 0.3-2, -10: 2. ->
		eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009327000008250403000"; clusterCutArray[3] = "1111111063032230000"; mesonCutArray[3] = "0163103100000010"; // dEdx pi: 2: 0.4-3, 0.5: 3. ->
	} else if (trainConfig == 129){ // EMCAL clusters No NonLinearity + Calo/ConvCalo kSDM
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "1111100063032230000"; mesonCutArray[0] = "0163103100000010"; // NonLinearity none
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "1111115063032230000"; mesonCutArray[1] = "0163103100000010"; // NonLinearity ConvCalo kSDM
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "1111116063032230000"; mesonCutArray[2] = "0163103100000010"; // NonLinearity Calo kSDM

	// ************************************* PHOS cuts ****************************************************
	// LHC12
	} else if (trainConfig == 131){ // PHOS clusters 8 TeV LHC12
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "2444400048033200000"; mesonCutArray[0] = "0163103100000010"; // 400 MeV cluster min energy
	} else if (trainConfig == 142){ // With/without Added Signals
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "2444400048033200000"; mesonCutArray[0] = "0163103100000010"; //
		eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[1] = "0163103100000010"; //

	
	// ************************************* PHOS cuts ****************************************************
	// LHC11a	
	} else if (trainConfig == 301) { //PHOS clusters
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "2444400047033200000"; mesonCutArray[0] = "0163103100000010";
		eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[1] = "0163103100000010";
		eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "2444400049033200000"; mesonCutArray[2] = "0163103100000010";
		eventCutArray[ 3] = "00061113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "2444400047033200000"; mesonCutArray[3] = "0163103100000010";
		eventCutArray[ 4] = "00061113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "2444400048033200000"; mesonCutArray[4] = "0163103100000010";
		eventCutArray[ 5] = "00061113"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "2444400049033200000"; mesonCutArray[5] = "0163103100000010";
	// LHC13g & LHC12x
	} else if (trainConfig == 302) { //PHOS clusters
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "2444400047033200000"; mesonCutArray[0] = "0163103100000010";
		eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[1] = "0163103100000010";
		eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "2444400049033200000"; mesonCutArray[2] = "0163103100000010";
		eventCutArray[ 3] = "00062113"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "2444400047033200000"; mesonCutArray[3] = "0163103100000010";
		eventCutArray[ 4] = "00062113"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "2444400048033200000"; mesonCutArray[4] = "0163103100000010";
		eventCutArray[ 5] = "00062113"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "2444400049033200000"; mesonCutArray[5] = "0163103100000010";
	// LHC11a
	} else if (trainConfig == 303) { //PHOS clusters without and with added signals
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "2444400047033200000"; mesonCutArray[0] = "0163103100000010";
		eventCutArray[ 1] = "00003123"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "2444400047033200000"; mesonCutArray[1] = "0163103100000010";
	
	} else {
		Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
	TList *ConvCutList = new TList();
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
	ConvCutList->SetOwner(kTRUE);
	AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
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
		
		analysisCuts[i] = new AliConversionPhotonCuts();
		analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());
		analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
		ConvCutList->Add(analysisCuts[i]);
		analysisCuts[i]->SetFillCutHistograms("",kFALSE);
	
		analysisClusterCuts[i] = new AliCaloPhotonCuts();
		analysisClusterCuts[i]->InitializeCutsFromCutString(clusterCutArray[i].Data());
		ClusterCutList->Add(analysisClusterCuts[i]);
		analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
		analysisClusterCuts[i]->SetFillCutHistograms("");
		
		analysisMesonCuts[i] = new AliConversionMesonCuts();
		analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
		MesonCutList->Add(analysisMesonCuts[i]);
		analysisMesonCuts[i]->SetFillCutHistograms("");
		analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
	}

	task->SetEventCutList(numberOfCuts,EventCutList);
	task->SetConversionCutList(numberOfCuts,ConvCutList);
	task->SetCaloCutList(numberOfCuts,ClusterCutList);
	task->SetMesonCutList(numberOfCuts,MesonCutList);
	task->SetMoveParticleAccordingToVertex(kTRUE);
	task->SetDoMesonAnalysis(kTRUE);
	task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
	task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
	task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
	task->SetUseTHnSparse(isUsingTHnSparse);
	//if(trainConfig==29) task->SetDoTreeConvGammaShowerShape(kTRUE);
	if(enableExtMatchAndQA == 2 || enableExtMatchAndQA == 3){ task->SetPlotHistsExtQA(kTRUE);}

	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer(Form("GammaConvCalo_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaConvCalo_%i.root",trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
