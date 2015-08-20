void AddTask_GammaCalo_pp(  Int_t 		trainConfig 				= 1,  								// change different set of cuts
							Int_t 		isMC   						= 0, 							// run MC
							Int_t 		enableQAMesonTask 			= 0, 								// enable QA in AliAnalysisTaskGammaConvV1
							Int_t 		enableQAClusterTask 		= 0, 								// enable additional QA task
							TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
                            TString 	cutnumberAODBranch 			= "000000006008400001001500000",
							TString 	periodname 					= "LHC12f1x", 						// period name
							Bool_t 		doWeighting 				= kFALSE,							// enables weighting
							Bool_t 		isUsingTHnSparse 			= kTRUE,							// enable or disable usage of THnSparses for background estimation
							Int_t 		enableExtQA					= 0,								// enable QA(3), disabled (0)
							Bool_t 		enableTriggerMimicking		= kFALSE							// enable trigger mimicking
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
	TString cutnumberEvent = "0000000";
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
	if (trainConfig == 2 || trainConfig == 3) numberOfCuts = 5;
	if (trainConfig == 5) numberOfCuts = 6;
	if (trainConfig == 31 || trainConfig == 51 || trainConfig == 53  || trainConfig == 55   || trainConfig == 57 ) numberOfCuts = 4;
	if (trainConfig == 32 || trainConfig == 101) numberOfCuts = 1;
	if (trainConfig == 59 || trainConfig == 60 || trainConfig == 61 || trainConfig == 62 || trainConfig == 112 || trainConfig == 113 || trainConfig == 114) numberOfCuts = 7;
	if (trainConfig == 111 || trainConfig == 121 || trainConfig == 52 || trainConfig == 54 || trainConfig == 56 || trainConfig == 58) numberOfCuts = 3;

	TString *eventCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// cluster cuts
	// 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
	// 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
	
	// ************************************* EMCAL cuts ****************************************************
	// LHC11a
	if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "11111050032230000"; mesonCutArray[0] = "0163103100000000"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "0005111"; clusterCutArray[1] = "11111050032230000"; mesonCutArray[1] = "0163103100000000"; // 400 MeV cluster min energy
	} else if (trainConfig == 2){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "11111050012230000"; mesonCutArray[0] = "0163103100000000"; //0.2 GeV/c
		eventCutArray[ 1] = "0000311"; clusterCutArray[1] = "11111050022230000"; mesonCutArray[1] = "0163103100000000"; //0.3 GeV/c
		eventCutArray[ 2] = "0000311"; clusterCutArray[2] = "11111050032230000"; mesonCutArray[2] = "0163103100000000"; //0.4 GeV/c default
		eventCutArray[ 3] = "0000311"; clusterCutArray[3] = "11111050042230000"; mesonCutArray[3] = "0163103100000000"; //0.5 GeV/c
		eventCutArray[ 4] = "0000311"; clusterCutArray[4] = "11111050052230000"; mesonCutArray[4] = "0163103100000000"; //0.6 GeV/c
	} else if (trainConfig == 3){ //EMCAL minNCells variation
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "11111050031230000"; mesonCutArray[0] = "0163103100000000"; //n cells >= 1
		eventCutArray[ 1] = "0000311"; clusterCutArray[1] = "11111050033230000"; mesonCutArray[1] = "0163103100000000"; //n cells >= 3
		eventCutArray[ 2] = "0000311"; clusterCutArray[2] = "11111050032000000"; mesonCutArray[2] = "0163103100000000"; //no M02 cut
		eventCutArray[ 3] = "0000311"; clusterCutArray[3] = "11131050032230000"; mesonCutArray[3] = "0163103100000000"; //only modules with TRD infront
		eventCutArray[ 4] = "0000311"; clusterCutArray[4] = "11112050032230000"; mesonCutArray[4] = "0163103100000000"; //no modules with TRD infront
	// LHC13g	
	} else if (trainConfig == 5){  // EMCAL clusters, EMCEGA triggers
		eventCutArray[ 0] = "0008311"; clusterCutArray[0] = "11111050032220000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1,
		eventCutArray[ 1] = "0008511"; clusterCutArray[1] = "11111050032220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG2,
		eventCutArray[ 2] = "0009311"; clusterCutArray[2] = "11111050032220000"; mesonCutArray[2] = "0163103100000000"; // EMCEJ1,
		eventCutArray[ 3] = "0009511"; clusterCutArray[3] = "11111050032220000"; mesonCutArray[3] = "0163103100000000"; // EMCEJ2,
		eventCutArray[ 4] = "0000011"; clusterCutArray[4] = "11111050032220000"; mesonCutArray[4] = "0163103100000000"; // INT7
		eventCutArray[ 5] = "0005211"; clusterCutArray[5] = "11111050032220000"; mesonCutArray[5] = "0163103100000000"; // EMC7
	} else if (trainConfig == 12){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0) without and with added signals
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "11111050032220000"; mesonCutArray[0] = "0163103100000000"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "0000312"; clusterCutArray[1] = "11111050032220000"; mesonCutArray[1] = "0163103100000000"; // 400 MeV cluster min energy
	} else if (trainConfig == 13){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
		eventCutArray[ 0] = "0000312"; clusterCutArray[0] = "11111050032230000"; mesonCutArray[0] = "0163103100000000"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "0005112"; clusterCutArray[1] = "11111050032230000"; mesonCutArray[1] = "0163103100000000"; // 400 MeV cluster min energy
	// ************************************* PHOS cuts ****************************************************
	} else if (trainConfig == 31) { //PHOS clusters
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "24444040033200000"; mesonCutArray[0] = "0163103100000000"; //pp LHC11a with SDD, PHOS
		eventCutArray[ 1] = "0000011"; clusterCutArray[1] = "24444040033200000"; mesonCutArray[1] = "0163103100000000"; //pp LHC13g default MB
		eventCutArray[ 2] = "0006111"; clusterCutArray[2] = "24444040033200000"; mesonCutArray[2] = "0163103100000000"; //pp LHC11a PHI1
		eventCutArray[ 3] = "0006211"; clusterCutArray[3] = "24444040033200000"; mesonCutArray[3] = "0163103100000000"; //pp LHC11a PHI7
	} else if (trainConfig == 32){ // Validation PHOS
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "24444040033200000"; mesonCutArray[0] = "0163003100900000";
	} else if (trainConfig == 33){ // PHOS clusters, without and with added signals
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "24444040033200000"; mesonCutArray[0] = "0163003100900000";
		eventCutArray[ 1] = "0000312"; clusterCutArray[1] = "24444040033200000"; mesonCutArray[1] = "0163003100900000";
	  // LHC13g cut studies
	} else if (trainConfig == 51){  // EMCAL clusters, EMCEG1 trigger
		eventCutArray[ 0] = "0008311"; clusterCutArray[0] = "11111050032220000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1, 400 MeV min energy, NCells >=2, M02 default cut, 50ns timing
		eventCutArray[ 1] = "0008311"; clusterCutArray[1] = "11111050052220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG1, 600 MeV min energy
		eventCutArray[ 2] = "0008311"; clusterCutArray[2] = "11111050031220000"; mesonCutArray[2] = "0163103100000000"; // EMCEG1,                     NCells >=1
		eventCutArray[ 3] = "0008311"; clusterCutArray[3] = "11111050033220000"; mesonCutArray[3] = "0163103100000000"; // EMCEG1,                     NCells >=3
	} else if (trainConfig == 52){  // EMCAL clusters, EMCEG1 trigger
		eventCutArray[ 0] = "0008311"; clusterCutArray[0] = "11111050032000000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1,                                 no M02 cut
		eventCutArray[ 1] = "0008311"; clusterCutArray[1] = "11111020032220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG1, 						 	500ns timing
		eventCutArray[ 2] = "0008311"; clusterCutArray[2] = "11111040032220000"; mesonCutArray[2] = "0163103100000000"; // EMCEG1,							100ns timing
	} else if (trainConfig == 53){  // EMCAL clusters, EMCEG2 trigger
		eventCutArray[ 0] = "0008511"; clusterCutArray[0] = "11111050032220000"; mesonCutArray[0] = "0163103100000000"; // EMCEG2, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "0008511"; clusterCutArray[1] = "11111050052220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG2, 600 MeV min energy
		eventCutArray[ 2] = "0008511"; clusterCutArray[2] = "11111050031220000"; mesonCutArray[2] = "0163103100000000"; // EMCEG2,                     NCells >=1
		eventCutArray[ 3] = "0008511"; clusterCutArray[3] = "11111050033220000"; mesonCutArray[3] = "0163103100000000"; // EMCEG2,                     NCells >=3
	} else if (trainConfig == 54){  // EMCAL clusters, EMCEG2 trigger	
		eventCutArray[ 0] = "0008511"; clusterCutArray[0] = "11111050032000000"; mesonCutArray[0] = "0163103100000000"; // EMCEG2,                                 no M02 cut
		eventCutArray[ 1] = "0008511"; clusterCutArray[1] = "11111020032220000"; mesonCutArray[1] = "0163103100000000"; // EMCEG2,                                                 500ns timing
		eventCutArray[ 2] = "0008511"; clusterCutArray[2] = "11111040032220000"; mesonCutArray[2] = "0163103100000000"; // EMCEG2,                                                 100ns timing
	} else if (trainConfig == 55){  // EMCAL clusters, INT7 trigger
		eventCutArray[ 0] = "0000011"; clusterCutArray[0] = "11111050032220000"; mesonCutArray[0] = "0163103100000000"; // INT7, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "0000011"; clusterCutArray[1] = "11111050052220000"; mesonCutArray[1] = "0163103100000000"; // INT7, 600 MeV min energy
		eventCutArray[ 2] = "0000011"; clusterCutArray[2] = "11111050031220000"; mesonCutArray[2] = "0163103100000000"; // INT7,                       NCells >=1
		eventCutArray[ 3] = "0000011"; clusterCutArray[3] = "11111050033220000"; mesonCutArray[3] = "0163103100000000"; // INT7,                       NCells >=3
	} else if (trainConfig == 56){  // EMCAL clusters, INT7 trigger
		eventCutArray[ 0] = "0000011"; clusterCutArray[0] = "11111050032000000"; mesonCutArray[0] = "0163103100000000"; // INT7,                                   no M02 cut
		eventCutArray[ 1] = "0000011"; clusterCutArray[1] = "11111020032220000"; mesonCutArray[1] = "0163103100000000"; // INT7,                                                   500ns timing
		eventCutArray[ 2] = "0000011"; clusterCutArray[2] = "11111040032220000"; mesonCutArray[2] = "0163103100000000"; // INT7,                                                   100ns timing
	} else if (trainConfig == 57){  // EMCAL clusters, EMC7 trigger
		eventCutArray[ 0] = "0005211"; clusterCutArray[0] = "11111050032220000"; mesonCutArray[0] = "0163103100000000"; // EMC7, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "0005211"; clusterCutArray[1] = "11111050052220000"; mesonCutArray[1] = "0163103100000000"; // EMC7, 600 MeV min energy
		eventCutArray[ 2] = "0005211"; clusterCutArray[2] = "11111050031220000"; mesonCutArray[2] = "0163103100000000"; // EMC7,                     NCells >=1
		eventCutArray[ 3] = "0005211"; clusterCutArray[3] = "11111050033220000"; mesonCutArray[3] = "0163103100000000"; // EMC7,                     NCells >=3
	} else if (trainConfig == 58){  // EMCAL clusters, EMC7 trigger
		eventCutArray[ 0] = "0005211"; clusterCutArray[0] = "11111050032000000"; mesonCutArray[0] = "0163103100000000"; // EMC7,                                 no M02 cut
		eventCutArray[ 1] = "0005211"; clusterCutArray[1] = "11111020032220000"; mesonCutArray[1] = "0163103100000000"; // EMC7,                                                 500ns timing
		eventCutArray[ 2] = "0005211"; clusterCutArray[2] = "11111040032220000"; mesonCutArray[2] = "0163103100000000"; // EMC7,                                                 100ns timing
	  // LHC11a cut studies
	} else if (trainConfig == 59){  // EMCAL clusters, MB (INT1) trigger
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "11111050032220000"; mesonCutArray[0] = "0163103100000000"; // MB, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "0000311"; clusterCutArray[1] = "11111050052220000"; mesonCutArray[1] = "0163103100000000"; // MB, 600 MeV min energy
		eventCutArray[ 2] = "0000311"; clusterCutArray[2] = "11111050031220000"; mesonCutArray[2] = "0163103100000000"; // MB,                     NCells >=1
		eventCutArray[ 3] = "0000311"; clusterCutArray[3] = "11111050033220000"; mesonCutArray[3] = "0163103100000000"; // MB,                     NCells >=3
		eventCutArray[ 4] = "0000311"; clusterCutArray[4] = "11111050032000000"; mesonCutArray[4] = "0163103100000000"; // MB,                                 no M02 cut
		eventCutArray[ 5] = "0000311"; clusterCutArray[5] = "11111020032220000"; mesonCutArray[5] = "0163103100000000"; // MB,                                                 500ns timing
		eventCutArray[ 6] = "0000311"; clusterCutArray[6] = "11111040032220000"; mesonCutArray[6] = "0163103100000000"; // MB,                                                 100ns timing
	} else if (trainConfig == 60){  // EMCAL clusters, EMC1 trigger
		eventCutArray[ 0] = "0005111"; clusterCutArray[0] = "11111050032220000"; mesonCutArray[0] = "0163103100000000"; // EMC1, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "0005111"; clusterCutArray[1] = "11111050052220000"; mesonCutArray[1] = "0163103100000000"; // EMC1, 600 MeV min energy
		eventCutArray[ 2] = "0005111"; clusterCutArray[2] = "11111050031220000"; mesonCutArray[2] = "0163103100000000"; // EMC1,                     NCells >=1
		eventCutArray[ 3] = "0005111"; clusterCutArray[3] = "11111050033220000"; mesonCutArray[3] = "0163103100000000"; // EMC1,                     NCells >=3
		eventCutArray[ 4] = "0005111"; clusterCutArray[4] = "11111050032000000"; mesonCutArray[4] = "0163103100000000"; // EMC1,                                 no M02 cut
		eventCutArray[ 5] = "0005111"; clusterCutArray[5] = "11111020032220000"; mesonCutArray[5] = "0163103100000000"; // EMC1,                                                 500ns timing
		eventCutArray[ 6] = "0005111"; clusterCutArray[6] = "11111040032220000"; mesonCutArray[6] = "0163103100000000"; // EMC1,                                                 100ns timing
	} else if (trainConfig == 61){  // EMCAL clusters, MB (INT1) trigger, for added signals
		eventCutArray[ 0] = "0000312"; clusterCutArray[0] = "11111050032220000"; mesonCutArray[0] = "0163103100000000"; // MB, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "0000312"; clusterCutArray[1] = "11111050052220000"; mesonCutArray[1] = "0163103100000000"; // MB, 600 MeV min energy
		eventCutArray[ 2] = "0000312"; clusterCutArray[2] = "11111050031220000"; mesonCutArray[2] = "0163103100000000"; // MB,                     NCells >=1
		eventCutArray[ 3] = "0000312"; clusterCutArray[3] = "11111050033220000"; mesonCutArray[3] = "0163103100000000"; // MB,                     NCells >=3
		eventCutArray[ 4] = "0000312"; clusterCutArray[4] = "11111050032000000"; mesonCutArray[4] = "0163103100000000"; // MB,                                 no M02 cut
		eventCutArray[ 5] = "0000312"; clusterCutArray[5] = "11111020032220000"; mesonCutArray[5] = "0163103100000000"; // MB,                                                 500ns timing
		eventCutArray[ 6] = "0000312"; clusterCutArray[6] = "11111040032220000"; mesonCutArray[6] = "0163103100000000"; // MB,                                                 100ns timing
	} else if (trainConfig == 62){  // EMCAL clusters, EMC1 trigger
		eventCutArray[ 0] = "0005112"; clusterCutArray[0] = "11111050032220000"; mesonCutArray[0] = "0163103100000000"; // EMC1, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "0005112"; clusterCutArray[1] = "11111050052220000"; mesonCutArray[1] = "0163103100000000"; // EMC1, 600 MeV min energy
		eventCutArray[ 2] = "0005112"; clusterCutArray[2] = "11111050031220000"; mesonCutArray[2] = "0163103100000000"; // EMC1,                     NCells >=1
		eventCutArray[ 3] = "0005112"; clusterCutArray[3] = "11111050033220000"; mesonCutArray[3] = "0163103100000000"; // EMC1,                     NCells >=3
		eventCutArray[ 4] = "0005112"; clusterCutArray[4] = "11111050032000000"; mesonCutArray[4] = "0163103100000000"; // EMC1,                                 no M02 cut
		eventCutArray[ 5] = "0005112"; clusterCutArray[5] = "11111020032220000"; mesonCutArray[5] = "0163103100000000"; // EMC1,                                                 500ns timing
		eventCutArray[ 6] = "0005112"; clusterCutArray[6] = "11111040032220000"; mesonCutArray[6] = "0163103100000000"; // EMC1,                                                 100ns timing
	} else if (trainConfig == 101){ // EMCAL clusters pp 8 TeV / 7 TeV
		eventCutArray[ 0] = "0000011"; clusterCutArray[0] = "11111050032230000"; mesonCutArray[0] = "0163103100000000"; // 400 MeV cluster min energy
   // LHC12fa-i and MC
	} else if (trainConfig == 111){  // EMCAL clusters, EMCEGA triggers
		eventCutArray[ 0] = "0008111"; clusterCutArray[0] = "11111050032230000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1,
		eventCutArray[ 1] = "0000011"; clusterCutArray[1] = "11111050032230000"; mesonCutArray[1] = "0163103100000000"; // INT7
		eventCutArray[ 2] = "0005211"; clusterCutArray[2] = "11111050032230000"; mesonCutArray[2] = "0163103100000000"; // EMC7
	} else if (trainConfig == 112){  // EMCAL clusters, EMCEG2 trigger
		eventCutArray[ 0] = "0008111"; clusterCutArray[0] = "11111050032230000"; mesonCutArray[0] = "0163103100000000"; // EMCEGA, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "0008111"; clusterCutArray[1] = "11111050052230000"; mesonCutArray[1] = "0163103100000000"; // EMCEGA, 600 MeV min energy
		eventCutArray[ 2] = "0008111"; clusterCutArray[2] = "11111050031230000"; mesonCutArray[2] = "0163103100000000"; // EMCEGA,                     NCells >=1
		eventCutArray[ 3] = "0008111"; clusterCutArray[3] = "11111050033230000"; mesonCutArray[3] = "0163103100000000"; // EMCEGA,                     NCells >=3
		eventCutArray[ 4] = "0008111"; clusterCutArray[4] = "11111050032000000"; mesonCutArray[4] = "0163103100000000"; // EMCEGA,                                 no M02 cut
		eventCutArray[ 5] = "0008111"; clusterCutArray[5] = "11111020032230000"; mesonCutArray[5] = "0163103100000000"; // EMCEGA,                                                 500ns timing
		eventCutArray[ 6] = "0008511"; clusterCutArray[6] = "11111040032230000"; mesonCutArray[6] = "0163103100000000"; // EMCEGA,                                                 100ns timing
	} else if (trainConfig == 113){  // EMCAL clusters, INT7 trigger
		eventCutArray[ 0] = "0000011"; clusterCutArray[0] = "11111050032230000"; mesonCutArray[0] = "0163103100000000"; // INT7, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "0000011"; clusterCutArray[1] = "11111050052230000"; mesonCutArray[1] = "0163103100000000"; // INT7, 600 MeV min energy
		eventCutArray[ 2] = "0000011"; clusterCutArray[2] = "11111050031230000"; mesonCutArray[2] = "0163103100000000"; // INT7,                       NCells >=1
		eventCutArray[ 3] = "0000011"; clusterCutArray[3] = "11111050033230000"; mesonCutArray[3] = "0163103100000000"; // INT7,                       NCells >=3
		eventCutArray[ 4] = "0000011"; clusterCutArray[4] = "11111050032000000"; mesonCutArray[4] = "0163103100000000"; // INT7,                                   no M02 cut
		eventCutArray[ 5] = "0000011"; clusterCutArray[5] = "11111020032230000"; mesonCutArray[5] = "0163103100000000"; // INT7,                                                   500ns timing
		eventCutArray[ 6] = "0000011"; clusterCutArray[6] = "11111040032230000"; mesonCutArray[6] = "0163103100000000"; // INT7,                                                   100ns timing
	} else if (trainConfig == 114){  // EMCAL clusters, EMC7 trigger
		eventCutArray[ 0] = "0005211"; clusterCutArray[0] = "11111050032230000"; mesonCutArray[0] = "0163103100000000"; // EMC7, 400 MeV min energy, NCells >=2, M02 default cut
		eventCutArray[ 1] = "0005211"; clusterCutArray[1] = "11111050052230000"; mesonCutArray[1] = "0163103100000000"; // EMC7, 600 MeV min energy
		eventCutArray[ 2] = "0005211"; clusterCutArray[2] = "11111050031230000"; mesonCutArray[2] = "0163103100000000"; // EMC7,                     NCells >=1
		eventCutArray[ 3] = "0005211"; clusterCutArray[3] = "11111050033230000"; mesonCutArray[3] = "0163103100000000"; // EMC7,                     NCells >=3
		eventCutArray[ 4] = "0005211"; clusterCutArray[4] = "11111050032000000"; mesonCutArray[4] = "0163103100000000"; // EMC7,                                 no M02 cut
		eventCutArray[ 5] = "0005211"; clusterCutArray[5] = "11111020032230000"; mesonCutArray[5] = "0163103100000000"; // EMC7,                                                 500ns timing
		eventCutArray[ 6] = "0005211"; clusterCutArray[6] = "11111040032230000"; mesonCutArray[6] = "0163103100000000"; // EMC7,                                                 100ns timing
	} else if (trainConfig == 121){ // EMCAL clusters, EMCEGA triggers added signals
		eventCutArray[ 0] = "0008112"; clusterCutArray[0] = "11111050032230000"; mesonCutArray[0] = "0163103100000000"; // EMCEG1,
		eventCutArray[ 1] = "0000012"; clusterCutArray[1] = "11111050032230000"; mesonCutArray[1] = "0163103100000000"; // INT7
		eventCutArray[ 2] = "0005212"; clusterCutArray[2] = "11111050032230000"; mesonCutArray[2] = "0163103100000000"; // EMC7
	}
	else {
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
		Bool_t fAddedSignal = eventCutArray[i].EndsWith("2");
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

		analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
		EventCutList->Add(analysisEventCuts[i]);
		analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
		analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
		
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
