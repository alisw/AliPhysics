void AddTask_GammaConvCalo_pp(  	Int_t 		trainConfig 				= 1,  								//change different set of cuts
									Bool_t	 	isMC 						= kFALSE, 							//run MC
									Int_t 		enableQAMesonTask 			= 1, 								//enable QA in AliAnalysisTaskGammaConvV1
									Int_t 		enableQAPhotonTask 			= 1, 								// enable additional QA task
									TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
									TString 	cutnumberAODBranch 			= "000000006008400001001500000",
									Bool_t 		enableExtendedMatching 		= kFALSE, 							//enable or disable extended matching histograms for conversion electrons <-> cluster
									TString 	periodname 					= "LHC12f1x", 						// period name
									Bool_t 		doWeighting 				= kFALSE,							// enables weighting
									Bool_t 		enableV0findingEffi 		= kFALSE							// enables V0finding efficiency histograms
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
		Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
		return ;
	}

	// ================== GetInputEventHandler =============================
	AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
	
	//========= Add PID Reponse to ANALYSIS manager ====
	if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
		gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
		AddTaskPIDResponse(isMC);
	}
	
	Printf("here \n");
	
	//=========  Set Cutnumber for V0Reader ================================
	TString cutnumberPhoton = "00000000600840001001500000";
	TString cutnumberEvent = "0000000";
	Bool_t doEtaShift = kFALSE;
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

	//========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
	if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
		AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
		
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
    if (trainConfig==8 || trainConfig==10) {numberOfCuts = 4;}
    if (trainConfig==2 || trainConfig==3 || trainConfig==5 || trainConfig==6 || trainConfig==7 ) {numberOfCuts = 5;}
    if (trainConfig==4 || trainConfig==11 || trainConfig==31 || trainConfig==32 ) {numberOfCuts = 6;}

	TString *eventCutArray = new TString[numberOfCuts];
	TString *photonCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// cluster cuts
	// 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
	// 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"

	// ************************************* EMCAL cuts ****************************************************
	// LHC11a
	if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
        eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "10000043032230000"; mesonCutArray[0] = "01631031000000"; // 400 MeV cluster min energy
        eventCutArray[ 1] = "0005111"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "10000043032230000"; mesonCutArray[1] = "01631031000000"; // 400 MeV cluster min energy
	} else if (trainConfig == 2){ //EMCAL minEnergy variation
        eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "10000043012230000"; mesonCutArray[0] = "01631031000000"; //0.2 GeV/c
        eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "10000043022230000"; mesonCutArray[1] = "01631031000000"; //0.3 GeV/c
        eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "10000043032230000"; mesonCutArray[2] = "01631031000000"; //0.4 GeV/c default
        eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "10000043042230000"; mesonCutArray[3] = "01631031000000"; //0.5 GeV/c
        eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "10000043052230000"; mesonCutArray[4] = "01631031000000"; //0.6 GeV/c
	} else if (trainConfig == 3){ //EMCAL minNCells variation
        eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "10000043031230000"; mesonCutArray[0] = "01631031000000"; //n cells >= 1
        eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "10000043033230000"; mesonCutArray[1] = "01631031000000"; //n cells >= 3
        eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "10000043032200000"; mesonCutArray[2] = "01631031000000"; //no M02 cut
        eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "10031043032230000"; mesonCutArray[3] = "01631031000000"; //only modules with TRD infront
        eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "10012043032230000"; mesonCutArray[4] = "01631031000000"; //no modules with TRD infront
	} else if (trainConfig == 4){ // EMCAL track matching variations 
        eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "10000041032230000"; mesonCutArray[0] = "01631031000000"; //
        eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "10000042032230000"; mesonCutArray[1] = "01631031000000"; //
        eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "10000043032230000"; mesonCutArray[2] = "01631031000000"; //
        eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "10000044032230000"; mesonCutArray[3] = "01631031000000"; //
        eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "10000045032230000"; mesonCutArray[4] = "01631031000000"; //
        eventCutArray[ 5] = "0000311"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "10000046032230000"; mesonCutArray[5] = "01631031000000"; //
	} else if (trainConfig == 5){ // PCM variations
        eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00200009227000008250400000"; clusterCutArray[0] = "10000043032230000"; mesonCutArray[0] = "01631031000000"; // dEdx e -3, 5
        eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "00200009127000008250400000"; clusterCutArray[1] = "10000043032230000"; mesonCutArray[1] = "01631031000000"; // dEdx e -5, 5
        eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "00200009357000008250400000"; clusterCutArray[2] = "10000043032230000"; mesonCutArray[2] = "01631031000000"; // dEdx pi 2
        eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "00200009317000008250400000"; clusterCutArray[3] = "10000043032230000"; mesonCutArray[3] = "01631031000000"; // dEdx pi 0
        eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "00200009387300008250400000"; clusterCutArray[4] = "10000043032230000"; mesonCutArray[4] = "01631031000000"; // dEdx pi 2 high 1 (> 3.5 GeV)
	} else if (trainConfig == 6){ // PCM variations
        eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00200009327000009250400000"; clusterCutArray[0] = "10000043032230000"; mesonCutArray[0] = "01631031000000"; // qt 2D 0.03
        eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "00200009327000003250400000"; clusterCutArray[1] = "10000043032230000"; mesonCutArray[1] = "01631031000000"; // qt 1D 0.05
        eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "00200009327000002250400000"; clusterCutArray[2] = "10000043032230000"; mesonCutArray[2] = "01631031000000"; // qt 1D 0.07
        eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "00200049327000008250400000"; clusterCutArray[3] = "10000043032230000"; mesonCutArray[3] = "01631031000000"; // single pt > 0.075
        eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "00200019327000008250400000"; clusterCutArray[4] = "10000043032230000"; mesonCutArray[4] = "01631031000000"; // single pt > 0.1
	} else if (trainConfig == 7){ // PCM variations
        eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00200009327000008850400000"; clusterCutArray[0] = "10000043032230000"; mesonCutArray[0] = "01631031000000"; // 2D psi pair chi2 var
        eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "00200009327000008260400000"; clusterCutArray[1] = "10000043032230000"; mesonCutArray[1] = "01631031000000"; // 2D psi pair chi2 var
        eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "00200009327000008860400000"; clusterCutArray[2] = "10000043032230000"; mesonCutArray[2] = "01631031000000"; // 2D psi pair chi2 var
        eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "00200009327000008280400000"; clusterCutArray[3] = "10000043032230000"; mesonCutArray[3] = "01631031000000"; // 2D psi pair chi2 var
        eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "00200009327000008880400000"; clusterCutArray[4] = "10000043032230000"; mesonCutArray[4] = "01631031000000"; // 2D psi pair chi2 var
	} else if (trainConfig == 8){ // PCM variations
        eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00200006327000008250400000"; clusterCutArray[0] = "10000043032230000"; mesonCutArray[0] = "01631031000000"; // min TPC cl > 0.7
        eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "00200008327000008250400000"; clusterCutArray[1] = "10000043032230000"; mesonCutArray[1] = "01631031000000"; // min TPC cl > 0.35
        eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "10000043032230000"; mesonCutArray[2] = "01631061000000"; // alpha < 0.8
        eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "10000043032230000"; mesonCutArray[3] = "01631051000000"; // alpha < 0.75
	} else if (trainConfig == 9){ // PCM variations
        eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00202209327000008250400000"; clusterCutArray[0] = "10000043032230000"; mesonCutArray[0] = "01631031000000"; // restrict acceptance to EMCAL loose
        eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "00204409327000008250400000"; clusterCutArray[1] = "10000043032230000"; mesonCutArray[1] = "01631031000000"; // restrict acceptance to EMCAL tight
	} else if (trainConfig == 10){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
        eventCutArray[ 0] = "0005111"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "10000043062230000"; mesonCutArray[0] = "01631031000000"; // min Energy cluster = 4.5 GeV
        eventCutArray[ 1] = "0005111"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "10000043072230000"; mesonCutArray[1] = "01631031000000"; // min Energy cluster = 5.0 GeV
        eventCutArray[ 2] = "0005111"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "10000043082230000"; mesonCutArray[2] = "01631031000000"; // min Energy cluster = 5.5 GeV
        eventCutArray[ 3] = "0005111"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "10000043092230000"; mesonCutArray[3] = "01631031000000"; // min Energy cluster = 6.0 GeV
		// LHC13g	
	} else if (trainConfig == 11){  // EMCAL clusters, EMCEGA triggers, track matching 0.035
        eventCutArray[ 0] = "0008311"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "10000043032230000"; mesonCutArray[0] = "01631031000000"; // EMCEG1,
        eventCutArray[ 1] = "0008511"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "10000043032230000"; mesonCutArray[1] = "01631031000000"; // EMCEG2,
        eventCutArray[ 2] = "0009311"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "10000043032230000"; mesonCutArray[2] = "01631031000000"; // EMCEJ1,
        eventCutArray[ 3] = "0009511"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "10000043032230000"; mesonCutArray[3] = "01631031000000"; // EMCEJ2,
        eventCutArray[ 4] = "0000011"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "10000043032230000"; mesonCutArray[4] = "01631031000000"; // INT7
        eventCutArray[ 5] = "0005211"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "10000043032230000"; mesonCutArray[5] = "01631031000000"; // EMC7
	// LHC11a	
	} else if (trainConfig == 12){ // EMCAL clusters 2.76 TeV LHC11a, with SDD without and with added signals
        eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "10000043032230000"; mesonCutArray[0] = "01631031000000"; // 400 MeV cluster min energy
        eventCutArray[ 1] = "0000312"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "10000043032230000"; mesonCutArray[1] = "01631031000000"; // 400 MeV cluster min energy
		
	// ************************************* PHOS cuts ****************************************************
	// LHC11a	
	} else if (trainConfig == 31) { //PHOS clusters
        eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "20000047033200000"; mesonCutArray[0] = "01631031000000";
        eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "20000048033200000"; mesonCutArray[1] = "01631031000000";
        eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "20000049033200000"; mesonCutArray[2] = "01631031000000";
        eventCutArray[ 3] = "0006111"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "20000047033200000"; mesonCutArray[3] = "01631031000000";
        eventCutArray[ 4] = "0006111"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "20000048033200000"; mesonCutArray[4] = "01631031000000";
        eventCutArray[ 5] = "0006111"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "20000049033200000"; mesonCutArray[5] = "01631031000000";
	// LHC13g & LHC12x
	} else if (trainConfig == 32) { //PHOS clusters
        eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "20000047033200000"; mesonCutArray[0] = "01631031000000";
        eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "20000048033200000"; mesonCutArray[1] = "01631031000000";
        eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "00200009327000008250400000"; clusterCutArray[2] = "20000049033200000"; mesonCutArray[2] = "01631031000000";
        eventCutArray[ 3] = "0006211"; photonCutArray[ 3] = "00200009327000008250400000"; clusterCutArray[3] = "20000047033200000"; mesonCutArray[3] = "01631031000000";
        eventCutArray[ 4] = "0006211"; photonCutArray[ 4] = "00200009327000008250400000"; clusterCutArray[4] = "20000048033200000"; mesonCutArray[4] = "01631031000000";
        eventCutArray[ 5] = "0006211"; photonCutArray[ 5] = "00200009327000008250400000"; clusterCutArray[5] = "20000049033200000"; mesonCutArray[5] = "01631031000000";
	// LHC11a
	} else if (trainConfig == 33) { //PHOS clusters without and with added signals
		eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "00200009327000008250400000"; clusterCutArray[0] = "20000047033200000"; mesonCutArray[0] = "01631031000000";
		eventCutArray[ 1] = "0000312"; photonCutArray[ 1] = "00200009327000008250400000"; clusterCutArray[1] = "20000047033200000"; mesonCutArray[1] = "01631031000000";
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
		
		analysisCuts[i] = new AliConversionPhotonCuts();
		analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());
		analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
		ConvCutList->Add(analysisCuts[i]);
		analysisCuts[i]->SetFillCutHistograms("",kFALSE);
	
		analysisClusterCuts[i] = new AliCaloPhotonCuts();
		analysisClusterCuts[i]->InitializeCutsFromCutString(clusterCutArray[i].Data());
		ClusterCutList->Add(analysisClusterCuts[i]);
        analysisClusterCuts[i]->SetExtendedMatching(enableExtendedMatching);
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

	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer(Form("GammaConvCalo_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaConvCalo_%i.root",trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
