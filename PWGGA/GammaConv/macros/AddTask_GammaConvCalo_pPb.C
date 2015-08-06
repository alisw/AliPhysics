void AddTask_GammaConvCalo_pPb(  	Int_t 		trainConfig 				= 1,  								// change different set of cuts
									Int_t 		isMC   						= 0, 							// run MC
									Int_t 		enableQAMesonTask 			= 0, 								// enable QA in AliAnalysisTaskGammaConvV1
									Int_t 		enableQAPhotonTask 			= 0, 								// enable additional QA task
									TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
									Int_t 		doWeightingPart 			= 0,  								// enable Weighting
									TString 	generatorName 				= "DPMJET",							// generator Name	
									TString 	cutnumberAODBranch 			= "800000006008400000001500000",	// cutnumber for AOD branch
									Int_t 		enableExtMatchAndQA 		= 0,								// enable matching histograms (1) and extended QA (2), only QA(3), all disabled (0)
									Bool_t 		isUsingTHnSparse 			= kTRUE, 							// enable or disable usage of THnSparses for background estimation
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

	Int_t isHeavyIon = 2;
	
	// ================== GetAnalysisManager ===============================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error(Form("AddTask_GammaConvCalo_pPb_%i",trainConfig), "No analysis manager found.");
		return ;
	}

	// ================== GetInputEventHandler =============================
	AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
	
	Bool_t isMCForOtherTasks = kFALSE;
	if (isMC > 0) isMCForOtherTasks = kTRUE;

	
	//========= Add PID Reponse to ANALYSIS manager ====
	if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
		gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
		AddTaskPIDResponse(isMCForOtherTasks);
	}
	
	Printf("here \n");
	
	//=========  Set Cutnumber for V0Reader ================================
	TString cutnumberPhoton = "00000008400100001500000000";
	TString cutnumberEvent = "80000003";
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

    if ( trainConfig==32 || trainConfig==33 ) {numberOfCuts = 3;}
    if ( trainConfig==10 ) {numberOfCuts = 4;}
    if ( trainConfig==7 || trainConfig==8 ) {numberOfCuts = 5;}
    if ( trainConfig==5 || trainConfig==6 ) {numberOfCuts = 6;}
	
	TString *eventCutArray = new TString[numberOfCuts];
	TString *photonCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// cluster cuts
	// 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
	// 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
	
	//************************************************ EMCAL clusters **********************************************************
	if (trainConfig == 1){ // min energy = 0.3 GeV/c
		eventCutArray[ 0] = "80000013"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "1111100053022230000"; mesonCutArray[0] = "0163103100000000"; //standart cut, kINT7
		eventCutArray[ 1] = "80052013"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "1111100053022230000"; mesonCutArray[1] = "0163103100000000"; //standard cut, kEMC7
	} else if (trainConfig == 2){  // min energy = 0.3 GeV/c
		eventCutArray[ 0] = "80083013"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "1111100053022230000"; mesonCutArray[0] = "0163103100000000"; //standard cut, kEMCEG1 based on INT7
		eventCutArray[ 1] = "80085013"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "1111100053022230000"; mesonCutArray[1] = "0163103100000000"; //standard cut, kEMCEG2 based on INT7
	} else if (trainConfig == 3){ // min energy = 0.4 GeV/c
		eventCutArray[ 0] = "80000013"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "1111100053032230000"; mesonCutArray[0] = "0163103100000000"; //standart cut, kINT7
		eventCutArray[ 1] = "80052013"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "1111100053032230000"; mesonCutArray[1] = "0163103100000000"; //standard cut, kEMC7
    } else if (trainConfig == 4){ // min energy = 0.4 GeV/
		eventCutArray[ 0] = "80083013"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "1111100053032230000"; mesonCutArray[0] = "0163103100000000"; //standard cut, kEMCEG1 based on INT7
		eventCutArray[ 1] = "80085013"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "1111100053032230000"; mesonCutArray[1] = "0163103100000000"; //standard cut, kEMCEG2 based on INT7
	} else if (trainConfig == 5){ //EMCAL variation of track matching
		eventCutArray[ 0] = "80000013"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "1111100051032230000"; mesonCutArray[0] = "0163103100000000"; //
		eventCutArray[ 1] = "80000013"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "1111100052032230000"; mesonCutArray[1] = "0163103100000000";
		eventCutArray[ 2] = "80000013"; photonCutArray[ 2] = "00200009327002008250400000"; clusterCutArray[2] = "1111100053032230000"; mesonCutArray[2] = "0163103100000000";
		eventCutArray[ 3] = "80000013"; photonCutArray[ 3] = "00200009327002008250400000"; clusterCutArray[3] = "1111100054032230000"; mesonCutArray[3] = "0163103100000000";
		eventCutArray[ 4] = "80000013"; photonCutArray[ 4] = "00200009327002008250400000"; clusterCutArray[4] = "1111100055032230000"; mesonCutArray[4] = "0163103100000000";
		eventCutArray[ 5] = "80000013"; photonCutArray[ 5] = "00200009327002008250400000"; clusterCutArray[5] = "1111100056032230000"; mesonCutArray[5] = "0163103100000000";
	} else if (trainConfig == 6){ //EMCAL added signal
		eventCutArray[ 0] = "80000023"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "1111100051032230000"; mesonCutArray[0] = "0163103100000000";
		eventCutArray[ 1] = "80000023"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "1111100052032230000"; mesonCutArray[1] = "0163103100000000";
		eventCutArray[ 2] = "80000023"; photonCutArray[ 2] = "00200009327002008250400000"; clusterCutArray[2] = "1111100053032230000"; mesonCutArray[2] = "0163103100000000";
		eventCutArray[ 3] = "80000023"; photonCutArray[ 3] = "00200009327002008250400000"; clusterCutArray[3] = "1111100054032230000"; mesonCutArray[3] = "0163103100000000";
		eventCutArray[ 4] = "80000023"; photonCutArray[ 4] = "00200009327002008250400000"; clusterCutArray[4] = "1111100055032230000"; mesonCutArray[4] = "0163103100000000";
		eventCutArray[ 5] = "80000023"; photonCutArray[ 5] = "00200009327002008250400000"; clusterCutArray[5] = "1111100056032230000"; mesonCutArray[5] = "0163103100000000";
	} else if (trainConfig == 7){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "80000013"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "1111100053012230000"; mesonCutArray[0] = "0163103100000000"; //0.2 GeV/c
		eventCutArray[ 1] = "80000013"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "1111100053022230000"; mesonCutArray[1] = "0163103100000000"; //0.3 GeV/c
		eventCutArray[ 2] = "80000013"; photonCutArray[ 2] = "00200009327002008250400000"; clusterCutArray[2] = "1111100053032230000"; mesonCutArray[2] = "0163103100000000"; //0.4 GeV/c default
		eventCutArray[ 3] = "80000013"; photonCutArray[ 3] = "00200009327002008250400000"; clusterCutArray[3] = "1111100053042230000"; mesonCutArray[3] = "0163103100000000"; //0.5 GeV/c
		eventCutArray[ 4] = "80000013"; photonCutArray[ 4] = "00200009327002008250400000"; clusterCutArray[4] = "1111100053052230000"; mesonCutArray[4] = "0163103100000000"; //0.6 GeV/c
	} else if (trainConfig == 8){ //EMCAL minNCells variation
		eventCutArray[ 0] = "80000013"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "1111100053031230000"; mesonCutArray[0] = "0163103100000000"; //n cells >= 1
		eventCutArray[ 1] = "80000013"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "1111100053033230000"; mesonCutArray[1] = "0163103100000000"; //n cells >= 3
		eventCutArray[ 2] = "80000013"; photonCutArray[ 2] = "00200009327002008250400000"; clusterCutArray[2] = "1111100053032200000"; mesonCutArray[2] = "0163103100000000"; //no M02 cut
		eventCutArray[ 3] = "80000013"; photonCutArray[ 3] = "00200009327002008250400000"; clusterCutArray[3] = "1112100053032230000"; mesonCutArray[3] = "0163103100000000"; //only modules with TRD infront
		eventCutArray[ 4] = "80000013"; photonCutArray[ 4] = "00200009327002008250400000"; clusterCutArray[4] = "1111300053032230000"; mesonCutArray[4] = "0163103100000000"; //no modules with TRD infront
	} else if (trainConfig == 9){ //PCM restriction in acceptance 
		eventCutArray[ 0] = "80000013"; photonCutArray[ 0] = "00202209327002008250400000"; clusterCutArray[0] = "1111100053032230000"; mesonCutArray[0] = "0163103100000000"; // PCM photons pointing to EMCAL loose
		eventCutArray[ 1] = "80000013"; photonCutArray[ 1] = "00204409327002008250400000"; clusterCutArray[1] = "1111100053032230000"; mesonCutArray[1] = "0163103100000000"; // PCM photons pointing to EMCAL tight
	} else if (trainConfig == 10){ 
		eventCutArray[ 0] = "80052013"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "1111100053062230000"; mesonCutArray[0] = "0163103100000000";
		eventCutArray[ 1] = "80052013"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "1111100053072230000"; mesonCutArray[1] = "0163103100000000";
		eventCutArray[ 2] = "80052013"; photonCutArray[ 2] = "00200009327002008250400000"; clusterCutArray[2] = "1111100053082230000"; mesonCutArray[2] = "0163103100000000";
		eventCutArray[ 3] = "80052013"; photonCutArray[ 3] = "00200009327002008250400000"; clusterCutArray[3] = "1111100053092230000"; mesonCutArray[3] = "0163103100000000";
		
	//************************************************ PHOS clusters **********************************************************	
	} else if (trainConfig == 31) {	// min energy = 0.3 GeV/c
		eventCutArray[ 0] = "80000013"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "2444400048033200000"; mesonCutArray[0] = "0163103100000000"; //standart cut, kINT7
		eventCutArray[ 1] = "80062013"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[1] = "0163103100000000"; //standard cut, kPHI7
	} else if (trainConfig == 32) { //PHOS
		eventCutArray[ 0] = "80000013"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "2444400047033200000"; mesonCutArray[0] = "0163103100000000";
		eventCutArray[ 1] = "80000013"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[1] = "0163103100000000";
		eventCutArray[ 2] = "80000013"; photonCutArray[ 2] = "00200009327002008250400000"; clusterCutArray[2] = "2444400049033200000"; mesonCutArray[2] = "0163103100000000";
	} else if (trainConfig == 33) { //PHOS
		eventCutArray[ 0] = "80000023"; photonCutArray[ 0] = "00200009327002008250400000"; clusterCutArray[0] = "2444400047033200000"; mesonCutArray[0] = "0163103100000000";
		eventCutArray[ 1] = "80000023"; photonCutArray[ 1] = "00200009327002008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[1] = "0163103100000000";
		eventCutArray[ 2] = "80000023"; photonCutArray[ 2] = "00200009327002008250400000"; clusterCutArray[2] = "2444400049033200000"; mesonCutArray[2] = "0163103100000000";
	} else {
		Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
	TList *ConvCutList = new TList();
	TList *ClusterCutList = new TList();
	TList *MesonCutList = new TList();

	TList *HeaderList = new TList();
	if (doWeightingPart==1) {
		TObjString *Header1 = new TObjString("pi0_1");
		HeaderList->Add(Header1);
	}
	if (doWeightingPart==2){
		TObjString *Header3 = new TObjString("eta_2");
		HeaderList->Add(Header3);
	}
	if (doWeightingPart==3) {
		TObjString *Header1 = new TObjString("pi0_1");
		HeaderList->Add(Header1);
		TObjString *Header3 = new TObjString("eta_2");
		HeaderList->Add(Header3);
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
