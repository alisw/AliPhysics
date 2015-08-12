void AddTask_GammaCalo_pPb(  
							Int_t 		trainConfig 				= 1,  								// change different set of cuts
							Int_t 		isMC   						= 0, 							// run MC
							Int_t 		enableQAMesonTask 			= 0, 								// enable QA in AliAnalysisTaskGammaConvV1
							Int_t 		enableQAClusterTask 		= 0, 								// enable additional QA task
							TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
							Int_t 		doWeightingPart 			= 0,  								// enable Weighting
							TString 	generatorName 				= "DPMJET",
                            TString 	cutnumberAODBranch 			= "800000006008400000001500000", 	// cutnumber for AOD branch
							Bool_t 		isUsingTHnSparse 			= kTRUE, 							// enable or disable usage of THnSparses for background estimation
							Int_t 		enableExtQA					= 0,								// enable QA(3), disabled (0)
							Bool_t 		enableTriggerMimicking		= kFALSE,							// enable trigger mimicking
							Bool_t 		enableTriggerOverlapRej		= kFALSE,							// enable trigger overlap rejection
							Float_t		maxFacPtHard				= 3									// maximum factor between hardest jet and ptHard generated
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
		Error(Form("AddTask_GammaCalo_pPb_%i",trainConfig), "No analysis manager found.");
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
	TString cutnumberEvent = "80000003";
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
	if (trainConfig == 5 || trainConfig == 6 ){ numberOfCuts = 5;}
	if (trainConfig == 7 || trainConfig == 8 || trainConfig == 32 || trainConfig == 33){ numberOfCuts = 1;}
	
	TString *eventCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// cluster cuts
	// 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
	// 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
	
	//************************************************ EMCAL clusters *************************************************
	if (trainConfig == 1){ // min energy = 0.3 GeV/c
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111100050022230000"; mesonCutArray[0] = "0163103100000000"; //standart cut, kINT7 // EMCAL clusters
        eventCutArray[ 1] = "80052013"; clusterCutArray[1] = "1111100050022230000"; mesonCutArray[1] = "0163103100000000"; //standard cut, kEMC7 // EMCAL clusters
	} else if (trainConfig == 2){  // min energy = 0.3 GeV/c
        eventCutArray[ 0] = "80083013"; clusterCutArray[0] = "1111100050022230000"; mesonCutArray[0] = "0163103100000000"; //standard cut, kEMCEG1 based on INT7 // EMCAL clusters
        eventCutArray[ 1] = "80085013"; clusterCutArray[1] = "1111100050022230000"; mesonCutArray[1] = "0163103100000000"; //standard cut, kEMCEG2 based on INT7 // EMCAL clusters
	} else if (trainConfig == 3){ // min energy = 0.4 GeV/c
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[0] = "0163103100000000"; //standart cut, kINT7 // EMCAL clusters
        eventCutArray[ 1] = "80052013"; clusterCutArray[1] = "1111100050032230000"; mesonCutArray[1] = "0163103100000000"; //standard cut, kEMC7 // EMCAL clusters
	} else if (trainConfig == 4){ // min energy = 0.4 GeV/c
        eventCutArray[ 0] = "80083013"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[0] = "0163103100000000"; //standard cut, kEMCEG1 based on INT7 // EMCAL clusters
        eventCutArray[ 1] = "80085013"; clusterCutArray[1] = "1111100050032230000"; mesonCutArray[1] = "0163103100000000"; //standard cut, kEMCEG2 based on INT7 // EMCAL clusters
	} else if (trainConfig == 5){ //EMCAL minEnergy variation
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111100050012230000"; mesonCutArray[0] = "0163103100000000"; //0.2 GeV/c
        eventCutArray[ 1] = "80000013"; clusterCutArray[1] = "1111100050022230000"; mesonCutArray[1] = "0163103100000000"; //0.3 GeV/c
        eventCutArray[ 2] = "80000013"; clusterCutArray[2] = "1111100050032230000"; mesonCutArray[2] = "0163103100000000"; //0.4 GeV/c default
        eventCutArray[ 3] = "80000013"; clusterCutArray[3] = "1111100050042230000"; mesonCutArray[3] = "0163103100000000"; //0.5 GeV/c
        eventCutArray[ 4] = "80000013"; clusterCutArray[4] = "1111100050052230000"; mesonCutArray[4] = "0163103100000000"; //0.6 GeV/c
	} else if (trainConfig == 6){ //EMCAL minNCells variation
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111100050031230000"; mesonCutArray[0] = "0163103100000000"; //n cells >= 1
        eventCutArray[ 1] = "80000013"; clusterCutArray[1] = "1111100050033230000"; mesonCutArray[1] = "0163103100000000"; //n cells >= 3
        eventCutArray[ 2] = "80000013"; clusterCutArray[2] = "1111100050032000000"; mesonCutArray[2] = "0163103100000000"; //no M02 cut
        eventCutArray[ 3] = "80000013"; clusterCutArray[3] = "1112100050032230000"; mesonCutArray[3] = "0163103100000000"; //only modules with TRD infront
        eventCutArray[ 4] = "80000013"; clusterCutArray[4] = "1111300050032230000"; mesonCutArray[4] = "0163103100000000"; //no modules with TRD infront
	} else if (trainConfig == 7){ // Validation EMCAL
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111100050032000000"; mesonCutArray[0] = "0163103100000000";
	} else if (trainConfig == 8){ // Validation EMCAL, only added signals
        eventCutArray[ 0] = "80000023"; clusterCutArray[0] = "1111100050032000000"; mesonCutArray[0] = "0163103100000000";
		
	//************************************************ PHOS clusters *************************************************
	} else if (trainConfig == 31) {	// min energy = 0.3 GeV/c
		eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "2444400040033200000"; mesonCutArray[0] = "0163103100000000"; //standart cut, kINT7 // PHOS clusters
		eventCutArray[ 1] = "80062013"; clusterCutArray[1] = "2444400040033200000"; mesonCutArray[1] = "0163103100000000"; //standard cut, kPHI7	// PHOS clusters	
	} else if (trainConfig == 32){ // Validation PHOS
		eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "2444400040053200000"; mesonCutArray[0] = "0163103100000000"; 	
	} else if (trainConfig == 33){ // Validation PHOS, only added signals
		eventCutArray[ 0] = "80000023"; clusterCutArray[0] = "2444400040053200000"; mesonCutArray[0] = "0163103100000000"; 	

		
	} else {
		Error(Form("GammaCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
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
	ClusterCutList->SetOwner(kTRUE);
	AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

	for(Int_t i = 0; i<numberOfCuts; i++){
		analysisEventCuts[i] = new AliConvEventCuts();   
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
