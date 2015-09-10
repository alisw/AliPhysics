void AddTask_GammaCalo_pPb(  
							Int_t 		trainConfig 				= 1,  								// change different set of cuts
							Int_t 		isMC   						= 0, 								// run MC
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
							Float_t		maxFacPtHard				= 3,									// maximum factor between hardest jet and ptHard generated
							TString		periodNameV0Reader			= ""
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
		if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
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
	// change to 1 cuts per cutselection
	if (trainConfig == 7 	|| trainConfig == 8 	|| trainConfig == 32 	|| trainConfig == 33 || trainConfig == 44)
		numberOfCuts = 1;
	// change to 3 cuts per cutselection
	if (trainConfig == 14 || trainConfig == 42)
		numberOfCuts = 3;
	// change to 4 cuts per cutselection
	if (trainConfig == 9 	|| trainConfig == 10	|| trainConfig == 11	|| trainConfig == 12 	|| trainConfig == 13 || trainConfig == 40 || trainConfig == 41 || trainConfig == 43)
		numberOfCuts = 4;
	// change to 5 cuts per cutselection
	if (trainConfig == 5 	|| trainConfig == 6 )
		numberOfCuts = 5;
	
	TString *eventCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// cluster cuts
	// 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
	// 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
	
	//************************************************ EMCAL clusters *************************************************
	if (trainConfig == 1){ // min energy = 0.3 GeV/c
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111181050022230000"; mesonCutArray[0] = "0163103100000050"; //standart cut, kINT7 // EMCAL clusters
        eventCutArray[ 1] = "80052013"; clusterCutArray[1] = "1111181050022230000"; mesonCutArray[1] = "0163103100000050"; //standard cut, kEMC7 // EMCAL clusters
	} else if (trainConfig == 2){  // min energy = 0.3 GeV/c
        eventCutArray[ 0] = "80083013"; clusterCutArray[0] = "1111181050022230000"; mesonCutArray[0] = "0163103100000050"; //standard cut, kEMCEG1 based on INT7 // EMCAL clusters
        eventCutArray[ 1] = "80085013"; clusterCutArray[1] = "1111181050022230000"; mesonCutArray[1] = "0163103100000050"; //standard cut, kEMCEG2 based on INT7 // EMCAL clusters
	} else if (trainConfig == 3){ // min energy = 0.4 GeV/c
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111181050032230000"; mesonCutArray[0] = "0163103100000050"; //standart cut, kINT7 // EMCAL clusters
        eventCutArray[ 1] = "80052013"; clusterCutArray[1] = "1111181050032230000"; mesonCutArray[1] = "0163103100000050"; //standard cut, kEMC7 // EMCAL clusters
	} else if (trainConfig == 4){ // min energy = 0.4 GeV/c
        eventCutArray[ 0] = "80083013"; clusterCutArray[0] = "1111181050032230000"; mesonCutArray[0] = "0163103100000050"; //standard cut, kEMCEG1 based on INT7 // EMCAL clusters
        eventCutArray[ 1] = "80085013"; clusterCutArray[1] = "1111181050032230000"; mesonCutArray[1] = "0163103100000050"; //standard cut, kEMCEG2 based on INT7 // EMCAL clusters
	} else if (trainConfig == 5){ //EMCAL minEnergy variation
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111181050012230000"; mesonCutArray[0] = "0163103100000050"; //0.2 GeV/c
        eventCutArray[ 1] = "80000013"; clusterCutArray[1] = "1111181050022230000"; mesonCutArray[1] = "0163103100000050"; //0.3 GeV/c
        eventCutArray[ 2] = "80000013"; clusterCutArray[2] = "1111181050032230000"; mesonCutArray[2] = "0163103100000050"; //0.4 GeV/c default
        eventCutArray[ 3] = "80000013"; clusterCutArray[3] = "1111181050042230000"; mesonCutArray[3] = "0163103100000050"; //0.5 GeV/c
        eventCutArray[ 4] = "80000013"; clusterCutArray[4] = "1111181050052230000"; mesonCutArray[4] = "0163103100000050"; //0.6 GeV/c
	} else if (trainConfig == 6){ //EMCAL minNCells variation
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111181050031230000"; mesonCutArray[0] = "0163103100000050"; //n cells >= 1
        eventCutArray[ 1] = "80000013"; clusterCutArray[1] = "1111181050033230000"; mesonCutArray[1] = "0163103100000050"; //n cells >= 3
        eventCutArray[ 2] = "80000013"; clusterCutArray[2] = "1111181050032000000"; mesonCutArray[2] = "0163103100000050"; //no M02 cut
        eventCutArray[ 3] = "80000013"; clusterCutArray[3] = "1112181050032230000"; mesonCutArray[3] = "0163103100000050"; //only modules with TRD infront
        eventCutArray[ 4] = "80000013"; clusterCutArray[4] = "1111381050032230000"; mesonCutArray[4] = "0163103100000050"; //no modules with TRD infront
	} else if (trainConfig == 7){ // Validation EMCAL
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111181050032000000"; mesonCutArray[0] = "0163103100000050";
	} else if (trainConfig == 8){ // Validation EMCAL, only added signals
        eventCutArray[ 0] = "80000023"; clusterCutArray[0] = "1111181050032000000"; mesonCutArray[0] = "0163103100000050";
	} else if (trainConfig == 9){ // non linearity variations INT7
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[0] = "0163103100000050"; // non nonlinearity
        eventCutArray[ 1] = "80000013"; clusterCutArray[1] = "1111101050032230000"; mesonCutArray[1] = "0163103100000050"; // kSDM
        eventCutArray[ 2] = "80000013"; clusterCutArray[2] = "1111181050032230000"; mesonCutArray[2] = "0163103100000050"; // conv calo
        eventCutArray[ 3] = "80000013"; clusterCutArray[3] = "1111182050032230000"; mesonCutArray[3] = "0163103100000050"; // calo
	} else if (trainConfig == 10){ // non linearity variations EMC7
        	eventCutArray[ 0] = "80052013"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[0] = "0163103100000050"; // non nonlinearity
		eventCutArray[ 1] = "80052013"; clusterCutArray[1] = "1111101050032230000"; mesonCutArray[1] = "0163103100000050"; // kSDM
		eventCutArray[ 2] = "80052013"; clusterCutArray[2] = "1111181050032230000"; mesonCutArray[2] = "0163103100000050"; // conv calo
		eventCutArray[ 3] = "80052013"; clusterCutArray[3] = "1111182050032230000"; mesonCutArray[3] = "0163103100000050"; // calo
	} else if (trainConfig == 11){ /// non linearity variations EG2
        eventCutArray[ 0] = "80085013"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[0] = "0163103100000050"; // non nonlinearity
		eventCutArray[ 1] = "80085013"; clusterCutArray[1] = "1111101050032230000"; mesonCutArray[1] = "0163103100000050"; // kSDM
		eventCutArray[ 2] = "80085013"; clusterCutArray[2] = "1111181050032230000"; mesonCutArray[2] = "0163103100000050"; // conv calo
		eventCutArray[ 3] = "80085013"; clusterCutArray[3] = "1111182050032230000"; mesonCutArray[3] = "0163103100000050"; // calo
	} else if (trainConfig == 12){ // non linearity variations EG1
		eventCutArray[ 0] = "80083013"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[0] = "0163103100000050"; // non nonlinearity
		eventCutArray[ 1] = "80083013"; clusterCutArray[1] = "1111101050032230000"; mesonCutArray[1] = "0163103100000050"; // kSDM
		eventCutArray[ 2] = "80083013"; clusterCutArray[2] = "1111181050032230000"; mesonCutArray[2] = "0163103100000050"; // conv calo
		eventCutArray[ 3] = "80083013"; clusterCutArray[3] = "1111182050032230000"; mesonCutArray[3] = "0163103100000050"; // calo
	} else if (trainConfig == 13){ // no non linearity
        eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[0] = "0163103100000050"; // kINT7 // EMCAL clusters
        eventCutArray[ 1] = "80052013"; clusterCutArray[1] = "1111100050032230000"; mesonCutArray[1] = "0163103100000050"; // kEMC7 // EMCAL clusters
        eventCutArray[ 2] = "80083013"; clusterCutArray[2] = "1111100050032230000"; mesonCutArray[2] = "0163103100000050"; // kEMCEG1 based on INT7 // EMCAL clusters
        eventCutArray[ 3] = "80085013"; clusterCutArray[3] = "1111100050032230000"; mesonCutArray[3] = "0163103100000050"; // kEMCEG2 based on INT7 // EMCAL clusters
	} else 	if(trainConfig == 14){ // variation opening angle
		eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111181050022230000"; mesonCutArray[0] = "0163103100000050"; // standard
		eventCutArray[ 1] = "80000013"; clusterCutArray[1] = "1111181050022230000"; mesonCutArray[1] = "0163103100000060"; // 2 EMCal cell diagonals
		eventCutArray[ 2] = "80000013"; clusterCutArray[2] = "1111181050022230000"; mesonCutArray[2] = "0163103100000040"; // 0.75 EMCal cell diagonals

	// SYSTEMATIC STUDY NEUTRAl MESON MEASUREMENTS MIKE SAS 10-09-2015
	} else 	if(trainConfig == 40){ // default cutstring and first set of variations nonlinearity
		eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111181050032230000"; mesonCutArray[0] = "0163403100000050"; // default
		eventCutArray[ 1] = "80000013"; clusterCutArray[1] = "1111182050032230000"; mesonCutArray[1] = "0163403100000050"; // calo nonlinearity variation
		eventCutArray[ 2] = "80000013"; clusterCutArray[2] = "1111183050032230000"; mesonCutArray[2] = "0163403100000050"; // calo nonlinearity variation
		eventCutArray[ 3] = "80000013"; clusterCutArray[3] = "1111184050032230000"; mesonCutArray[3] = "0163403100000050"; // calo nonlinearity variation
	} else 	if(trainConfig == 41){ // second set of variations CLUSTER ====>> SETTING TENDER: M_seed=100MeV
		eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111181050012230000"; mesonCutArray[0] = "0163403100000050"; // min energy cluster variation 1	200 MeV
		eventCutArray[ 1] = "80000013"; clusterCutArray[1] = "1111181050022230000"; mesonCutArray[1] = "0163403100000050"; // min energy cluster variation 2	300 MeV
		eventCutArray[ 2] = "80000013"; clusterCutArray[2] = "1111181050042230000"; mesonCutArray[2] = "0163403100000050"; // min energy cluster variation 3	500 MeV
		eventCutArray[ 3] = "80000013"; clusterCutArray[3] = "1111181050052230000"; mesonCutArray[3] = "0163403100000050"; // min energy cluster variation 4	600 MeV
	} else 	if(trainConfig == 42){ // third set of variations CLUSTER
		eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111181050032030000"; mesonCutArray[0] = "0163403100000050"; // min/max M02	0<M<0.5
		eventCutArray[ 1] = "80000013"; clusterCutArray[1] = "1111181050032200000"; mesonCutArray[1] = "0163403100000050"; // min/max M02	0.1<M<100
		eventCutArray[ 2] = "80000013"; clusterCutArray[2] = "1111181050031230000"; mesonCutArray[2] = "0163403100000050"; // min number of cells variation 1	1 cell
	} else 	if(trainConfig == 43){ // third set of variations MESON
		eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111181050032230000"; mesonCutArray[0] = "0163303100000050"; // rapidity variation	y<0.6
		eventCutArray[ 1] = "80000013"; clusterCutArray[1] = "1111181050032230000"; mesonCutArray[1] = "0163103100000050"; // rapidity variation	y<0.8	
		eventCutArray[ 2] = "80000013"; clusterCutArray[2] = "1111181050032230000"; mesonCutArray[2] = "0163406100000050"; // alpha meson variation 1 	0<alpha<0.8
		eventCutArray[ 3] = "80000013"; clusterCutArray[3] = "1111181050032230000"; mesonCutArray[3] = "0163405100000050"; // alpha meson variation 2	0<alpha<0.75
	} else 	if(trainConfig == 44){ // default cutstring for different tender settings
		eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "1111181050032230000"; mesonCutArray[0] = "0163403100000050"; // default


	//************************************************ PHOS clusters *************************************************
	} else if (trainConfig == 31) {	// min energy = 0.3 GeV/c
		eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "2444400040033200000"; mesonCutArray[0] = "0163103100000050"; //standart cut, kINT7 // PHOS clusters
		eventCutArray[ 1] = "80062013"; clusterCutArray[1] = "2444400040033200000"; mesonCutArray[1] = "0163103100000050"; //standard cut, kPHI7	// PHOS clusters	
	} else if (trainConfig == 32){ // Validation PHOS
		eventCutArray[ 0] = "80000013"; clusterCutArray[0] = "2444400040053200000"; mesonCutArray[0] = "0163103100000050"; 	
	} else if (trainConfig == 33){ // Validation PHOS, only added signals
		eventCutArray[ 0] = "80000023"; clusterCutArray[0] = "2444400040053200000"; mesonCutArray[0] = "0163103100000050"; 	

		
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
