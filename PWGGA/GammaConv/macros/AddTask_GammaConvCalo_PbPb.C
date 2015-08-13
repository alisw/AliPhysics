void AddTask_GammaConvCalo_PbPb(   	Int_t 		trainConfig 				= 1,  								// change different set of cuts
									Int_t 		isMC 						= 0, 							// run MC 
									Int_t 		enableQAMesonTask 			= 0, 								// enable QA in AliAnalysisTaskGammaConvV1
									Int_t 		enableQAPhotonTask 			= 0, 								// enable additional QA task
									TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
									Int_t 		headerSelectionInt 			= 0,  								// 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
									TString 	cutnumberAODBranch 			= "100000006008400000001500000",
									TString 	periodName 					= "LHC13d2",  						// name of the period for added signals and weighting
									Bool_t 		doWeighting 				= kFALSE,  							// enable Weighting
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

	Int_t isHeavyIon = 1;
	
	// ================== GetAnalysisManager ===============================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error(Form("AddTask_GammaConvCalo_PbPb_%i",trainConfig), "No analysis manager found.");
		return ;
	}

	// ================== GetInputEventHandler =============================
	AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
	//**********************************************************************************************
	//**********************************************************************************************
		inputHandler->SetInactiveBranches("AliESDFMD"); // Disable FMD branch, see ALIROOT-6222
	//**********************************************************************************************
	//**********************************************************************************************
	Bool_t isMCForOtherTasks = kFALSE;
	if (isMC > 0) isMCForOtherTasks = kTRUE;
	
	//========= Add PID Reponse to ANALYSIS manager ====
	if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
		gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
		AddTaskPIDResponse(isMCForOtherTasks);
	}
	
	//=========  Set Cutnumber for V0Reader ================================
	TString cutnumberPhoton = "00000008400100001500000000";
	TString cutnumberEvent = "10000003";
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
	Int_t numberOfCuts = 5;

	TString *eventCutArray = new TString[numberOfCuts];
	TString *photonCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];
  
	if (trainConfig == 1){ // EMCAL clusters
		eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[ 0] = "0163103100000000"; // 0-5%
		eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "1111100050032230000"; mesonCutArray[ 1] = "0163103100000000"; // 5-10%
		eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "1111100050032230000"; mesonCutArray[ 2] = "0163103100000000"; // 0-10%
		eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "1111100050032230000"; mesonCutArray[ 3] = "0163103100000000"; // 20-40%
		eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "1111100050032230000"; mesonCutArray[ 4] = "0163103100000000"; // 20-50%
	} else if (trainConfig == 2){ // EMCAL clusters
		eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[ 0] = "0163103100000000"; // 0-5%
		eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "1111100050032230000"; mesonCutArray[ 1] = "0163103100000000"; // 5-10%
		eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "1111100050032230000"; mesonCutArray[ 2] = "0163103100000000"; // 0-10%
		eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "1111100050032230000"; mesonCutArray[ 3] = "0163103100000000"; // 10-20%
		eventCutArray[ 4] = "52400013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "1111100050032230000"; mesonCutArray[ 4] = "0163103100000000"; // 20-40%
	} else if (trainConfig == 3){ // EMCAL clusters
		eventCutArray[ 0] = "54600013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[ 0] = "0163103100000000"; // 40-60%
		eventCutArray[ 1] = "56800013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "1111100050032230000"; mesonCutArray[ 1] = "0163103100000000"; // 60-80%
		eventCutArray[ 2] = "52600013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "1111100050032230000"; mesonCutArray[ 2] = "0163103100000000"; // 20-60%
		eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "1111100050032230000"; mesonCutArray[ 3] = "0163103100000000"; // 40-80%
		eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "1111100050032230000"; mesonCutArray[ 4] = "0163103100000000"; // 20-50%
	} else if (trainConfig == 4){ // EMCAL clusters
		eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "1111100053032230000"; mesonCutArray[ 0] = "0163103100000000"; // 0-5%
		eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "1111100053032230000"; mesonCutArray[ 1] = "0163103100000000"; // 5-10%
		eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "1111100053032230000"; mesonCutArray[ 2] = "0163103100000000"; // 0-10%
		eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "1111100053032230000"; mesonCutArray[ 3] = "0163103100000000"; // 20-40%
		eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "1111100053032230000"; mesonCutArray[ 4] = "0163103100000000"; // 20-50%
	} else if (trainConfig == 5){ // EMCAL clusters
		eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "1111100053032230000"; mesonCutArray[ 0] = "0163103100000000"; // 0-5%
		eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "1111100053032230000"; mesonCutArray[ 1] = "0163103100000000"; // 5-10%
		eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "1111100053032230000"; mesonCutArray[ 2] = "0163103100000000"; // 0-10%
		eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "1111100053032230000"; mesonCutArray[ 3] = "0163103100000000"; // 10-20%
		eventCutArray[ 4] = "52400013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "1111100053032230000"; mesonCutArray[ 4] = "0163103100000000"; // 20-40%
	} else if (trainConfig == 6){ // EMCAL clusters
		eventCutArray[ 0] = "54600013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "1111100053032230000"; mesonCutArray[ 0] = "0163103100000000"; // 40-60%
		eventCutArray[ 1] = "56800013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "1111100053032230000"; mesonCutArray[ 1] = "0163103100000000"; // 60-80%
		eventCutArray[ 2] = "52600013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "1111100053032230000"; mesonCutArray[ 2] = "0163103100000000"; // 20-60%
		eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "1111100053032230000"; mesonCutArray[ 3] = "0163103100000000"; // 40-80%
		eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "1111100053032230000"; mesonCutArray[ 4] = "0163103100000000"; // 20-50%
	} else if (trainConfig == 31){ // PHOS clusters
		eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "2444400048033200000"; mesonCutArray[ 0] = "0163103100000000"; // 0-5%
		eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[ 1] = "0163103100000000"; // 5-10%
		eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "2444400048033200000"; mesonCutArray[ 2] = "0163103100000000"; // 0-10%
		eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "2444400048033200000"; mesonCutArray[ 3] = "0163103100000000"; // 20-40%
		eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "2444400048033200000"; mesonCutArray[ 4] = "0163103100000000"; // 20-50%
	} else if (trainConfig == 32){ // PHOS clusters
		eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "2444400048033200000"; mesonCutArray[ 0] = "0163103100000000"; // 0-5%
		eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[ 1] = "0163103100000000"; // 5-10%
		eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "2444400048033200000"; mesonCutArray[ 2] = "0163103100000000"; // 0-10%
		eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "2444400048033200000"; mesonCutArray[ 3] = "0163103100000000"; // 10-20%
		eventCutArray[ 4] = "52400013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "2444400048033200000"; mesonCutArray[ 4] = "0163103100000000"; // 20-40%
	} else if (trainConfig == 33){ // PHOS clusters
		eventCutArray[ 0] = "54600013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "2444400048033200000"; mesonCutArray[ 0] = "0163103100000000"; // 40-60%
		eventCutArray[ 1] = "56800013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[ 1] = "0163103100000000"; // 60-80%
		eventCutArray[ 2] = "52600013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "2444400048033200000"; mesonCutArray[ 2] = "0163103100000000"; // 20-60%
		eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "2444400048033200000"; mesonCutArray[ 3] = "0163103100000000"; // 40-80%
		eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "2444400048033200000"; mesonCutArray[ 4] = "0163103100000000"; // 20-50%
	} else if (trainConfig == 34){ // PHOS clusters
		eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "2444400048033200000"; mesonCutArray[ 0] = "0163103100000000"; // 0-5%
		eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[ 1] = "0163103100000000"; // 5-10%
		eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "2444400048033200000"; mesonCutArray[ 2] = "0163103100000000"; // 0-10%
		eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "2444400048033200000"; mesonCutArray[ 3] = "0163103100000000"; // 20-40%
		eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "2444400048033200000"; mesonCutArray[ 4] = "0163103100000000"; // 20-50%
	} else if (trainConfig == 35){ // PHOS clusters
		eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "2444400048033200000"; mesonCutArray[ 0] = "0163103100000000"; // 0-5%
		eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[ 1] = "0163103100000000"; // 5-10%
		eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "2444400048033200000"; mesonCutArray[ 2] = "0163103100000000"; // 0-10%
		eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "2444400048033200000"; mesonCutArray[ 3] = "0163103100000000"; // 10-20%
		eventCutArray[ 4] = "52400013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "2444400048033200000"; mesonCutArray[ 4] = "0163103100000000"; // 20-40%
	} else if (trainConfig == 36){ // PHOS clusters
		eventCutArray[ 0] = "54600013"; photonCutArray[ 0] = "00200009297002008250400000"; clusterCutArray[0] = "2444400048033200000"; mesonCutArray[ 0] = "0163103100000000"; // 40-60%
		eventCutArray[ 1] = "56800013"; photonCutArray[ 1] = "00200009297002008250400000"; clusterCutArray[1] = "2444400048033200000"; mesonCutArray[ 1] = "0163103100000000"; // 60-80%
		eventCutArray[ 2] = "52600013"; photonCutArray[ 2] = "00200009297002008250400000"; clusterCutArray[2] = "2444400048033200000"; mesonCutArray[ 2] = "0163103100000000"; // 20-60%
		eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "00200009297002008250400000"; clusterCutArray[3] = "2444400048033200000"; mesonCutArray[ 3] = "0163103100000000"; // 40-80%
		eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008250400000"; clusterCutArray[4] = "2444400048033200000"; mesonCutArray[ 4] = "0163103100000000"; // 20-50%
	} else {
		Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
	TList *ConvCutList = new TList();
	TList *ClusterCutList = new TList();
	TList *MesonCutList = new TList();

	TList *HeaderList = new TList();
	if (periodName.CompareTo("LHC13d2")==0){
		TObjString *Header1 = new TObjString("pi0_1");
		HeaderList->Add(Header1);
	//    TObjString *Header3 = new TObjString("eta_2");
	//    HeaderList->Add(Header3);

	} else if (periodName.CompareTo("LHC12a17x_fix")==0){
		TObjString *Header1 = new TObjString("PARAM");
		HeaderList->Add(Header1);
	} else if (periodName.CompareTo("LHC14a1a")==0){
		if (headerSelectionInt == 1){ 
			TObjString *Header1 = new TObjString("pi0_1");
			HeaderList->Add(Header1);
		} else if (headerSelectionInt == 2){
			TObjString *Header1 = new TObjString("eta_2");
			HeaderList->Add(Header1);
		} else {
			TObjString *Header1 = new TObjString("pi0_1");
			HeaderList->Add(Header1);
			TObjString *Header2 = new TObjString("eta_2");
			HeaderList->Add(Header2);
		}  
	} else if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0){
		TObjString *Header1 = new TObjString("BOX");
		HeaderList->Add(Header1);
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
// 		if ( trainConfig == 1){
// 			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
// 				if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
// 				if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
// 				if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
// 				if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
// 				if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
// 			}	
// 		} 
		analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
		if (periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
			if (headerSelectionInt == 1) analysisEventCuts[i]->SetAddedSignalPDGCode(111);
			if (headerSelectionInt == 2) analysisEventCuts[i]->SetAddedSignalPDGCode(221);
		}
		EventCutList->Add(analysisEventCuts[i]);
		analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

		analysisCuts[i] = new AliConversionPhotonCuts();
		analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());
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
