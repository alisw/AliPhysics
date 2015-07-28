void AddTask_GammaConvV1_PbPb2(  	Int_t 		trainConfig 				= 1,  								// change different set of cuts
									Bool_t 		isMC   						= kFALSE, 							// run MC 
									Int_t 		enableQAMesonTask 			= 0, 								// enable QA in AliAnalysisTaskGammaConvV1
									Int_t 		enableQAPhotonTask 			= 0, 								// enable additional QA task
									TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
									Bool_t 		doWeighting 				= kFALSE,  							// enable Weighting
									TString 	cutnumberAODBranch 			= "100000006008400000001500000", 	// cutnumber with which AODs have been filtered
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
	
	//=========  Set Cutnumber for V0Reader ================================
	TString cutnumberPhoton = "00000008400100001500000000";
	TString cutnumberEvent = "1000000";
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

		AliLog::SetGlobalLogLevel(AliLog::kInfo);

		//connect input V0Reader
		mgr->AddTask(fV0ReaderV1);
		mgr->ConnectInput(fV0ReaderV1,0,cinput);

	}

	//================================================
	//========= Add task to the ANALYSIS manager =====
	//================================================
	AliAnalysisTaskGammaConvV1 *task=NULL;
	task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 1;

	TString *eventCutArray = new TString[numberOfCuts];
	TString *photonCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	if (trainConfig == 1){ 
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; 
	} else if (trainConfig == 2) { 
		eventCutArray[ 0] = "6120001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; 
	} else if (trainConfig == 3) { 
		eventCutArray[ 0] = "5010001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; 
	} else if (trainConfig == 4) { 
		eventCutArray[ 0] = "5020001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000";    
	} else if (trainConfig == 5) { 
		eventCutArray[ 0] = "5120001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000";    
	} else if (trainConfig == 6) { 
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000";       
	} else if (trainConfig == 7) {    
		eventCutArray[ 0] = "5460001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000"; 
	} else if (trainConfig == 8) {    
		eventCutArray[ 0] = "5480001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";    
	} else if (trainConfig == 9) {    
		eventCutArray[ 0] = "5450001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000"; 
	} else if (trainConfig == 10) { 
		eventCutArray[ 0] = "5560001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";
	} else if (trainConfig == 11) { 
		eventCutArray[ 0] = "5680001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";    
	} else if (trainConfig == 12) { 
		eventCutArray[ 0] = "5670001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000"; 
	} else if (trainConfig == 13) { 
		eventCutArray[ 0] = "5780001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000"; 
	} else if (trainConfig == 14) { 
		eventCutArray[ 0] = "4690001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";
	} else if (trainConfig == 15) { 
		eventCutArray[ 0] = "5890001"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";    
	} else  if (trainConfig == 16){ 
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; 
	} else if (trainConfig == 17) { 
		eventCutArray[ 0] = "6120001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; 
	} else if (trainConfig == 18) { 
		eventCutArray[ 0] = "5010001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; 
	} else if (trainConfig == 19) { 
		eventCutArray[ 0] = "5020001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";    
	} else if (trainConfig == 20) { 
		eventCutArray[ 0] = "5120001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";    
	} else if (trainConfig == 21) { 
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";       
	} else if (trainConfig == 22) {    
		eventCutArray[ 0] = "5460001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; 
	} else if (trainConfig == 23) {    
		eventCutArray[ 0] = "5480001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";    
	} else if (trainConfig == 24) {    
		eventCutArray[ 0] = "5450001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; 
	} else if (trainConfig == 25) { 
		eventCutArray[ 0] = "5560001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";
	} else if (trainConfig == 26) { 
		eventCutArray[ 0] = "5680001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";    
	} else if (trainConfig == 27) { 
		eventCutArray[ 0] = "5670001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; 
	} else if (trainConfig == 28) { 
		eventCutArray[ 0] = "5780001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; 
	} else if (trainConfig == 29) { 
		eventCutArray[ 0] = "4690001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";
	} else if (trainConfig == 30) { 
		eventCutArray[ 0] = "5890001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";    
	} else if (trainConfig == 31) { 
		eventCutArray[ 0] = "5080001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";    
	} else if (trainConfig == 32) { 
		eventCutArray[ 0] = "5250001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";       
	} else if (trainConfig == 33) { 
		eventCutArray[ 0] = "5350001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";       
	} else if (trainConfig == 34) { 
		eventCutArray[ 0] = "5450001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";       
	} else if (trainConfig == 35) { 
		eventCutArray[ 0] = "5340001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";       
	} else if (trainConfig == 36) { 
		eventCutArray[ 0] = "5230001"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000";       
	} else {
		Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
	TList *ConvCutList = new TList();
	TList *MesonCutList = new TList();

	TList *HeaderList = new TList();
	TObjString *Header1 = new TObjString("BOX");
	HeaderList->Add(Header1);
	//    TObjString *Header3 = new TObjString("eta_2");
	//    HeaderList->Add(Header3);
	
	EventCutList->SetOwner(kTRUE);
	AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
	ConvCutList->SetOwner(kTRUE);
	AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

	for(Int_t i = 0; i<numberOfCuts; i++){
		analysisEventCuts[i] = new AliConvEventCuts();
		analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
		EventCutList->Add(analysisEventCuts[i]);
		analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
		
		analysisCuts[i] = new AliConversionPhotonCuts();
		analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());
		ConvCutList->Add(analysisCuts[i]);
		analysisCuts[i]->SetFillCutHistograms("",kFALSE);

		analysisMesonCuts[i] = new AliConversionMesonCuts();
		analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
		MesonCutList->Add(analysisMesonCuts[i]);
		analysisMesonCuts[i]->SetFillCutHistograms("");
		
		analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
	}

	task->SetEventCutList(numberOfCuts,EventCutList);
	task->SetConversionCutList(numberOfCuts,ConvCutList);
	task->SetMesonCutList(numberOfCuts,MesonCutList);
	task->SetMoveParticleAccordingToVertex(kTRUE);
	task->SetDoMesonAnalysis(kTRUE);
	task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
	task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
	task->SetDoPlotVsCentrality(kTRUE);

	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
