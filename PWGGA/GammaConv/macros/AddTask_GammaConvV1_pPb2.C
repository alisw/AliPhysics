void AddTask_GammaConvV1_pPb2(  Int_t 		trainConfig 				= 1,  								// change different set of cuts
								Int_t 		isMC   						= 0,								// run MC
								Int_t 		enableQAMesonTask 			= 0, 								// enable QA in AliAnalysisTaskGammaConvV1
								Int_t 		enableQAPhotonTask 			= 0, 								// enable additional QA task
								TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
								Bool_t 		doWeightingPart 			= kFALSE,  							// enable Weighting
								TString 	generatorName 				= "DPMJET",							// generator Name	
								TString 	cutnumberAODBranch 			= "800000006008400000001500000",	// cutnumber for AOD branch
								Bool_t 		enableV0findingEffi 		= kFALSE,							// enables V0finding efficiency histograms
								Bool_t		enablePlotVsCentrality		= kFALSE
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
		Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
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
	
	//=========  Set Cutnumber for V0Reader ================================
	TString cutnumberPhoton = "06000008400100001500000000";
	TString cutnumberEvent = "80000103";
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
			if(trainConfig==15 ||trainConfig==16 ||trainConfig==17  ||trainConfig==18  ){
			  fCuts->SetDodEdxSigmaCut(kFALSE);
			}
		}
		if(inputHandler->IsA()==AliAODInputHandler::Class()){
		// AOD mode
			cout << "AOD handler: adding " << cutnumberAODBranch.Data() << " as conversion branch" << endl;
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
	//            find input container
	AliAnalysisTaskGammaConvV1 *task=NULL;
	task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 1;
	
	TString *eventCutArray = new TString[numberOfCuts];
	TString *photonCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];
	Bool_t doEtaShiftIndCuts = kFALSE;
	TString stringShift = "";
	
	if(trainConfig == 1){
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009317200003290000000"; mesonCutArray[ 0] = "0162103500900000"; 
	} else if (trainConfig == 2) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";
	} else if (trainConfig == 3) {   
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";      
	} else if (trainConfig == 4) {   
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";         
	} else if (trainConfig == 5) {   
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";         
	} else if (trainConfig == 6) {   
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";      
	} else if (trainConfig == 7) {   
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";         
	} else if (trainConfig == 8) {   
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00900009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";    //RCut 7.5cm   
	} else if (trainConfig == 9) {   
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00500009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";    //RCut 10cm     
	} else if (trainConfig == 10) {   
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00800009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";    //RCut 12.5cm    
	} else if (trainConfig == 11) {   
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00600009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";    //RCut 20cm    
	} else if (trainConfig == 12) {   
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00700009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";    //RCut 35cm    
	} else if (trainConfig == 13) {   
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";    //standard, opening angle =0.0   
	} else if (trainConfig == 14) {   
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";   //standard, opening angle =0.0  
	} else if (trainConfig == 15) {   
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009000000008260404000"; mesonCutArray[ 0] = "0162101500900000";   //standard, nodEdx, alpha 1 MB
	} else if (trainConfig == 16) {   
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009000000008260404000"; mesonCutArray[ 0] = "0162101500900000";   //standard, nodEdx, alpha 1 AddSignal Pi0
	} else if (trainConfig == 17) {   
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009000000008260404000"; mesonCutArray[ 0] = "0162101500000000";   //standard, nodEdx, alpha 1 MB Eta 
	} else if (trainConfig == 18) {   
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009000000008260404000"; mesonCutArray[ 0] = "0162101500000000";   //standard, nodEdx, alpha 1 AddSignal Eta
	} else if (trainConfig == 19) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260404000"; mesonCutArray[ 0] = "0162101500000000";   //standard, alpha 1
	} else {
		Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
	TList *ConvCutList = new TList();
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
	
	
	Bool_t doWeighting = kFALSE;
	if (doWeightingPart == 1 || doWeightingPart == 2 || doWeightingPart == 3) doWeighting = kTRUE;
	
	EventCutList->SetOwner(kTRUE);
	AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
	ConvCutList->SetOwner(kTRUE);
	AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];
	
	if (doWeighting) Printf("weighting has been switched on");
	
	for(Int_t i = 0; i<numberOfCuts; i++){
		
		analysisEventCuts[i] = new AliConvEventCuts();
		if ( trainConfig == 13 || trainConfig == 15 || trainConfig == 17 || trainConfig == 19){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A",
																				 "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A",
																				 "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
				}
			}
		}   
		if ( trainConfig == 14 || trainConfig == 16 || trainConfig == 18){
			if (doWeighting){
				analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A",
																			 "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
			}
			
		}   
 
		analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
		if (doEtaShiftIndCuts) {
			analysisEventCuts[i]->DoEtaShift(doEtaShiftIndCuts);
			analysisEventCuts[i]->SetEtaShift(stringShift);
		}
		
		EventCutList->Add(analysisEventCuts[i]);
		analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
		
		analysisCuts[i] = new AliConversionPhotonCuts();
		analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());
		analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
		if (trainConfig == 15 || trainConfig==16 || trainConfig==17  || trainConfig==18) {
		        analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);
		}
		ConvCutList->Add(analysisCuts[i]);
		analysisCuts[i]->SetFillCutHistograms("",kFALSE);
		
		analysisMesonCuts[i] = new AliConversionMesonCuts();
		if (trainConfig ==13 || trainConfig ==14){
		  analysisMesonCuts[i]->SetOpeningAngleCut(0.000);		  
		}
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
	if (trainConfig ==13 || trainConfig ==14){
	        task->SetDoTHnSparse(0);
	}
	task->SetDoPlotVsCentrality(enablePlotVsCentrality);

	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));
	
	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);
	
	return;
}
