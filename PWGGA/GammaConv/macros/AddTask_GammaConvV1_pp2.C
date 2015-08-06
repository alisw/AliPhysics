void AddTask_GammaConvV1_pp2(  	Int_t 		trainConfig 				= 1,  								//change different set of cuts
								Bool_t 		isMC   						= kFALSE, 							//run MC 
								Int_t 		enableQAMesonTask 			= 0, 								//enable QA in AliAnalysisTaskGammaConvV1
								Int_t 		enableQAPhotonTask 			= 0, 								// enable additional QA task
								TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
								TString 	cutnumberAODBranch 			= "000000006008400001001500000", 	// cutnumber with which AODs have been filtered 
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
	
	//=========  Set Cutnumber for V0Reader ================================
		TString cutnumberPhoton = "00200008400000002200000000";
		TString cutnumberEvent = "00000003"; 
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
				if (trainConfig == 21){
				  fCuts->SetDodEdxSigmaCut(kFALSE);
				}
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
		eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , only boxes
	} else if (trainConfig == 2) { 
		eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD, V0AND , only boxes
	} else if (trainConfig == 3) { 
		eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009326000003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Gamma pp 2-76TeV , only boxes
	} else if (trainConfig == 4) { 
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , only Minbias MC
	} else if (trainConfig == 5) { 
		eventCutArray[ 0] = "00010113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD, V0AND
	} else if (trainConfig == 6) { 
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009326000003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Gamma pp 2-76TeV
	} else if (trainConfig == 7) {    
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , only Minbias MC
	} else if (trainConfig == 8) {    
		eventCutArray[ 0] = "00013113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND , only Minbias MC
	} else if (trainConfig == 9) {    
		eventCutArray[ 0] = "00003123"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , only Boxes MC
	} else if (trainConfig == 10) { 
		eventCutArray[ 0] = "00013123"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND, only Boxes MC
	} else if (trainConfig == 11) { 
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , all photon qualities
	} else if (trainConfig == 12) { 
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00700009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , all photon qualities, min R = 35 cm
	} else if (trainConfig == 13) { 
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , all photon qualities
	} else if (trainConfig == 14) { 
		eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00700009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , all photon qualities, min R = 35 cm
	} else if (trainConfig == 15) { 
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; //standard cut LHC11h pp 2.76TeV 
	} else if (trainConfig == 16) { 
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; //standard cut pp 8 TeV		
	} else if (trainConfig == 17){
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500000000"; //New standard cut for eta analysis 8 TeV
	} else if (trainConfig == 18){
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152101500000000"; //New standard cut for eta analysis 8 TeV
	} else if (trainConfig == 19){
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152101500000002"; //New standard cut for eta analysis 8 TeV
	} else if (trainConfig == 20){
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152101500000022"; //New standard cut for eta analysis 8 TeV	
	} else if (trainConfig == 21){
		eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250404000"; mesonCutArray[0] = "0152101500000000"; //New standard cut for eta analysis 8 TeV
	} else {
		Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
	TList *ConvCutList = new TList();
	TList *MesonCutList = new TList();

	TList *HeaderList = new TList();
	TObjString *Header2 = new TObjString("BOX");
	HeaderList->Add(Header2);

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
		if (trainConfig == 21){
		  analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);
		}
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

	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
