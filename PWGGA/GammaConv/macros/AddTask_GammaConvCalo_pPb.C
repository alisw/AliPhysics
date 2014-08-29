void AddTask_GammaConvCalo_pPb(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC
                              Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvV1
                              Int_t enableQAPhotonTask = 0, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                              Int_t doWeightingPart = 0,  //enable Weighting
                              TString generatorName = "DPMJET",
                              TString cutnumberAODBranch = "8000000060084000001500000" // cutnumber for AOD branch
							) {

	// ================= Load Librariers =================================
	gSystem->Load("libCore.so");  
	gSystem->Load("libTree.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libVMC.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libMinuit");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");  
	gSystem->Load("libPWGGAGammaConv.so");
	gSystem->Load("libCDB.so");
	gSystem->Load("libSTEER.so");
	gSystem->Load("libSTEERBase.so");
	gSystem->Load("libTENDER.so");
	gSystem->Load("libTENDERSupplies.so");
		
	Int_t isHeavyIon = 2;
	
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
	TString cutnumberPhoton = "060084001001500000000";
	TString cutnumberEvent = "8000000";
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
	AliAnalysisTaskGammaConvCalo *task=NULL;
	task= new AliAnalysisTaskGammaConvCalo(Form("GammaConvCalo_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 2;
	if (trainConfig==3 || trainConfig==5  ){ numberOfCuts =6;}
	if (trainConfig==6   ){ numberOfCuts =5;}
	
	TString *eventCutArray = new TString[numberOfCuts];
	TString *photonCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// cluster cuts
	// 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
	// 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
	
	if (trainConfig == 1){ 
		eventCutArray[ 0] = "8000001"; photonCutArray[ 0] = "002092970028250400000"; clusterCutArray[0] = "10000040022030000"; mesonCutArray[0] = "01525065000000"; //standart cut, kINT7 // EMCAL clusters
		eventCutArray[ 1] = "8005201"; photonCutArray[ 1] = "002092970028250400000"; clusterCutArray[1] = "10000040022030000"; mesonCutArray[1] = "01525065000000"; //standard cut, kEMC7 // EMCAL clusters
	} else if (trainConfig == 2) {	
		eventCutArray[ 0] = "8000001"; photonCutArray[ 0] = "002092970028250400000"; clusterCutArray[0] = "20000030022000000"; mesonCutArray[0] = "01525065000000"; //standart cut, kINT7 // PHOS clusters
		eventCutArray[ 1] = "8006201"; photonCutArray[ 1] = "002092970028250400000"; clusterCutArray[1] = "20000030022000000"; mesonCutArray[1] = "01525065000000"; //standard cut, kPHI7	// PHOS clusters	
	} else if (trainConfig == 3){ 
		eventCutArray[ 0] = "8008101"; photonCutArray[ 0] = "002092970028250400000"; clusterCutArray[0] = "10000040022030000"; mesonCutArray[0] = "01525065000000"; //standart cut, kEMCEGA based on INT7 // EMCAL clusters
		eventCutArray[ 1] = "8008301"; photonCutArray[ 1] = "002092970028250400000"; clusterCutArray[1] = "10000040022030000"; mesonCutArray[1] = "01525065000000"; //standard cut, kEMCEG1 based on INT7 // EMCAL clusters
		eventCutArray[ 2] = "8008501"; photonCutArray[ 2] = "002092970028250400000"; clusterCutArray[2] = "10000040022030000"; mesonCutArray[2] = "01525065000000"; //standard cut, kEMCEG2 based on INT7 // EMCAL clusters
		eventCutArray[ 3] = "8009101"; photonCutArray[ 3] = "002092970028250400000"; clusterCutArray[3] = "10000040022030000"; mesonCutArray[3] = "01525065000000"; //standard cut, kEMCEJE based on INT7 // EMCAL clusters
		eventCutArray[ 4] = "8009301"; photonCutArray[ 4] = "002092970028250400000"; clusterCutArray[4] = "10000040022030000"; mesonCutArray[4] = "01525065000000"; //standard cut, kEMCEJ1 based on INT7 // EMCAL clusters
		eventCutArray[ 5] = "8009501"; photonCutArray[ 5] = "002092970028250400000"; clusterCutArray[5] = "10000040022030000"; mesonCutArray[5] = "01525065000000"; //standard cut, kEMCEG2 based on INT7 // EMCAL clusters
	} else if (trainConfig == 4){ 
		eventCutArray[ 0] = "8000001"; photonCutArray[ 0] = "002092970028250400000"; clusterCutArray[0] = "10000040052030000"; mesonCutArray[0] = "01525065000000"; //standart cut, kINT7 // EMCAL clusters
		eventCutArray[ 1] = "8005201"; photonCutArray[ 1] = "002092970028250400000"; clusterCutArray[1] = "10000040052030000"; mesonCutArray[1] = "01525065000000"; //standard cut, kEMC7 // EMCAL clusters	
	} else if (trainConfig == 5){ 
		eventCutArray[ 0] = "8008101"; photonCutArray[ 0] = "002092970028250400000"; clusterCutArray[0] = "10000040052030000"; mesonCutArray[0] = "01525065000000"; //standart cut, kEMCEGA based on INT7 // EMCAL clusters
		eventCutArray[ 1] = "8008301"; photonCutArray[ 1] = "002092970028250400000"; clusterCutArray[1] = "10000040052030000"; mesonCutArray[1] = "01525065000000"; //standard cut, kEMCEG1 based on INT7 // EMCAL clusters
		eventCutArray[ 2] = "8008501"; photonCutArray[ 2] = "002092970028250400000"; clusterCutArray[2] = "10000040052030000"; mesonCutArray[2] = "01525065000000"; //standard cut, kEMCEG2 based on INT7 // EMCAL clusters
		eventCutArray[ 3] = "8009101"; photonCutArray[ 3] = "002092970028250400000"; clusterCutArray[3] = "10000040052030000"; mesonCutArray[3] = "01525065000000"; //standard cut, kEMCEJE based on INT7 // EMCAL clusters
		eventCutArray[ 4] = "8009301"; photonCutArray[ 4] = "002092970028250400000"; clusterCutArray[4] = "10000040052030000"; mesonCutArray[4] = "01525065000000"; //standard cut, kEMCEJ1 based on INT7 // EMCAL clusters
		eventCutArray[ 5] = "8009501"; photonCutArray[ 5] = "002092970028250400000"; clusterCutArray[5] = "10000040052030000"; mesonCutArray[5] = "01525065000000"; //standard cut, kEMCEG2 based on INT7 // EMCAL clusters		
	} else if (trainConfig == 6){ 
		eventCutArray[ 0] = "8000001"; photonCutArray[ 0] = "002092970028250400000"; clusterCutArray[0] = "10000043052030000"; mesonCutArray[0] = "01525065000000"; //standart cut, kINT7 // EMCAL clusters
		eventCutArray[ 1] = "8000001"; photonCutArray[ 1] = "002092970028250400000"; clusterCutArray[1] = "10000044052030000"; mesonCutArray[1] = "01525065000000"; //standard cut, kEMC7 // EMCAL clusters	
		eventCutArray[ 2] = "8000001"; photonCutArray[ 2] = "002092970028250400000"; clusterCutArray[2] = "10000045052030000"; mesonCutArray[2] = "01525065000000"; //standard cut, kEMC7 // EMCAL clusters	
		eventCutArray[ 3] = "8000001"; photonCutArray[ 3] = "002092970028250400000"; clusterCutArray[3] = "10000046052030000"; mesonCutArray[3] = "01525065000000"; //standard cut, kEMC7 // EMCAL clusters	
		eventCutArray[ 4] = "8000001"; photonCutArray[ 4] = "002092970028250400000"; clusterCutArray[4] = "10000047052030000"; mesonCutArray[4] = "01525065000000"; //standard cut, kEMC7 // EMCAL clusters	
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
