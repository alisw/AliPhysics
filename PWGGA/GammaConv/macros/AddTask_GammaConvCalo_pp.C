void AddTask_GammaConvCalo_pp(  Int_t trainConfig = 1,  //change different set of cuts
								Bool_t isMC = kFALSE, //run MC
								Int_t enableQAMesonTask = 1, //enable QA in AliAnalysisTaskGammaConvV1
								Int_t enableQAPhotonTask = 1, // enable additional QA task
								TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                                TString cutnumberAODBranch = "0000000060084001001500000",
                                Bool_t enableExtendedMatching = kFALSE //enable or disable extended matching histograms for conversion electrons <-> cluster
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
	gSystem->Load("libCDB.so");
	gSystem->Load("libSTEER.so");
	gSystem->Load("libSTEERBase.so");
	gSystem->Load("libTENDER.so");
	gSystem->Load("libTENDERSupplies.so");
	gSystem->Load("libPWGflowBase.so");
	gSystem->Load("libPWGflowTasks.so");
	gSystem->Load("libPWGGAGammaConv.so");
	
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
	TString cutnumberPhoton = "060084001001500000000";
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
	AliAnalysisTaskGammaConvCalo *task=NULL;
	task= new AliAnalysisTaskGammaConvCalo(Form("GammaConvCalo_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 2;
	if (trainConfig==2 || trainConfig==3 || trainConfig==4 || trainConfig==5 || trainConfig==6 || trainConfig==7 ){ numberOfCuts =5;}
	if (trainConfig==8 ){ numberOfCuts =4;}
	if (trainConfig==10 || trainConfig==31 || trainConfig==32 ){ numberOfCuts =6;}

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
		eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "002093270008250400000"; clusterCutArray[0] = "10000042032030000"; mesonCutArray[0] = "01631031000000"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "0005111"; photonCutArray[ 1] = "002093270008250400000"; clusterCutArray[1] = "10000042032030000"; mesonCutArray[1] = "01631031000000"; // 400 MeV cluster min energy
	} else if (trainConfig == 2){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "002093270008250400000"; clusterCutArray[0] = "10000042012030000"; mesonCutArray[0] = "01631031000000"; //0.2 GeV/c
		eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "002093270008250400000"; clusterCutArray[1] = "10000042022030000"; mesonCutArray[1] = "01631031000000"; //0.3 GeV/c
		eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "002093270008250400000"; clusterCutArray[2] = "10000042032030000"; mesonCutArray[2] = "01631031000000"; //0.4 GeV/c default
		eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "002093270008250400000"; clusterCutArray[3] = "10000042042030000"; mesonCutArray[3] = "01631031000000"; //0.5 GeV/c
		eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "002093270008250400000"; clusterCutArray[4] = "10000042052030000"; mesonCutArray[4] = "01631031000000"; //0.6 GeV/c
	} else if (trainConfig == 3){ //EMCAL minNCells variation
		eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "002093270008250400000"; clusterCutArray[0] = "10000042031030000"; mesonCutArray[0] = "01631031000000"; //n cells >= 1
		eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "002093270008250400000"; clusterCutArray[1] = "10000042033030000"; mesonCutArray[1] = "01631031000000"; //n cells >= 3
		eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "002093270008250400000"; clusterCutArray[2] = "10000042032000000"; mesonCutArray[2] = "01631031000000"; //no M02 cut
		eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "002093270008250400000"; clusterCutArray[3] = "10031042032030000"; mesonCutArray[3] = "01631031000000"; //only modules with TRD infront
		eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "002093270008250400000"; clusterCutArray[4] = "10012042032030000"; mesonCutArray[4] = "01631031000000"; //no modules with TRD infront		
	} else if (trainConfig == 4){ // EMCAL track matching variations 
		eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "002093270008250400000"; clusterCutArray[0] = "10000041032030000"; mesonCutArray[0] = "01631031000000"; // 
		eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "002093270008250400000"; clusterCutArray[1] = "10000042032030000"; mesonCutArray[1] = "01631031000000"; // 
		eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "002093270008250400000"; clusterCutArray[2] = "10000043032030000"; mesonCutArray[2] = "01631031000000"; // 
		eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "002093270008250400000"; clusterCutArray[3] = "10000044032030000"; mesonCutArray[3] = "01631031000000"; // 
		eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "002093270008250400000"; clusterCutArray[4] = "10000045032030000"; mesonCutArray[4] = "01631031000000"; // 
	} else if (trainConfig == 5){ // PCM variations
		eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "002092270008250400000"; clusterCutArray[0] = "10000042032030000"; mesonCutArray[0] = "01631031000000"; // dEdx e -3, 5
		eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "002091270008250400000"; clusterCutArray[1] = "10000042032030000"; mesonCutArray[1] = "01631031000000"; // dEdx e -5, 5
		eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "002093570008250400000"; clusterCutArray[2] = "10000042032030000"; mesonCutArray[2] = "01631031000000"; // dEdx pi 2
		eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "002093170008250400000"; clusterCutArray[3] = "10000042032030000"; mesonCutArray[3] = "01631031000000"; // dEdx pi 0
		eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "002093873008250400000"; clusterCutArray[4] = "10000042032030000"; mesonCutArray[4] = "01631031000000"; // dEdx pi 2 high 1 (> 3.5 GeV)
	} else if (trainConfig == 6){ // PCM variations
		eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "002093270009250400000"; clusterCutArray[0] = "10000042032030000"; mesonCutArray[0] = "01631031000000"; // qt 2D 0.03
		eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "002093270003250400000"; clusterCutArray[1] = "10000042032030000"; mesonCutArray[1] = "01631031000000"; // qt 1D 0.05
		eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "002093270002250400000"; clusterCutArray[2] = "10000042032030000"; mesonCutArray[2] = "01631031000000"; // qt 1D 0.07
		eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "002493270008250400000"; clusterCutArray[3] = "10000042032030000"; mesonCutArray[3] = "01631031000000"; // single pt > 0.075
		eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "002193270008250400000"; clusterCutArray[4] = "10000042032030000"; mesonCutArray[4] = "01631031000000"; // single pt > 0.1
	} else if (trainConfig == 7){ // PCM variations
		eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "002093270008850400000"; clusterCutArray[0] = "10000042032030000"; mesonCutArray[0] = "01631031000000"; // 2D psi pair chi2 var
		eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "002093270008260400000"; clusterCutArray[1] = "10000042032030000"; mesonCutArray[1] = "01631031000000"; // 2D psi pair chi2 var
		eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "002093270008860400000"; clusterCutArray[2] = "10000042032030000"; mesonCutArray[2] = "01631031000000"; // 2D psi pair chi2 var
		eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "002093270008280400000"; clusterCutArray[3] = "10000042032030000"; mesonCutArray[3] = "01631031000000"; // 2D psi pair chi2 var
		eventCutArray[ 4] = "0000311"; photonCutArray[ 4] = "002093270008880400000"; clusterCutArray[4] = "10000042032030000"; mesonCutArray[4] = "01631031000000"; // 2D psi pair chi2 var
	} else if (trainConfig == 8){ // PCM variations
		eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "002063270008250400000"; clusterCutArray[0] = "10000042032030000"; mesonCutArray[0] = "01631031000000"; // min TPC cl > 0.7
		eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "002083270008250400000"; clusterCutArray[1] = "10000042032030000"; mesonCutArray[1] = "01631031000000"; // min TPC cl > 0.35
		eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "002093270008250400000"; clusterCutArray[2] = "10000042032030000"; mesonCutArray[2] = "01631061000000"; // alpha < 0.8
		eventCutArray[ 3] = "0000311"; photonCutArray[ 3] = "002093270008250400000"; clusterCutArray[3] = "10000042032030000"; mesonCutArray[3] = "01631051000000"; // alpha < 0.75

		// LHC13g	
	} else if (trainConfig == 10){  // EMCAL clusters, EMCEGA triggers, track matching 0.035
		eventCutArray[ 0] = "0008311"; photonCutArray[ 0] = "002093270008250400000"; clusterCutArray[0] = "10000042032030000"; mesonCutArray[0] = "01631031000000"; // EMCEG1, 
		eventCutArray[ 1] = "0008511"; photonCutArray[ 1] = "002093270008250400000"; clusterCutArray[1] = "10000042032030000"; mesonCutArray[1] = "01631031000000"; // EMCEG2, 
		eventCutArray[ 2] = "0009311"; photonCutArray[ 2] = "002093270008250400000"; clusterCutArray[2] = "10000042032030000"; mesonCutArray[2] = "01631031000000"; // EMCEJ1, 
		eventCutArray[ 3] = "0009511"; photonCutArray[ 3] = "002093270008250400000"; clusterCutArray[3] = "10000042032030000"; mesonCutArray[3] = "01631031000000"; // EMCEJ2, 
		eventCutArray[ 4] = "0000011"; photonCutArray[ 4] = "002093270008250400000"; clusterCutArray[4] = "10000042032030000"; mesonCutArray[4] = "01631031000000"; // INT7
		eventCutArray[ 5] = "0005211"; photonCutArray[ 5] = "002093270008250400000"; clusterCutArray[5] = "10000042032030000"; mesonCutArray[5] = "01631031000000"; // EMC7
		
	// ************************************* PHOS cuts ****************************************************
	// LHC11a	
	} else if (trainConfig == 31) { //PHOS clusters
		eventCutArray[ 0] = "0000311"; photonCutArray[ 0] = "002093270008250400000"; clusterCutArray[0] = "20000041033200000"; mesonCutArray[0] = "01631031000000"; 
		eventCutArray[ 1] = "0000311"; photonCutArray[ 1] = "002093270008250400000"; clusterCutArray[1] = "20000042033200000"; mesonCutArray[1] = "01631031000000"; 
		eventCutArray[ 2] = "0000311"; photonCutArray[ 2] = "002093270008250400000"; clusterCutArray[2] = "20000043033200000"; mesonCutArray[2] = "01631031000000"; 
		eventCutArray[ 3] = "0006111"; photonCutArray[ 3] = "002093270008250400000"; clusterCutArray[3] = "20000041033200000"; mesonCutArray[3] = "01631031000000"; 
		eventCutArray[ 4] = "0006111"; photonCutArray[ 4] = "002093270008250400000"; clusterCutArray[4] = "20000042033200000"; mesonCutArray[4] = "01631031000000"; 
		eventCutArray[ 5] = "0006111"; photonCutArray[ 5] = "002093270008250400000"; clusterCutArray[5] = "20000043033200000"; mesonCutArray[5] = "01631031000000"; 
	// LHC13g & LHC12x
	} else if (trainConfig == 32) { //PHOS clusters
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002093270008250400000"; clusterCutArray[0] = "20000041033200000"; mesonCutArray[0] = "01631031000000"; 
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002093270008250400000"; clusterCutArray[1] = "20000042033200000"; mesonCutArray[1] = "01631031000000"; 
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002093270008250400000"; clusterCutArray[2] = "20000043033200000"; mesonCutArray[2] = "01631031000000"; 
		eventCutArray[ 3] = "0006211"; photonCutArray[ 3] = "002093270008250400000"; clusterCutArray[3] = "20000041033200000"; mesonCutArray[3] = "01631031000000"; 
		eventCutArray[ 4] = "0006211"; photonCutArray[ 4] = "002093270008250400000"; clusterCutArray[4] = "20000042033200000"; mesonCutArray[4] = "01631031000000"; 
		eventCutArray[ 5] = "0006211"; photonCutArray[ 5] = "002093270008250400000"; clusterCutArray[5] = "20000043033200000"; mesonCutArray[5] = "01631031000000"; 
	} else {
		Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
	TList *ConvCutList = new TList();
	TList *ClusterCutList = new TList();
	TList *MesonCutList = new TList();

// 	TList *HeaderList = new TList();
// 	if (doWeightingPart==1) {
// 		TObjString *Header1 = new TObjString("pi0_1");
// 		HeaderList->Add(Header1);
// 	}
// 	if (doWeightingPart==2){
// 		TObjString *Header3 = new TObjString("eta_2");
// 		HeaderList->Add(Header3);
// 	}
// 	if (doWeightingPart==3) {
// 		TObjString *Header1 = new TObjString("pi0_1");
// 		HeaderList->Add(Header1);
// 		TObjString *Header3 = new TObjString("eta_2");
// 		HeaderList->Add(Header3);
// 	}

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
        analysisClusterCuts[i]->SetExtendedMatching(enableExtendedMatching);
		analysisClusterCuts[i]->SetFillCutHistograms("");
		
		analysisMesonCuts[i] = new AliConversionMesonCuts();
		analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
		MesonCutList->Add(analysisMesonCuts[i]);
		analysisMesonCuts[i]->SetFillCutHistograms("");
// 		analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
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
