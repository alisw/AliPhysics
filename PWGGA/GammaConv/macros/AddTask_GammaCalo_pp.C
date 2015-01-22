void AddTask_GammaCalo_pp(  Int_t trainConfig = 1,  //change different set of cuts
							Bool_t isMC   = kFALSE, //run MC
							Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvV1
							Int_t enableQAClusterTask = 0, // enable additional QA task
							TString fileNameInputForWeighting = "MCSpectraInput.root", 	// path to file for weigting input
                            TString cutnumberAODBranch = "0000000060084001001500000",
							TString periodname = "LHC12f1x", 							// period name
							Bool_t doWeighting = kFALSE									// enables weighting
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
	TString cutnumberPhoton = "060000084001001500000000";
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
	AliAnalysisTaskGammaCalo *task=NULL;
	task= new AliAnalysisTaskGammaCalo(Form("GammaCalo_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 2;
	if (trainConfig == 2 || trainConfig == 3) numberOfCuts = 5;
	if (trainConfig == 5) numberOfCuts = 6;
	if (trainConfig == 31) numberOfCuts = 4;
	if (trainConfig == 32) numberOfCuts = 1;
	
	TString *eventCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// cluster cuts
	// 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
	// 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
	
	// ************************************* EMCAL cuts ****************************************************
	// LHC11a
	if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "10000040032030000"; mesonCutArray[0] = "01631031000000"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "0005111"; clusterCutArray[1] = "10000040032030000"; mesonCutArray[1] = "01631031000000"; // 400 MeV cluster min energy
	} else if (trainConfig == 2){ //EMCAL minEnergy variation
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "10000040012030000"; mesonCutArray[0] = "01631031000000"; //0.2 GeV/c
		eventCutArray[ 1] = "0000311"; clusterCutArray[1] = "10000040022030000"; mesonCutArray[1] = "01631031000000"; //0.3 GeV/c
		eventCutArray[ 2] = "0000311"; clusterCutArray[2] = "10000040032030000"; mesonCutArray[2] = "01631031000000"; //0.4 GeV/c default
		eventCutArray[ 3] = "0000311"; clusterCutArray[3] = "10000040042030000"; mesonCutArray[3] = "01631031000000"; //0.5 GeV/c
		eventCutArray[ 4] = "0000311"; clusterCutArray[4] = "10000040052030000"; mesonCutArray[4] = "01631031000000"; //0.6 GeV/c
	} else if (trainConfig == 3){ //EMCAL minNCells variation
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "10000040031030000"; mesonCutArray[0] = "01631031000000"; //n cells >= 1
		eventCutArray[ 1] = "0000311"; clusterCutArray[1] = "10000040033030000"; mesonCutArray[1] = "01631031000000"; //n cells >= 3
		eventCutArray[ 2] = "0000311"; clusterCutArray[2] = "10000040032000000"; mesonCutArray[2] = "01631031000000"; //no M02 cut
		eventCutArray[ 3] = "0000311"; clusterCutArray[3] = "10031040032030000"; mesonCutArray[3] = "01631031000000"; //only modules with TRD infront
		eventCutArray[ 4] = "0000311"; clusterCutArray[4] = "10012040032030000"; mesonCutArray[4] = "01631031000000"; //no modules with TRD infront		
	// LHC13g	
	} else if (trainConfig == 5){  // EMCAL clusters, EMCEGA triggers
		eventCutArray[ 0] = "0008311"; clusterCutArray[0] = "10000040032030000"; mesonCutArray[0] = "01631031000000"; // EMCEG1, 
		eventCutArray[ 1] = "0008511"; clusterCutArray[1] = "10000040032030000"; mesonCutArray[1] = "01631031000000"; // EMCEG2, 
		eventCutArray[ 2] = "0009311"; clusterCutArray[2] = "10000040032030000"; mesonCutArray[2] = "01631031000000"; // EMCEJ1, 
		eventCutArray[ 3] = "0009511"; clusterCutArray[3] = "10000040032030000"; mesonCutArray[3] = "01631031000000"; // EMCEJ2, 
		eventCutArray[ 4] = "0000011"; clusterCutArray[4] = "10000040032030000"; mesonCutArray[4] = "01631031000000"; // INT7
		eventCutArray[ 5] = "0005211"; clusterCutArray[5] = "10000040032030000"; mesonCutArray[5] = "01631031000000"; // EMC7
	} else if (trainConfig == 12){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0) without and with added signals
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "10000040032030000"; mesonCutArray[0] = "01631031000000"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "0000312"; clusterCutArray[1] = "10000040032030000"; mesonCutArray[1] = "01631031000000"; // 400 MeV cluster min energy

	// ************************************* PHOS cuts ****************************************************
	} else if (trainConfig == 31) { //PHOS clusters
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "20000040033200000"; mesonCutArray[0] = "01631031000000"; //pp LHC11a with SDD, PHOS
		eventCutArray[ 1] = "0000011"; clusterCutArray[1] = "20000040033200000"; mesonCutArray[1] = "01631031000000"; //pp LHC13g default MB
		eventCutArray[ 2] = "0006111"; clusterCutArray[2] = "20000040033200000"; mesonCutArray[2] = "01631031000000"; //pp LHC11a PHI1
		eventCutArray[ 3] = "0006211"; clusterCutArray[3] = "20000040033200000"; mesonCutArray[3] = "01631031000000"; //pp LHC11a PHI7
	} else if (trainConfig == 32){ // Validation PHOS
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "20000040033200000"; mesonCutArray[0] = "01630031009000"; 		
	} else if (trainConfig == 33){ // PHOS clusters, without and with added signals
		eventCutArray[ 0] = "0000311"; clusterCutArray[0] = "20000040033200000"; mesonCutArray[0] = "01630031009000"; 		
		eventCutArray[ 1] = "0000312"; clusterCutArray[1] = "20000040033200000"; mesonCutArray[1] = "01630031009000"; 		
	} else {
		Error(Form("GammaCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
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
	task->SetCaloCutList(numberOfCuts,ClusterCutList);
	task->SetMesonCutList(numberOfCuts,MesonCutList);
	task->SetDoMesonAnalysis(kTRUE);
	task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
	task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA

	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer(Form("GammaCalo_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaCalo_%i.root",trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
