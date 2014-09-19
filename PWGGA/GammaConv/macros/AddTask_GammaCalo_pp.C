void AddTask_GammaCalo_pp(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC
                              Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvV1
                              Int_t enableQAClusterTask = 0, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                              TString cutnumberAODBranch = "0000000060084000001500000" // cutnumber for AOD branch
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
	AliAnalysisTaskGammaCalo *task=NULL;
	task= new AliAnalysisTaskGammaCalo(Form("GammaCalo_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 4;
	
	TString *eventCutArray = new TString[numberOfCuts];
	TString *clusterCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// cluster cuts
	// 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
	// 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
	if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0,1), kEMC1 (2,3)
		eventCutArray[ 0] = "0002011"; clusterCutArray[0] = "10000040022030000"; mesonCutArray[0] = "01631031009000"; // 100 MeV cluster min energy
		eventCutArray[ 1] = "0002011"; clusterCutArray[1] = "10000040052030000"; mesonCutArray[1] = "01631031009000"; // 300 MeV cluster min energy
		eventCutArray[ 2] = "0005111"; clusterCutArray[2] = "10000040022030000"; mesonCutArray[2] = "01631031009000"; // 100 MeV cluster min energy
		eventCutArray[ 3] = "0005111"; clusterCutArray[3] = "10000040052030000"; mesonCutArray[3] = "01631031009000"; // 300 MeV cluster min energy
	} else if (trainConfig == 2){  // EMCAL clusters, EMCEGA triggers
		eventCutArray[ 0] = "0008311"; clusterCutArray[0] = "10000040022030000"; mesonCutArray[0] = "01631031009000"; // EMCEG1, 100 MeV cluster min energy
		eventCutArray[ 1] = "0008311"; clusterCutArray[1] = "10000040052030000"; mesonCutArray[1] = "01631031009000"; // EMCEG1, 300 MeV cluster min energy
		eventCutArray[ 2] = "0008511"; clusterCutArray[2] = "10000040022030000"; mesonCutArray[2] = "01631031009000"; // EMCEG2, 100 MeV cluster min energy
		eventCutArray[ 3] = "0008511"; clusterCutArray[3] = "10000040052030000"; mesonCutArray[3] = "01631031009000"; // EMCEG2, 300 MeV cluster min energy
	} else if (trainConfig == 3){  // EMCAL clusters, EMCEJE triggers
		eventCutArray[ 0] = "0009311"; clusterCutArray[0] = "10000040022030000"; mesonCutArray[0] = "01631031009000"; // EMCEJ1, 100 MeV cluster min energy
		eventCutArray[ 1] = "0009311"; clusterCutArray[1] = "10000040052030000"; mesonCutArray[1] = "01631031009000"; // EMCEJ1, 300 MeV cluster min energy
		eventCutArray[ 2] = "0009511"; clusterCutArray[2] = "10000040022030000"; mesonCutArray[2] = "01631031009000"; // EMCEJ2, 100 MeV cluster min energy
		eventCutArray[ 3] = "0009511"; clusterCutArray[3] = "10000040052030000"; mesonCutArray[3] = "01631031009000"; // EMCEJ2, 300 MeV cluster min energy
	} else if (trainConfig == 4){ // EMCAL clusters 2.76 TeV LHC13g, kINT7 (0,1), kEMC7 (2,3), track 
		eventCutArray[ 0] = "0000011"; clusterCutArray[0] = "10000040022030000"; mesonCutArray[0] = "01631031009000"; // 100 MeV cluster min energy
		eventCutArray[ 1] = "0000011"; clusterCutArray[1] = "10000040052030000"; mesonCutArray[1] = "01631031009000"; // 300 MeV cluster min energy
		eventCutArray[ 2] = "0005211"; clusterCutArray[2] = "10000040022030000"; mesonCutArray[2] = "01631031009000"; // 100 MeV cluster min energy
		eventCutArray[ 3] = "0005211"; clusterCutArray[3] = "10000040052030000"; mesonCutArray[3] = "01631031009000"; // 300 MeV cluster min energy
	} else if (trainConfig == 5){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0,1), kEMC1 (2,3)
		eventCutArray[ 0] = "0002011"; clusterCutArray[0] = "10000040062030000"; mesonCutArray[0] = "01631031009000"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "0002011"; clusterCutArray[1] = "10000040062031000"; mesonCutArray[1] = "01631031009000"; // 400 MeV cluster min energy, min M20 > 0.02
		eventCutArray[ 2] = "0005111"; clusterCutArray[2] = "10000040062030000"; mesonCutArray[2] = "01631031009000"; // 500 MeV cluster min energy
		eventCutArray[ 3] = "0005111"; clusterCutArray[3] = "10000040062031000"; mesonCutArray[3] = "01631031009000"; // 500 MeV cluster min energy, min M20 > 0.02
	} else if (trainConfig == 6){  // EMCAL clusters, EMCEGA triggers
		eventCutArray[ 0] = "0008311"; clusterCutArray[0] = "10000040062030000"; mesonCutArray[0] = "01631031009000"; // EMCEG1, 400 MeV cluster min energy
		eventCutArray[ 1] = "0008311"; clusterCutArray[1] = "10000040062031000"; mesonCutArray[1] = "01631031009000"; // EMCEG1, 400 MeV cluster min energy, min M20 > 0.02
		eventCutArray[ 2] = "0008511"; clusterCutArray[2] = "10000040062030000"; mesonCutArray[2] = "01631031009000"; // EMCEG2, 400 MeV cluster min energy
		eventCutArray[ 3] = "0008511"; clusterCutArray[3] = "10000040062031000"; mesonCutArray[3] = "01631031009000"; // EMCEG2, 400 MeV cluster min energy, min M20 > 0.02
	} else if (trainConfig == 7){  // EMCAL clusters, EMCEJE triggers
		eventCutArray[ 0] = "0009311"; clusterCutArray[0] = "10000040062030000"; mesonCutArray[0] = "01631031009000"; // EMCEJ1, 400 MeV cluster min energy
		eventCutArray[ 1] = "0009311"; clusterCutArray[1] = "10000040062031000"; mesonCutArray[1] = "01631031009000"; // EMCEJ1, 400 MeV cluster min energy, min M20 > 0.02
		eventCutArray[ 2] = "0009511"; clusterCutArray[2] = "10000040062030000"; mesonCutArray[2] = "01631031009000"; // EMCEJ2, 400 MeV cluster min energy
		eventCutArray[ 3] = "0009511"; clusterCutArray[3] = "10000040062031000"; mesonCutArray[3] = "01631031009000"; // EMCEJ2, 400 MeV cluster min energy, min M20 > 0.02
	} else if (trainConfig == 8){ // EMCAL clusters 2.76 TeV LHC13g, kINT7 (0,1), kEMC7 (2,3), track 
		eventCutArray[ 0] = "0000011"; clusterCutArray[0] = "10000040062030000"; mesonCutArray[0] = "01631031009000"; // 400 MeV cluster min energy
		eventCutArray[ 1] = "0000011"; clusterCutArray[1] = "10000040062031000"; mesonCutArray[1] = "01631031009000"; // 400 MeV cluster min energy, min M20 > 0.02
		eventCutArray[ 2] = "0005211"; clusterCutArray[2] = "10000040062030000"; mesonCutArray[2] = "01631031009000"; // 400 MeV cluster min energy
		eventCutArray[ 3] = "0005211"; clusterCutArray[3] = "10000040062031000"; mesonCutArray[3] = "01631031009000"; // 400 MeV cluster min energy, min M20 > 0.02
	} else if (trainConfig == 31) { //PHOS clusters
		eventCutArray[ 0] = "0002011"; clusterCutArray[0] = "20000030022000000"; mesonCutArray[0] = "01631031009000"; //pp LHC11a with SDD, PHOS
		eventCutArray[ 1] = "0000011"; clusterCutArray[1] = "20000030022000000"; mesonCutArray[1] = "01631031009000"; //pp LHC13g default MB
		eventCutArray[ 2] = "0006111"; clusterCutArray[2] = "20000030022000000"; mesonCutArray[2] = "01631031009000"; //pp LHC11a PHI1
		eventCutArray[ 3] = "0006211"; clusterCutArray[3] = "20000030022000000"; mesonCutArray[3] = "01631031009000"; //pp LHC11a PHI7
	} else {
		Error(Form("GammaCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
	TList *ClusterCutList = new TList();
	TList *MesonCutList = new TList();


	EventCutList->SetOwner(kTRUE);
	AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
	ClusterCutList->SetOwner(kTRUE);
	AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

	for(Int_t i = 0; i<numberOfCuts; i++){
		analysisEventCuts[i] = new AliConvEventCuts();   
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
