void AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_MixedMode_pPb(    
										Int_t trainConfig = 1,
										Bool_t isMC       = kFALSE, //run MC 
										Bool_t enableQAMesonTask = kTRUE, //enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero
										TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
										Bool_t doWeighting = kFALSE,  //enable Weighting
										TString generatorName = "HIJING",				
										TString cutnumberAODBranch = "000000006008400001001500000"
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
	Int_t neutralPionMode = 1;
	
	// ================== GetAnalysisManager ===============================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error(Form("AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_MixedMode_pPb_%i",trainConfig), "No analysis manager found.");
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
	TString cutnumberPhoton = "06000008400100001500000000";
	TString cutnumberEvent = "80000003";
	TString PionCuts      = "000000200";            //Electron Cuts
	

	
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

		AliLog::SetGlobalLogLevel(AliLog::kInfo);

		//connect input V0Reader
		mgr->AddTask(fV0ReaderV1);
		mgr->ConnectInput(fV0ReaderV1,0,cinput);
	}

	//================================================
	//========= Add Electron Selector ================


	if( !(AliPrimaryPionSelector*)mgr->GetTask("PionSelector") ){

		AliPrimaryPionSelector *fPionSelector = new AliPrimaryPionSelector("PionSelector");
		// Set AnalysisCut Number

		AliPrimaryPionCuts *fPionCuts=0;
		if( PionCuts!=""){
			fPionCuts= new AliPrimaryPionCuts(PionCuts.Data(),PionCuts.Data());
			if(fPionCuts->InitializeCutsFromCutString(PionCuts.Data())){
				fPionSelector->SetPrimaryPionCuts(fPionCuts);
				fPionCuts->SetFillCutHistograms("",kTRUE);

			}
		}

		fPionSelector->Init();
		mgr->AddTask(fPionSelector);
		
		AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();

		//connect input V0Reader
		mgr->ConnectInput (fPionSelector,0,cinput1);

	}

	
	
	AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero *task=NULL;

	task= new AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero(Form("GammaConvNeutralMesonPiPlPiMiPiZero_%i_%i",neutralPionMode, trainConfig));

	task->SetIsHeavyIon(2);
	task->SetIsMC(isMC);

	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 1;

	TString *eventCutArray 			= new TString[numberOfCuts];
	TString *ClusterCutarray   		= new TString[numberOfCuts];
	TString *ConvCutarray   		= new TString[numberOfCuts];
	TString *PionCutarray    		= new TString[numberOfCuts];
	TString *NeutralPionCutarray   	= new TString[numberOfCuts];
	TString *MesonCutarray   		= new TString[numberOfCuts];
	
	Bool_t doEtaShiftIndCuts = kFALSE;
	TString stringShift = "";

	// Shifting in pPb direction

	doEtaShiftIndCuts = kTRUE;
	stringShift = "pPb";

	// EMCAL mode
	if( trainConfig == 1 ) {
		// everything open
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "1111100040022030000"; PionCutarray[0] = "000010400"; NeutralPionCutarray[0] = "01035030000000"; MesonCutarray[0] = "0103503000000000"; 
	} else if( trainConfig == 2 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "1111100040022030000"; PionCutarray[0] = "002010700"; NeutralPionCutarray[0] = "01035030000000"; MesonCutarray[0] = "0103503000000000"; 
	} else if( trainConfig == 3 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "1111100040022030000"; PionCutarray[0] = "002013700"; NeutralPionCutarray[0] = "01035030000000"; MesonCutarray[0] = "0103503000000000"; 
	} else if( trainConfig == 4 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "1111100040022030000"; PionCutarray[0] = "002016700"; NeutralPionCutarray[0] = "01035030000000"; MesonCutarray[0] = "0103503000000000"; 
	} else if( trainConfig == 5 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
		// closing neural pion cuts, 0.1 < M_gamma,gamma < 0.145
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "1111100040022030000"; PionCutarray[0] = "002016700"; NeutralPionCutarray[0] = "01035031000000"; MesonCutarray[0] = "0103503000000000"; 	
	} else if( trainConfig == 6 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
		// closing neural pion cuts, 0.11 < M_gamma,gamma < 0.145
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "1111100040022030000"; PionCutarray[0] = "002016700"; NeutralPionCutarray[0] = "01035032000000"; MesonCutarray[0] = "0103503000000000"; 	
	} else if( trainConfig == 7 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
		// closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "1111100040022030000"; PionCutarray[0] = "002016700"; NeutralPionCutarray[0] = "01035033000000"; MesonCutarray[0] = "0103503000000000"; 	
	} else if( trainConfig == 8 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
		// closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "1111100040022030000"; PionCutarray[0] = "002013700"; NeutralPionCutarray[0] = "01035033000000"; MesonCutarray[0] = "0103503000000000"; 			
	} else if( trainConfig == 9 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.75, min pt charged pi = 100 MeV
		// closing neural pion cuts, 0.1 < M_gamma,gamma < 0.145
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "1111100040022030000"; PionCutarray[0] = "002010702"; NeutralPionCutarray[0] = "01035031000000"; MesonCutarray[0] = "0103503000000000"; 
	}

	// PHOS mode
	if( trainConfig == 31 ) {
		// everything open
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "2444400030022000000"; PionCutarray[0] = "000010400"; NeutralPionCutarray[0] = "01035030000000"; MesonCutarray[0] = "0103503000000000"; 
	} else if( trainConfig == 32 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "2444400030022000000"; PionCutarray[0] = "002010700"; NeutralPionCutarray[0] = "01035030000000"; MesonCutarray[0] = "0103503000000000"; 
	} else if( trainConfig == 33 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "2444400030022000000"; PionCutarray[0] = "002013700"; NeutralPionCutarray[0] = "01035030000000"; MesonCutarray[0] = "0103503000000000"; 
	} else if( trainConfig == 34 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "2444400030022000000"; PionCutarray[0] = "002016700"; NeutralPionCutarray[0] = "01035030000000"; MesonCutarray[0] = "0103503000000000"; 
	} else if( trainConfig == 35 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
		// closing neural pion cuts, 0.1 < M_gamma,gamma < 0.145
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "2444400030022000000"; PionCutarray[0] = "002016700"; NeutralPionCutarray[0] = "01035031000000"; MesonCutarray[0] = "0103503000000000"; 	
	} else if( trainConfig == 36 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
		// closing neural pion cuts, 0.11 < M_gamma,gamma < 0.145
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "2444400030022000000"; PionCutarray[0] = "002016700"; NeutralPionCutarray[0] = "01035032000000"; MesonCutarray[0] = "0103503000000000"; 	
	} else if( trainConfig == 37 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
		// closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "2444400030022000000"; PionCutarray[0] = "002016700"; NeutralPionCutarray[0] = "01035033000000"; MesonCutarray[0] = "0103503000000000"; 	
	} else if( trainConfig == 38 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
		// closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "2444400030022000000"; PionCutarray[0] = "002013700"; NeutralPionCutarray[0] = "01035033000000"; MesonCutarray[0] = "0103503000000000"; 			
	} else if( trainConfig == 39 ) {
		// closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.75, min pt charged pi = 100 MeV
		// closing neural pion cuts, 0.1 < M_gamma,gamma < 0.145
		eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; ClusterCutarray[0] = "2444400030022000000"; PionCutarray[0] = "002010702"; NeutralPionCutarray[0] = "01035031000000"; MesonCutarray[0] = "0103503000000000"; 
	}

	TList *EventCutList = new TList();
	TList *ConvCutList  = new TList();
	TList *ClusterCutList  = new TList();
	TList *NeutralPionCutList = new TList();
	TList *MesonCutList = new TList();
	TList *PionCutList  = new TList();

	TList *HeaderList = new TList();
	TObjString *Header1 = new TObjString("pi0_1");
	HeaderList->Add(Header1);
	TObjString *Header3 = new TObjString("eta_2");
	HeaderList->Add(Header3);
	
	EventCutList->SetOwner(kTRUE);
	AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
	ConvCutList->SetOwner(kTRUE);
	AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
	ClusterCutList->SetOwner(kTRUE);
	AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
	NeutralPionCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisNeutralPionCuts   = new AliConversionMesonCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
	PionCutList->SetOwner(kTRUE);
	AliPrimaryPionCuts **analysisPionCuts     = new AliPrimaryPionCuts*[numberOfCuts];

	for(Int_t i = 0; i<numberOfCuts; i++){
		analysisEventCuts[i] = new AliConvEventCuts();   
		analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
		EventCutList->Add(analysisEventCuts[i]);
		analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

		analysisCuts[i] = new AliConversionPhotonCuts();
		if( ! analysisCuts[i]->InitializeCutsFromCutString(ConvCutarray[i].Data()) ) {
				cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
				return 0;
		} else {				
			ConvCutList->Add(analysisCuts[i]);
			analysisCuts[i]->SetFillCutHistograms("",kFALSE);	
		}

		analysisClusterCuts[i] = new AliCaloPhotonCuts();
		if( ! analysisClusterCuts[i]->InitializeCutsFromCutString(ClusterCutarray[i].Data()) ) {
				cout<<"ERROR: analysisClusterCuts [" <<i<<"]"<<endl;
				return 0;
		} else {				
			ClusterCutList->Add(analysisClusterCuts[i]);
			analysisClusterCuts[i]->SetFillCutHistograms("");			
		}

		analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
		if( ! analysisNeutralPionCuts[i]->InitializeCutsFromCutString(NeutralPionCutarray[i].Data()) ) {
			cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
			return 0;
		} else {
			NeutralPionCutList->Add(analysisNeutralPionCuts[i]);
			analysisNeutralPionCuts[i]->SetFillCutHistograms("");
		}
	
		analysisMesonCuts[i] = new AliConversionMesonCuts();
		if( ! analysisMesonCuts[i]->InitializeCutsFromCutString(MesonCutarray[i].Data()) ) {
			cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
			return 0;
		} else {
			MesonCutList->Add(analysisMesonCuts[i]);
			analysisMesonCuts[i]->SetFillCutHistograms("");
		}
		analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
		
		TString cutName( Form("%s_%s_%s_%s_%s_%s",eventCutArray[i].Data(), ConvCutarray[i].Data(), ClusterCutarray[i].Data(),PionCutarray[i].Data(),NeutralPionCutarray[i].Data(), MesonCutarray[i].Data() ) );
		analysisPionCuts[i] = new AliPrimaryPionCuts();
		if( !analysisPionCuts[i]->InitializeCutsFromCutString(PionCutarray[i].Data())) {
			cout<< "ERROR:  analysisPionCuts [ " <<i<<" ] "<<endl;
			return 0;
		} else { 
			PionCutList->Add(analysisPionCuts[i]);
			analysisPionCuts[i]->SetFillCutHistograms("",kFALSE,cutName); 
		}
	}

	task->SetNeutralPionMode(neutralPionMode);
	task->SetEventCutList(numberOfCuts,EventCutList);
	task->SetConversionCutList(ConvCutList);
	task->SetClusterCutList(ClusterCutList);
	task->SetNeutralPionCutList(NeutralPionCutList);
	task->SetMesonCutList(MesonCutList);
	task->SetPionCutList(PionCutList);

	task->SetMoveParticleAccordingToVertex(kTRUE);

	if(enableQAMesonTask) task->SetDoMesonQA(kTRUE);

	//connect containers
	AliAnalysisDataContainer *coutput =
	mgr->CreateContainer(Form("GammaConvNeutralMesonPiPlPiMiPiZero_%i_%i",neutralPionMode, trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiPiZero_%i_%i.root",neutralPionMode, trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
