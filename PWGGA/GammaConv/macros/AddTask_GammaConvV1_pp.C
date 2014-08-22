void AddTask_GammaConvV1_pp(  Int_t trainConfig = 1,  										// change different set of cuts
                              Bool_t isMC   = kFALSE, 										// run MC 
                              Int_t enableQAMesonTask = 0, 									// enable meson QA in AliAnalysisTaskGammaConvV1
                              Int_t enableQAPhotonTask = 0, 								// enable photon QA in AliAnalysisTaskGammaConvV1
                              TString fileNameInputForWeighting = "MCSpectraInput.root", 	// path to file for weigting input
                              TString cutnumberAODBranch = "0000000060084001001500000", 	// cutnumber for AOD branch
							  TString periodname = "LHC12f1x" 								// period name
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
	
	//=========  Set Cutnumber for V0Reader ================================
   	TString cutnumberPhoton = "002084000002200000000";
	TString cutnumberEvent = "0000000"; 
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
	//            find input container
	AliAnalysisTaskGammaConvV1 *task=NULL;
	task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 4;

	TString *eventCutArray = new TString[numberOfCuts];
	TString *photonCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	if(trainConfig == 1){
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , only boxes
		eventCutArray[ 1] = "0001012"; photonCutArray[ 1] = "002093663003800000000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD, V0AND , only boxes
		eventCutArray[ 2] = "0000012"; photonCutArray[ 2] = "002093260003800000000"; mesonCutArray[2] = "01631031009000"; //standard cut Gamma pp 2-76TeV , only boxes
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002093260003800000000"; mesonCutArray[3] = "01631031009000"; //standard cut Gamma pp 2-76TeV , only boxes
	} else if (trainConfig == 2) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , only Minbias MC
		eventCutArray[ 1] = "0001011"; photonCutArray[ 1] = "002093663003800000000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD, V0AND
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002093260003800000000"; mesonCutArray[2] = "01631031009000"; //standard cut Gamma pp 2-76TeV
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002093260003800000000"; mesonCutArray[3] = "01631031009000"; //standard cut Gamma pp 2-76TeV
	} else if (trainConfig == 3) {
		eventCutArray[ 0] = "0002011"; photonCutArray[ 0] = "002093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , only Minbias MC
		eventCutArray[ 1] = "0003011"; photonCutArray[ 1] = "002093663003800000000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND , only Minbias MC
		eventCutArray[ 2] = "0002012"; photonCutArray[ 2] = "002093663003800000000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , only Boxes MC
		eventCutArray[ 3] = "0003012"; photonCutArray[ 3] = "002093663003800000000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND, only Boxes MC
	} else if (trainConfig == 4) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , all photon qualities
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002093663003800020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 1
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002093663003800030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 2
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002093663003800040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 3
	} else if (trainConfig == 5) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "007093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , all photon qualities, min R = 35 cm
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "007093663003800020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 1, min R = 35 cm
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "007093663003800030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 2, min R = 35 cm
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "007093663003800040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 3, min R = 35 cm
	} else if (trainConfig == 6) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002083663003200000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, all photon qualities
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002083663003200020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 1
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002083663003200030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 2
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002083663003200040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 3   
	} else if (trainConfig == 7) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "007083663003200000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, all photon qualities, min R = 35 cm
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "007083663003200020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 1, min R = 35 cm
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "007083663003200030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 2, min R = 35 cm
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "007083663003200040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 3, min R = 35 cm
	} else if (trainConfig == 8) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002083663000200000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 7TeV, all photon qualities
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002083663000200020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 1
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002083663000200030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 2
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002083663000200040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 3   
	} else if (trainConfig == 9) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "007083663000200000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 7TeV, all photon qualities, min R = 35 cm
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "007083663000200020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 1, min R = 35 cm
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "007083663000200030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 2, min R = 35 cm
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "007083663000200040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 3, min R = 35 cm	   
	} else if (trainConfig == 10) {
		eventCutArray[ 0] = "0002011"; photonCutArray[ 0] = "002093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , all photon qualities
		eventCutArray[ 1] = "0002011"; photonCutArray[ 1] = "002093663003800020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 1
		eventCutArray[ 2] = "0002011"; photonCutArray[ 2] = "002093663003800030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 2
		eventCutArray[ 3] = "0002011"; photonCutArray[ 3] = "002093663003800040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 3
	} else if (trainConfig == 11) {
		eventCutArray[ 0] = "0002011"; photonCutArray[ 0] = "007093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , all photon qualities, min R = 35 cm
		eventCutArray[ 1] = "0002011"; photonCutArray[ 1] = "007093663003800020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 1, min R = 35 cm
		eventCutArray[ 2] = "0002011"; photonCutArray[ 2] = "007093663003800030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 2, min R = 35 cm
		eventCutArray[ 3] = "0002011"; photonCutArray[ 3] = "007093663003800040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 3, min R = 35 cm
	} else if (trainConfig == 12) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[0] = "01525065000000"; //standard cut LHC11h pp 2.76TeV 
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "032092970028250400000"; mesonCutArray[1] = "01525065000000"; //variation eta 0.65
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "042092970028250400000"; mesonCutArray[2] = "01525065000000"; //variation eta 0.75
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002092950028250400000"; mesonCutArray[3] = "01525065000000"; //variation pion p dEdx 0.3-5.
	} else if (trainConfig == 13) { //added signals
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[0] = "01525065000000"; //standard cut LHC11h pp 2.76TeV 
		eventCutArray[ 1] = "0000012"; photonCutArray[ 1] = "032092970028250400000"; mesonCutArray[1] = "01525065000000"; //variation eta 0.65
		eventCutArray[ 2] = "0000012"; photonCutArray[ 2] = "042092970028250400000"; mesonCutArray[2] = "01525065000000"; //variation eta 0.75
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002092950028250400000"; mesonCutArray[3] = "01525065000000"; //variation pion p dEdx 0.3-5.
	} else if (trainConfig == 14) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002492970028250400000"; mesonCutArray[0] = "01525065000000"; //variation pt 0.075 
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002192970028250400000"; mesonCutArray[1] = "01525065000000"; //variation pt 0.1
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002062970028250400000"; mesonCutArray[2] = "01525065000000"; //variation TPC cls 0.7
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002082970028250400000"; mesonCutArray[3] = "01525065000000"; //variation TPC cls 0.35 
	} else if (trainConfig == 15) { //added signals
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002492970028250400000"; mesonCutArray[0] = "01525065000000"; //variation pt 0.075 
		eventCutArray[ 1] = "0000012"; photonCutArray[ 1] = "002192970028250400000"; mesonCutArray[1] = "01525065000000"; //variation pt 0.1
		eventCutArray[ 2] = "0000012"; photonCutArray[ 2] = "002062970028250400000"; mesonCutArray[2] = "01525065000000"; //variation TPC cls 0.7
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002082970028250400000"; mesonCutArray[3] = "01525065000000"; //variation TPC cls 0.35 
	} else if (trainConfig == 16) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002093970028250400000"; mesonCutArray[0] = "01525065000000"; //variation edEdx -4,5
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002096970028250400000"; mesonCutArray[1] = "01525065000000"; //variation edEdx -2.5,4
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002092970038250400000"; mesonCutArray[2] = "01525065000000"; //variation TOF el. PID -3,5
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002092970048250400000"; mesonCutArray[3] = "01525065000000"; //variation TOF el. PID -2,3
	} else if (trainConfig == 17) { //added signals
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002093970028250400000"; mesonCutArray[0] = "01525065000000"; //variation edEdx -4,5
		eventCutArray[ 1] = "0000012"; photonCutArray[ 1] = "002096970028250400000"; mesonCutArray[1] = "01525065000000"; //variation edEdx -2.5,4
		eventCutArray[ 2] = "0000012"; photonCutArray[ 2] = "002092970038250400000"; mesonCutArray[2] = "01525065000000"; //variation TOF el. PID -3,5
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002092970048250400000"; mesonCutArray[3] = "01525065000000"; //variation TOF el. PID -2,3
	} else if (trainConfig == 18) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002092970029250400000"; mesonCutArray[0] = "01525065000000"; //variation qt 0.03
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002092970022250400000"; mesonCutArray[1] = "01525065000000"; //variation qt 0.07 no2D
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002092970028150400000"; mesonCutArray[2] = "01525065000000"; //variation chi2 50.
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002092970028850400000"; mesonCutArray[3] = "01525065000000"; //variation chi2 20.
	} else if (trainConfig == 19) { //added signals
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002092970029250400000"; mesonCutArray[0] = "01525065000000"; //variation qt 0.03
		eventCutArray[ 1] = "0000012"; photonCutArray[ 1] = "002092970022250400000"; mesonCutArray[1] = "01525065000000"; //variation qt 0.07 no2D
		eventCutArray[ 2] = "0000012"; photonCutArray[ 2] = "002092970028150400000"; mesonCutArray[2] = "01525065000000"; //variation chi2 50.
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002092970028850400000"; mesonCutArray[3] = "01525065000000"; //variation chi2 20.
	} else if (trainConfig == 20) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002092970028260400000"; mesonCutArray[0] = "01525065000000"; //variation psi pair 0.05
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002092970028280400000"; mesonCutArray[1] = "01525065000000"; //variation psi pair 0.2
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002092970028250000000"; mesonCutArray[2] = "01525065000000"; //variation cosPA -1
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[3] = "01525055000000"; //variation alpha 0.75
	} else if (trainConfig == 21) { //added signals
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002092970028260400000"; mesonCutArray[0] = "01525065000000"; //variation psi pair 0.05
		eventCutArray[ 1] = "0000012"; photonCutArray[ 1] = "002092970028280400000"; mesonCutArray[1] = "01525065000000"; //variation psi pair 0.2
		eventCutArray[ 2] = "0000012"; photonCutArray[ 2] = "002092970028250000000"; mesonCutArray[2] = "01525065000000"; //variation cosPA -1
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[3] = "01525055000000"; //variation alpha 0.75
	} else if (trainConfig == 22) {
		eventCutArray[ 0] = "0004011"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kTRD
		eventCutArray[ 1] = "0005011"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMC
		eventCutArray[ 2] = "0006011"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kPHI
		eventCutArray[ 3] = "0007011"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kHighMult
	} else if (trainConfig == 23) {
		eventCutArray[ 0] = "0008011"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMCEGA
		eventCutArray[ 1] = "0009011"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMCEJE
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[2] = "01525065000000"; // minimum bias
		eventCutArray[ 3] = "0001111"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kINT8
	} else if (trainConfig == 24) {
		eventCutArray[ 0] = "0004211"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kTRD CINT8 HEE
		eventCutArray[ 1] = "0004411"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kTRD CINT8 HSE
		eventCutArray[ 2] = "0004611"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kTRD CINT8 HJE
		eventCutArray[ 3] = "0004811"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kTRD CINT8 HQU
	} else if (trainConfig == 25) {
		eventCutArray[ 0] = "0004111"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kTRD CINT7 HEE
		eventCutArray[ 1] = "0004311"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kTRD CINT7 HSE
		eventCutArray[ 2] = "0004511"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kTRD CINT7 HJE
		eventCutArray[ 3] = "0004711"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kTRD CINT7 HQU
	} else if (trainConfig == 26) {
		eventCutArray[ 0] = "0005211"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMC7
		eventCutArray[ 1] = "0005311"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMC8
		eventCutArray[ 2] = "0006211"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kPHI7
		eventCutArray[ 3] = "0006311"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kPHI8
	} else if (trainConfig == 27) {
		eventCutArray[ 0] = "0005111"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMC1
		eventCutArray[ 1] = "0007111"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kSHM1
		eventCutArray[ 2] = "0007211"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kSHM7
		eventCutArray[ 3] = "0007311"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kSHM8
	} else if (trainConfig == 28) {
		eventCutArray[ 0] = "0008111"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMCEGA + CINT7
		eventCutArray[ 1] = "0008211"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMCEGA + CINT8
		eventCutArray[ 2] = "0008311"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kEMCEG1 + CINT7
		eventCutArray[ 3] = "0008411"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kEMCEG1 + CINT8
	} else if (trainConfig == 29) {
		eventCutArray[ 0] = "0008511"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMCEG2 + CINT7
		eventCutArray[ 1] = "0008611"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMCEG2 + CINT8
		eventCutArray[ 2] = "0009111"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kEMCEJE + CINT7
		eventCutArray[ 3] = "0009211"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kEMCEJE + CINT8
	} else if (trainConfig == 30) {
		eventCutArray[ 0] = "0009311"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMCEJ1 + CINT7
		eventCutArray[ 1] = "0009411"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMCEJ1 + CINT8
		eventCutArray[ 2] = "0009511"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kEMCEJ2 + CINT7
		eventCutArray[ 3] = "0009611"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kEMCEJ2 + CINT8		
	} else if (trainConfig == 31) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002092570028250400000"; mesonCutArray[0] = "01521065000000"; //new standard cut for pp 8 TeV
		eventCutArray[ 0] = "0000011"; photonCutArray[ 1] = "002093570028250400000"; mesonCutArray[0] = "01521065000000"; //variation edEdx -4,5
		eventCutArray[ 1] = "0000011"; photonCutArray[ 2] = "002096570028250400000"; mesonCutArray[1] = "01521065000000"; //variation edEdx -2.5,4
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002092550028250400000"; mesonCutArray[3] = "01521065000000"; //variation pion p dEdx 0.3-5.
	} else if (trainConfig == 32) { //added signals
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002092570028250400000"; mesonCutArray[0] = "01521065000000"; //new standard cut for pp 8 TeV
		eventCutArray[ 0] = "0000012"; photonCutArray[ 1] = "002093570028250400000"; mesonCutArray[0] = "01521065000000"; //variation edEdx -4,5
		eventCutArray[ 1] = "0000012"; photonCutArray[ 2] = "002096570028250400000"; mesonCutArray[1] = "01521065000000"; //variation edEdx -2.5,4
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002092550028250400000"; mesonCutArray[3] = "01521065000000"; //variation pion p dEdx 0.3-5.
	} else if (trainConfig == 33) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002492570028250400000"; mesonCutArray[0] = "01521065000000"; //variation pt 0.075 
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002192570028250400000"; mesonCutArray[1] = "01521065000000"; //variation pt 0.1
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002062570028250400000"; mesonCutArray[2] = "01521065000000"; //variation TPC cls 0.7
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002082570028250400000"; mesonCutArray[3] = "01521065000000"; //variation TPC cls 0.35 
	} else if (trainConfig == 34) { //added signals
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002492570028250400000"; mesonCutArray[0] = "01521065000000"; //variation pt 0.075 
		eventCutArray[ 1] = "0000012"; photonCutArray[ 1] = "002192570028250400000"; mesonCutArray[1] = "01521065000000"; //variation pt 0.1
		eventCutArray[ 2] = "0000012"; photonCutArray[ 2] = "002062570028250400000"; mesonCutArray[2] = "01521065000000"; //variation TPC cls 0.7
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002082570028250400000"; mesonCutArray[3] = "01521065000000"; //variation TPC cls 0.35 
	} else if (trainConfig == 35) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002092270028250400000"; mesonCutArray[0] = "01521065000000"; //variation pidEdx 1,-10
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002092370028250400000"; mesonCutArray[1] = "01521065000000"; //variation pidEdx 2.5,-10
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[2] = "01521065000000"; //variation pidEdx 3,-10
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002092500028250400000"; mesonCutArray[3] = "01521065000000"; //variation pion p dEdx 0.5-5
	} else if (trainConfig == 36) { //added signals
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002092270028250400000"; mesonCutArray[0] = "01521065000000"; //variation pidEdx 1,-10
		eventCutArray[ 1] = "0000012"; photonCutArray[ 1] = "002092370028250400000"; mesonCutArray[1] = "01521065000000"; //variation pidEdx 2.5,-10
		eventCutArray[ 2] = "0000012"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[2] = "01521065000000"; //variation pidEdx 3,-10
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002092500028250400000"; mesonCutArray[3] = "01521065000000"; //variation pion p dEdx 0.5-5
	} else if (trainConfig == 37) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002092570029250400000"; mesonCutArray[0] = "01521065000000"; //variation qt 0.03
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002092570022250400000"; mesonCutArray[1] = "01521065000000"; //variation qt 0.07 no2D
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002092570028150400000"; mesonCutArray[2] = "01521065000000"; //variation chi2 50.
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002092570028850400000"; mesonCutArray[3] = "01521065000000"; //variation chi2 20.
	} else if (trainConfig == 38) { //added signals
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002092570029250400000"; mesonCutArray[0] = "01521065000000"; //variation qt 0.03
		eventCutArray[ 1] = "0000012"; photonCutArray[ 1] = "002092570022250400000"; mesonCutArray[1] = "01521065000000"; //variation qt 0.07 no2D
		eventCutArray[ 2] = "0000012"; photonCutArray[ 2] = "002092570028150400000"; mesonCutArray[2] = "01521065000000"; //variation chi2 50.
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002092570028850400000"; mesonCutArray[3] = "01521065000000"; //variation chi2 20.
	} else if (trainConfig == 39) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002092570028260400000"; mesonCutArray[0] = "01521065000000"; //variation psi pair 0.05
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002092570028280400000"; mesonCutArray[1] = "01521065000000"; //variation psi pair 0.2
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002092570028250000000"; mesonCutArray[2] = "01521065000000"; //variation cosPA -1
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002092570028250600000"; mesonCutArray[3] = "01521065000000"; //variation cosPA 0.9
	} else if (trainConfig == 40) { //added signals
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002092570028260400000"; mesonCutArray[0] = "01521065000000"; //variation psi pair 0.05
		eventCutArray[ 1] = "0000012"; photonCutArray[ 1] = "002092570028280400000"; mesonCutArray[1] = "01521065000000"; //variation psi pair 0.2
		eventCutArray[ 2] = "0000012"; photonCutArray[ 2] = "002092570028250000000"; mesonCutArray[2] = "01521065000000"; //variation cosPA -1
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002092570028250600000"; mesonCutArray[3] = "01521065000000"; //variation cosPA 0.9
	} else if (trainConfig == 41) {
		eventCutArray[ 0] = "0000011"; photonCutArray[ 0] = "002092570028950400000"; mesonCutArray[0] = "01521065000000"; //variation chi2 15
		eventCutArray[ 1] = "0000011"; photonCutArray[ 1] = "002092570028230400000"; mesonCutArray[1] = "01521065000000"; //variation psi pair 0.035
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002092570028250400000"; mesonCutArray[2] = "01521055000000"; //variation alpha 0.75
		eventCutArray[ 3] = "0000011"; photonCutArray[ 3] = "002092570028250400000"; mesonCutArray[3] = "01521075000000"; //variation alpha 0.85
	} else if (trainConfig == 42) { //added signals
		eventCutArray[ 0] = "0000012"; photonCutArray[ 0] = "002092570028950400000"; mesonCutArray[0] = "01521065000000"; //variation chi2 15
		eventCutArray[ 1] = "0000012"; photonCutArray[ 1] = "002092570028230400000"; mesonCutArray[1] = "01521065000000"; //variation psi pair 0.035
		eventCutArray[ 2] = "0000012"; photonCutArray[ 2] = "002092570028250400000"; mesonCutArray[2] = "01521055000000"; //variation alpha 0.75
		eventCutArray[ 3] = "0000012"; photonCutArray[ 3] = "002092570028250400000"; mesonCutArray[3] = "01521075000000"; //variation alpha 0.85
	} else if (trainConfig == 43) {
		eventCutArray[ 0] = "0004011"; photonCutArray[ 0] = "002092570028250400000"; mesonCutArray[0] = "01521065000000"; // trigger kTRD with y 0.8
		eventCutArray[ 1] = "0005011"; photonCutArray[ 1] = "002092570028250400000"; mesonCutArray[1] = "01521065000000"; // trigger kEMC with y 0.8
		eventCutArray[ 2] = "0006011"; photonCutArray[ 2] = "002092570028250400000"; mesonCutArray[2] = "01521065000000"; // trigger kPHI with y 0.8
		eventCutArray[ 3] = "0007011"; photonCutArray[ 3] = "002092570028250400000"; mesonCutArray[3] = "01521065000000"; // trigger kHighMult with y 0.8
	} else if (trainConfig == 44) {
		eventCutArray[ 0] = "0008011"; photonCutArray[ 0] = "002092570028250400000"; mesonCutArray[0] = "01521065000000"; // trigger kEMCEGA with y 0.8
		eventCutArray[ 1] = "0009011"; photonCutArray[ 1] = "002092570028250400000"; mesonCutArray[1] = "01521065000000"; // trigger kEMCEJE with y 0.8
		eventCutArray[ 2] = "0000011"; photonCutArray[ 2] = "002092570028250400000"; mesonCutArray[2] = "01521065000000"; // minimum bias with y 0.8
		eventCutArray[ 3] = "0001111"; photonCutArray[ 3] = "002092570028250400000"; mesonCutArray[3] = "01521065000000"; // trigger kINT8 with y 0.8
	}
	
	 else {
			Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
			return;
	}

	TList *EventCutList = new TList();
	TList *ConvCutList = new TList();
	TList *MesonCutList = new TList();

	TList *HeaderList = new TList();
	if (periodname.CompareTo("LHC12i3") == 0){	
		TObjString *Header2 = new TObjString("BOX");
		HeaderList->Add(Header2);
	} else if (periodname.CompareTo("LHC14e2b")==0){
		TObjString *Header2 = new TObjString("pi0_1");
		HeaderList->Add(Header2);
		TObjString *Header3 = new TObjString("eta_2");
		HeaderList->Add(Header3);
	}	
		
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

	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
