void AddTask_GammaConvV1_pPb( 	Int_t 		trainConfig = 1,  												// change different set of cuts
								Bool_t	 	isMC   = kFALSE, 												// run MC
								Int_t 		enableQAMesonTask = 0, 											// enable QA in AliAnalysisTaskGammaConvV1
								Int_t 		enableQAPhotonTask = 0, 										// enable additional QA task
								TString 	fileNameInputForWeighting 	= "MCSpectraInput.root", 			// path to file for weigting input
								Int_t 		doWeightingPart 			= 0, 			 					// enable Weighting
								TString 	generatorName 				= "DPMJET",							// generator Name
                                TString 	cutnumberAODBranch 			= "800000006008400000150000000", 	// cutnumber for AOD branch
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

	//========= Add PID Reponse to ANALYSIS manager ====
	if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
		gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
		AddTaskPIDResponse(isMC);
	}

	//=========  Set Cutnumber for V0Reader ================================
	TString cutnumberPhoton = "06000008000100001500000000";
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
			if(trainConfig==193 || trainConfig==194 || trainConfig==195 || trainConfig==196){
			  fCuts->SetDodEdxSigmaCut(kFALSE);
			}
			if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
				fV0ReaderV1->SetConversionCuts(fCuts);
				fCuts->SetFillCutHistograms("",kTRUE);
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
	Int_t numberOfCuts = 4;

	TString *photonCutArray = new TString[numberOfCuts];
	TString *eventCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];
	Bool_t doEtaShiftIndCuts = kFALSE;
	TString stringShift = "";
	
    if(trainConfig == 1){
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 2) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 3) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 4) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 5) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts ///New STANDARD CUT
	} else if (trainConfig == 6) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts ///New STANDARD CUT
	} else if (trainConfig == 7) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 8) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 9) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 10) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 11) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 12) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 13) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 14) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 15) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 16) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 17) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 18) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 19) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 20) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 21) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 22) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 23){
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 24) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 25) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 26) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 27) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts
	} else if (trainConfig == 28) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts
	} else if (trainConfig == 29) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 30) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 31) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 32) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 33) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 34) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 35) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 36) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 37) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 38) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 39) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 40) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 41) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 42) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 43) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 44) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 45){
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "82400113"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "82400113"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "82400113"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 46) {
		eventCutArray[ 0] = "82400123"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "82400123"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "82400123"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "82400123"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 47) {
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "82400113"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "82400113"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "82400113"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 48) {
		eventCutArray[ 0] = "82400123"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "82400123"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "82400123"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "82400123"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 49) {
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "82400113"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "82400113"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "82400113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts
	} else if (trainConfig == 50) {
		eventCutArray[ 0] = "82400123"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "82400123"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "82400123"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "82400123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts
	} else if (trainConfig == 51) {
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "82400113"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "82400113"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "82400113"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 52) {
		eventCutArray[ 0] = "82400123"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "82400123"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "82400123"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "82400123"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 53) {
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "82400113"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "82400113"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "82400113"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 54) {
		eventCutArray[ 0] = "82400123"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "82400123"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "82400123"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "82400123"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 55) {
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "82400113"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "82400113"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "82400113"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 56) {
		eventCutArray[ 0] = "82400123"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "82400123"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "82400123"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "82400123"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 57) {
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "82400113"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "82400113"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "82400113"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 58) {
		eventCutArray[ 0] = "82400123"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "82400123"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "82400123"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "82400123"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 59) {
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "82400113"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "82400113"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "82400113"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 60) {
		eventCutArray[ 0] = "82400123"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "82400123"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "82400123"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "82400123"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 61) {
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "82400113"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "82400113"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "82400113"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 62) {
		eventCutArray[ 0] = "82400123"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "82400123"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "82400123"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "82400123"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 63) {
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "82400113"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "82400113"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "82400113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 64) {
		eventCutArray[ 0] = "82400123"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "82400123"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "82400123"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "82400123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 65) {
		eventCutArray[ 0] = "82400113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "82400113"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "82400113"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "82400113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 66) {
		eventCutArray[ 0] = "82400123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "82400123"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "82400123"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "82400123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 67){
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "84600113"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "84600113"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "84600113"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 68) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 69) {
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "84600113"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "84600113"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "84600113"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 70) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 71) {
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "84600113"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "84600113"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "84600113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts
	} else if (trainConfig == 72) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts
	} else if (trainConfig == 73) {
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "84600113"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "84600113"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "84600113"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 74) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 75) {
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "84600113"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "84600113"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "84600113"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 76) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 77) {
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "84600113"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "84600113"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "84600113"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 78) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 79) {
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "84600113"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "84600113"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "84600113"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 80) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 81) {
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "84600113"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "84600113"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "84600113"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 82) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 83) {
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "84600113"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "84600113"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "84600113"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 84) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 85) {
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "84600113"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "84600113"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "84600113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 86) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 87) {
		eventCutArray[ 0] = "84600113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "84600113"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "84600113"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "84600113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 88) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 89){
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "86800113"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "86800113"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "86800113"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 90) {
		eventCutArray[ 0] = "84600123"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "84600123"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "84600123"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "84600123"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 91) {
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "86800113"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "86800113"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "86800113"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 92) {
		eventCutArray[ 0] = "86800123"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "86800123"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "86800123"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "86800123"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 93) {
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "86800113"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "86800113"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "86800113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts
	} else if (trainConfig == 94) {
		eventCutArray[ 0] = "86800123"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "86800123"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "86800123"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "86800123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts
	} else if (trainConfig == 95) {
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "86800113"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "86800113"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "86800113"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 96) {
		eventCutArray[ 0] = "86800123"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "86800123"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "86800123"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "86800123"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 97) {
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "86800113"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "86800113"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "86800113"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 98) {
		eventCutArray[ 0] = "86800123"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "86800123"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "86800123"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "86800123"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 99) {
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "86800113"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "86800113"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "86800113"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 100) {
		eventCutArray[ 0] = "86800123"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "86800123"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "86800123"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "86800123"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 101) {
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "86800113"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "86800113"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "86800113"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 102) {
		eventCutArray[ 0] = "86800123"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "86800123"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "86800123"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "86800123"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 103) {
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "86800113"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "86800113"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "86800113"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 104) {
		eventCutArray[ 0] = "86800123"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "86800123"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "86800123"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "86800123"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 105) {
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "86800113"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "86800113"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "86800113"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 106) {
		eventCutArray[ 0] = "86800123"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "86800123"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "86800123"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "86800123"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 107) {
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "86800113"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "86800113"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "86800113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 108) {
		eventCutArray[ 0] = "86800123"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "86800123"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "86800123"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "86800123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 109) {
		eventCutArray[ 0] = "86800113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "86800113"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "86800113"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "86800113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 110) {
		eventCutArray[ 0] = "86800123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "86800123"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "86800123"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "86800123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 111){
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "86000113"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "86000113"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "86000113"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 112) {
		eventCutArray[ 0] = "86000123"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500900000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "86000123"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500900000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "86000123"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500900000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "86000123"; photonCutArray[ 3] = "00200009237002003220000000"; mesonCutArray[ 3] = "0162103500900000";
	} else if (trainConfig == 113) {
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "86000113"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "86000113"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "86000113"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 114) {
		eventCutArray[ 0] = "86000123"; photonCutArray[ 0] = "00200009277002003220000000"; mesonCutArray[ 0] = "0162103500900000";
		eventCutArray[ 1] = "86000123"; photonCutArray[ 1] = "00200009255102003220000000"; mesonCutArray[ 1] = "0162103500900000";
		eventCutArray[ 2] = "86000123"; photonCutArray[ 2] = "00200009217000003220000000"; mesonCutArray[ 2] = "0162103500900000"; //just tighten Psi pair
		eventCutArray[ 3] = "86000123"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500900000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 115) {
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "86000113"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "86000113"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "86000113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts
	} else if (trainConfig == 116) {
		eventCutArray[ 0] = "86000123"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500900000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "86000123"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500900000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "86000123"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500900000"; //clean cuts
		eventCutArray[ 3] = "86000123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //clean cuts
	} else if (trainConfig == 117) {
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "86000113"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "86000113"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "86000113"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 118) {
		eventCutArray[ 0] = "86000123"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "86000123"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "86000123"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "86000123"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 119) {
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "86000113"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "86000113"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "86000113"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 120) {
		eventCutArray[ 0] = "86000123"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "86000123"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "86000123"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "86000123"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 121) {
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "86000113"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "86000113"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "86000113"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 122) {
		eventCutArray[ 0] = "86000123"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "86000123"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "86000123"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "86000123"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 123) {
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "86000113"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "86000113"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "86000113"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 124) {
		eventCutArray[ 0] = "86000123"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "86000123"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "86000123"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "86000123"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 125) {
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "86000113"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "86000113"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "86000113"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 126) {
		eventCutArray[ 0] = "86000123"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "86000123"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "86000123"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "86000123"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 127) {
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "86000113"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "86000113"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "86000113"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 128) {
		eventCutArray[ 0] = "86000123"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "86000123"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "86000123"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "86000123"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500900000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 129) {
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "86000113"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "86000113"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "86000113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 130) {
		eventCutArray[ 0] = "86000123"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500900000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "86000123"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500900000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "86000123"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500900000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "86000123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400900000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 131) {
		eventCutArray[ 0] = "86000113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "86000113"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "86000113"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "86000113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 132) {
		eventCutArray[ 0] = "86000123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500900000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "86000123"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "86000123"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500900000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "86000123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 133) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500000000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500000000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500000000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500000000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 134) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500000000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500000000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500000000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500000000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 135) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500000000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500000000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500000000"; //clean cuts
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500000000"; //clean cuts
	} else if (trainConfig == 136) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500000000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500000000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500000000"; //clean cuts
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500000000"; //clean cuts
	} else if (trainConfig == 137) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 138) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 139) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 140) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 141) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 142) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 143) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 144) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 145) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 146) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 147) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 148) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 149) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500000000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400000000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 150) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500000000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400000000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 151) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500000000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500000000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 152) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500000000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500000000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 153) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500000000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500000000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500000000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500000000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 154) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "03200009217000008260400000"; mesonCutArray[ 0] = "0162303500000000"; //New STANDARD CUT |eta| < 0.65, |y| < 0.6
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "04200009217000008260400000"; mesonCutArray[ 1] = "0162203500000000"; //New STANDARD CUT |eta| < 0.75, |y| < 0.7
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "01200009217000008260400000"; mesonCutArray[ 2] = "0162403500000000"; //New STANDARD CUT |eta| < 0.6, |y| < 0.5
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000003260000000"; mesonCutArray[ 3] = "0162103500000000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 155) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500000000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500000000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500000000"; //clean cuts
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500000000"; //clean cuts
	} else if (trainConfig == 156) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217000008220000000"; mesonCutArray[ 0] = "0162103500000000"; //tighten psi pair and qt in 2D
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009217000008260000000"; mesonCutArray[ 1] = "0162103500000000"; //tighten psi pair and chi2 in 2D and qt in 2D
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009217000008220400000"; mesonCutArray[ 2] = "0162103500000000"; //clean cuts
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500000000"; //clean cuts
	} else if (trainConfig == 157) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 158) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00100009217000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8  // minR 2.8
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00900009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8  //minR 7.5
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200079217000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200019217000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 159) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 160) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200008217000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200006217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009317000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009617000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 161) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 162) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009227000008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009257000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009216000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009215000008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 163) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 164) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217200008260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009216200008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009226000008260400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009226200008260400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 165) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 166) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217000003260400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009217000009260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009217000002000260400"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000008220400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
	} else if (trainConfig == 167) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 168) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217000008160400000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009217000008860400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009217000008250400000"; mesonCutArray[ 2] = "0162103500000000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000008270400000"; mesonCutArray[ 3] = "0162103500000000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 169) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500000000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400000000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 170) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217000008260300000"; mesonCutArray[ 0] = "0162103500000000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009217000008260600000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162106500000000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103400000000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 171) {
		eventCutArray[ 0] = "80200113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500000000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "80200113"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "80200113"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500000000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "80200113"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 172) {
		eventCutArray[ 0] = "80200123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0262103500000000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		eventCutArray[ 1] = "80200123"; photonCutArray[ 1] = "00200009217000008260400000"; mesonCutArray[ 1] = "0162103500000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		eventCutArray[ 2] = "80200123"; photonCutArray[ 2] = "00200009317200003290000000"; mesonCutArray[ 2] = "0162103500000000";  //old standard eta=0.9 y=0.8
		eventCutArray[ 3] = "80200123"; photonCutArray[ 3] = "00200009217000008260400000"; mesonCutArray[ 3] = "0162103500800000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 173) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  // all Photon Qualities
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000008260420000"; mesonCutArray[ 1] = "0162103500900000";  // only photon quality 1
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000008260430000"; mesonCutArray[ 2] = "0162103500900000";  // only photon quality 2
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008260440000"; mesonCutArray[ 3] = "0162103500900000";  // only photon quality 3
	} else if (trainConfig == 174) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  // all Photon Qualities
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000008260420000"; mesonCutArray[ 1] = "0162103500900000";  // only photon quality 1
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000008260430000"; mesonCutArray[ 2] = "0162103500900000";  // only photon quality 2
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008260440000"; mesonCutArray[ 3] = "0162103500900000";  // only photon quality 3
	} else if (trainConfig == 175) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260410000"; mesonCutArray[ 0] = "0162103500900000";  //no shared electron cut on
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200019217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //single pt cut 0.1
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200029217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //single pt cut 0.15
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200049217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //single pt cut 0.075
	} else if (trainConfig == 176) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260410000"; mesonCutArray[ 0] = "0162103500900000";  //no shared electron cut on
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200019217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  //single pt cut 0.1
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200029217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  //single pt cut 0.15
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200049217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  //single pt cut 0.075
	} else if (trainConfig == 177) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00700009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  // min R =35 cm
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00700009217000008260420000"; mesonCutArray[ 1] = "0162103500900000";  // min R =35 cm & only photon quality 1
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00700009217000008260430000"; mesonCutArray[ 2] = "0162103500900000";  // min R =35 cm & only photon quality 2
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00700009217000008260440000"; mesonCutArray[ 3] = "0162103500900000";  // min R =35 cm & only photon quality 3
	} else if (trainConfig == 178) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00700009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  // min R =35 cm
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00700009217000008260420000"; mesonCutArray[ 1] = "0162103500900000";  // min R =35 cm & only photon quality 1
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00700009217000008260430000"; mesonCutArray[ 2] = "0162103500900000";  // min R =35 cm & only photon quality 2
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00700009217000008260440000"; mesonCutArray[ 3] = "0162103500900000";  // min R =35 cm & only photon quality 3
	} else if (trainConfig == 179) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162803500900000"; //new standard rapidity 0-0.5 in cms
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00500009217000008260400000"; mesonCutArray[ 1] = "0162103500900000"; //new standard RCut=10cm
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00800009217000008260400000"; mesonCutArray[ 2] = "0162103500900000"; //new standard RCut=12.5cm
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00600009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //new standard RCut=20cm
	} else if (trainConfig == 180) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162803500900000"; //new standard rapidity 0-0.5 in cms
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00500009217000008260400000"; mesonCutArray[ 1] = "0162103500900000"; //new standard RCut=10cm
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00800009217000008260400000"; mesonCutArray[ 2] = "0162103500900000"; //new standard RCut=12.5cm
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00600009217000008260400000"; mesonCutArray[ 3] = "0162103500900000"; //new standard RCut=20cm
	} else if (trainConfig == 181) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009297002008260400000"; mesonCutArray[ 1] = "0152506500000000";
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009297002008270400000"; mesonCutArray[ 2] = "0152506500000000";
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009297002008250000000"; mesonCutArray[ 3] = "0152506500000000";
	} else if (trainConfig == 182) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009297002008260400000"; mesonCutArray[ 1] = "0152506500000000";
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009297002008270400000"; mesonCutArray[ 2] = "0152506500000000";
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009297002008250000000"; mesonCutArray[ 3] = "0152506500000000";
	} else if (trainConfig == 183) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00000009217000008260440000"; mesonCutArray[ 0] = "0162103500900000";  // min R =0 cm & only photon quality 3
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00100009217000008260440000"; mesonCutArray[ 1] = "0162103500900000";  // min R =2.8 cm & only photon quality 3
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000008260440000"; mesonCutArray[ 2] = "0162103500900000";  // min R =5 cm & only photon quality 3
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00900009217000008260440000"; mesonCutArray[ 3] = "0162103500900000";  // min R =7.5 cm & only photon quality 3 
	} else if (trainConfig == 184) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00000009217000008260440000"; mesonCutArray[ 0] = "0162103500900000";  // min R =0 cm & only photon quality 3
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00100009217000008260440000"; mesonCutArray[ 1] = "0162103500900000";  // min R =2.8 cm & only photon quality 3
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000008260440000"; mesonCutArray[ 2] = "0162103500900000";  // min R =5 cm & only photon quality 3
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00900009217000008260440000"; mesonCutArray[ 3] = "0162103500900000";  // min R =7.5 cm & only photon quality 3
	} else if (trainConfig == 185) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00500009217000008260440000"; mesonCutArray[ 0] = "0162103500900000";  // min R =10 cm & only photon quality 3
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00800009217000008260440000"; mesonCutArray[ 1] = "0162103500900000";  // min R =12.5 cm & only photon quality 3
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00600009217000008260440000"; mesonCutArray[ 2] = "0162103500900000";  // min R =20 cm & only photon quality 3
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00700009217000008260440000"; mesonCutArray[ 3] = "0162103500900000";  // min R =35 cm & only photon quality 3 
	} else if (trainConfig == 186) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00500009217000008260440000"; mesonCutArray[ 0] = "0162103500900000";  // min R =10 cm & only photon quality 3
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00800009217000008260440000"; mesonCutArray[ 1] = "0162103500900000";  // min R =12.5 cm & only photon quality 3
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00600009217000008260440000"; mesonCutArray[ 2] = "0162103500900000";  // min R =20 cm & only photon quality 3
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00700009217000008260440000"; mesonCutArray[ 3] = "0162103500900000";  // min R =35 cm & only photon quality 3  
	} else if (trainConfig == 187) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00000009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  // min R =0 cm
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00100009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  // min R =2.8 cm
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  // min R =5 cm
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00900009217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  // min R =7.5 cm
	} else if (trainConfig == 188) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00000009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  // min R =0 cm
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00100009217000008260400000"; mesonCutArray[ 1] = "0162103500900000";  // min R =2.8 cm
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000008260400000"; mesonCutArray[ 2] = "0162103500900000";  // min R =5 cm
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00900009217000008260400000"; mesonCutArray[ 3] = "0162103500900000";  // min R =7.5 cm 
	} else if (trainConfig == 189) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  // standard
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000008260404000"; mesonCutArray[ 1] = "0162103500900000";  // double counting cut case 4
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000008260405000"; mesonCutArray[ 2] = "0162103500900000";  // double counting cut case 5
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008260406000"; mesonCutArray[ 3] = "0162103500900000";  // double counting cut case 6
	} else if (trainConfig == 190) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500900000";  // standard 
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000008260404000"; mesonCutArray[ 1] = "0162103500900000";  // double counting cut case 4
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000008260405000"; mesonCutArray[ 2] = "0162103500900000";  // double counting cut case 5
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008260406000"; mesonCutArray[ 3] = "0162103500900000";  // double counting cut case 6 
	} else if (trainConfig == 191) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500000000";  // standard //eta meson
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009217000008260404000"; mesonCutArray[ 1] = "0162103500000000";  // double counting cut case 4
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009217000008260405000"; mesonCutArray[ 2] = "0162103500000000";  // double counting cut case 5
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009217000008260406000"; mesonCutArray[ 3] = "0162103500000000";  // double counting cut case 6
	} else if (trainConfig == 192) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260400000"; mesonCutArray[ 0] = "0162103500000000";  // standard //eta meson
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009217000008260404000"; mesonCutArray[ 1] = "0162103500000000";  // double counting cut case 4
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009217000008260405000"; mesonCutArray[ 2] = "0162103500000000";  // double counting cut case 5
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009217000008260406000"; mesonCutArray[ 3] = "0162103500000000";  // double counting cut case 6 
	} else if (trainConfig == 193) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260404000"; mesonCutArray[ 0] = "0162103500900000";  // standard //pi0 dc cut4
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009000000008260404000"; mesonCutArray[ 1] = "0162103500900000";  // open dEdx Cut -10;10
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009400000008260404000"; mesonCutArray[ 2] = "0162103500900000";  // open dEdx cut; -6;7
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009000100008260404000"; mesonCutArray[ 3] = "0162103500900000";  // no dEdx Cut
	} else if (trainConfig == 194) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260404000"; mesonCutArray[ 0] = "0162103500900000";  // standard //pi0 dc cut4
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009000000008260404000"; mesonCutArray[ 1] = "0162103500900000";  //  open dEdx Cut -10;10
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009400000008260404000"; mesonCutArray[ 2] = "0162103500900000";  //  open dEdx cut; -6;7
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009000100008260404000"; mesonCutArray[ 3] = "0162103500900000";  //no dEdx Cut 
	} else if (trainConfig == 195) {
		eventCutArray[ 0] = "80000113"; photonCutArray[ 0] = "00200009217000008260404000"; mesonCutArray[ 0] = "0162103500000000";  // standard //eta dc cut4
		eventCutArray[ 1] = "80000113"; photonCutArray[ 1] = "00200009000000008260404000"; mesonCutArray[ 1] = "0162103500000000";  // open dEdx Cut -10;10
		eventCutArray[ 2] = "80000113"; photonCutArray[ 2] = "00200009400000008260404000"; mesonCutArray[ 2] = "0162103500000000";  // open dEdx cut; -6;7
		eventCutArray[ 3] = "80000113"; photonCutArray[ 3] = "00200009000100008260404000"; mesonCutArray[ 3] = "0162103500000000";  // no dEdx Cut
	} else if (trainConfig == 196) {
		eventCutArray[ 0] = "80000123"; photonCutArray[ 0] = "00200009217000008260404000"; mesonCutArray[ 0] = "0162103500000000";  // standard //eta dc cut4
		eventCutArray[ 1] = "80000123"; photonCutArray[ 1] = "00200009000000008260404000"; mesonCutArray[ 1] = "0162103500000000";  //  open dEdx Cut -10;10
		eventCutArray[ 2] = "80000123"; photonCutArray[ 2] = "00200009400000008260404000"; mesonCutArray[ 2] = "0162103500000000";  //  open dEdx cut; -6;7
		eventCutArray[ 3] = "80000123"; photonCutArray[ 3] = "00200009000100008260404000"; mesonCutArray[ 3] = "0162103500000000";  //no dEdx Cut
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
		if ( trainConfig == 1 || trainConfig == 3 || trainConfig == 5 || trainConfig == 7 || trainConfig == 9 || trainConfig == 11 || trainConfig == 13 || trainConfig == 15|| trainConfig == 17||
		     trainConfig == 19 || trainConfig == 21 || trainConfig == 133 || trainConfig == 135 || trainConfig == 137 || trainConfig == 139 || trainConfig == 141 || trainConfig == 143 ||
		     trainConfig == 145 || trainConfig == 147 || trainConfig == 149 || trainConfig == 151 || trainConfig == 173 || trainConfig == 175 || trainConfig == 177 || trainConfig == 179 ||
		     trainConfig == 181 || trainConfig == 183 || trainConfig == 185 || trainConfig == 187 || trainConfig == 189 || trainConfig == 191 || trainConfig == 193 || trainConfig == 195){
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
		if ( trainConfig == 2 || trainConfig == 4 || trainConfig == 6 || trainConfig == 8 || trainConfig == 10 || trainConfig == 12 || trainConfig == 14 || trainConfig == 16|| trainConfig == 18||
		     trainConfig == 20|| trainConfig == 22 || trainConfig == 134 || trainConfig == 136 || trainConfig == 138 || trainConfig == 140 || trainConfig == 142 || trainConfig == 144 || 
		     trainConfig == 146 || trainConfig == 148 || trainConfig == 150 || trainConfig == 152 || trainConfig == 174 || trainConfig == 176 || trainConfig == 178 || trainConfig == 180 ||
		     trainConfig == 182 || trainConfig == 184 || trainConfig == 186 || trainConfig == 188 || trainConfig == 190 || trainConfig == 192 || trainConfig == 194 || trainConfig == 196){
			if (doWeighting){
				analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A",
																			 "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
			}
			
		}   
		if ( trainConfig == 23 || trainConfig == 25 || trainConfig == 27 || trainConfig == 29 || trainConfig == 31 || trainConfig == 33 || trainConfig == 35 || trainConfig == 37|| trainConfig == 39||
			trainConfig == 41 || trainConfig == 43 || trainConfig == 153 || trainConfig == 155 || trainConfig == 157 || trainConfig == 159 || trainConfig == 161 || trainConfig == 163 ||
			trainConfig == 165 || trainConfig == 167 || trainConfig == 169 || trainConfig == 171 ){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_0020000V0A",
																				 "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_0020000V0A", "", "Pi0_Fit_Data_pPb_5023GeV_0020000V0A",
																				 "Eta_Fit_Data_pPb_5023GeV_0020000V0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_0020000V0A",
																				 "Eta_Hijing_LHC13e7_pPb_5023GeV_0020000V0A", "", "Pi0_Fit_Data_pPb_5023GeV_0020000V0A",
																				 "Eta_Fit_Data_pPb_5023GeV_0020000V0A");
				}
			}
		}   
		if ( trainConfig == 24 || trainConfig == 26 || trainConfig == 28 || trainConfig == 30 || trainConfig == 32 || trainConfig == 34 || trainConfig == 36 || trainConfig == 38|| trainConfig == 40||
			trainConfig == 42|| trainConfig == 44 || trainConfig == 154 || trainConfig == 156 || trainConfig == 158 || trainConfig == 160 || trainConfig == 162 || trainConfig == 164 || 
			trainConfig == 166 || trainConfig == 168 || trainConfig == 170 || trainConfig == 172){
			if (doWeighting){
				analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_0020000V0A",
																			 "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_0020000V0A", "", "Pi0_Fit_Data_pPb_5023GeV_0020000V0A",
																			 "Eta_Fit_Data_pPb_5023GeV_0020000V0A");
			}
		}   
		if ( trainConfig == 45 || trainConfig == 47 || trainConfig == 49 || trainConfig == 51 || trainConfig == 53 || trainConfig == 55 || trainConfig == 57 || trainConfig == 59|| trainConfig == 61||
			trainConfig == 63 || trainConfig == 65 ){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_2040V0A",
																				 "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_2040V0A", "", "Pi0_Fit_Data_pPb_5023GeV_2040V0A",
																				 "Eta_Fit_Data_pPb_5023GeV_2040V0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_2040V0A",
																				 "Eta_Hijing_LHC13e7_pPb_5023GeV_2040V0A", "", "Pi0_Fit_Data_pPb_5023GeV_2040V0A", "Eta_Fit_Data_pPb_5023GeV_2040V0A");
				}
			}
		}   
		if ( trainConfig == 46 || trainConfig == 48 || trainConfig == 50 || trainConfig == 52 || trainConfig == 54 || trainConfig == 56 || trainConfig == 58 || trainConfig == 60|| trainConfig == 62||
			trainConfig == 64|| trainConfig == 66){
			if (doWeighting){
				analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_2040V0A",
																			 "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_2040V0A", "","Pi0_Fit_Data_pPb_5023GeV_2040V0A","Eta_Fit_Data_pPb_5023GeV_2040V0A");
			}
		}   
		if ( trainConfig == 67 || trainConfig == 69 || trainConfig == 71 || trainConfig == 73 || trainConfig == 75 || trainConfig == 77 || trainConfig == 79 || trainConfig == 81 || trainConfig == 83 ||
			trainConfig == 85 || trainConfig == 87 ){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_4060V0A",
																				 "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_4060V0A", "", "Pi0_Fit_Data_pPb_5023GeV_4060V0A",
																				 "Eta_Fit_Data_pPb_5023GeV_4060V0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_4060V0A",
																				 "Eta_Hijing_LHC13e7_pPb_5023GeV_4060V0A", "","Pi0_Fit_Data_pPb_5023GeV_4060V0A","Eta_Fit_Data_pPb_5023GeV_4060V0A");
				}
			}
		}   
		if ( trainConfig == 68 || trainConfig == 70 || trainConfig == 72 || trainConfig == 74 || trainConfig == 76 || trainConfig == 78 || trainConfig == 80 || trainConfig == 82 || trainConfig == 84 ||
			trainConfig == 86 || trainConfig == 88){
			if (doWeighting){
				analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_4060V0A",
																			 "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_4060V0A", "","Pi0_Fit_Data_pPb_5023GeV_4060V0A","Eta_Fit_Data_pPb_5023GeV_4060V0A");
			}
		}
		if ( trainConfig == 89 || trainConfig == 91 || trainConfig == 93 || trainConfig == 95 || trainConfig == 97 || trainConfig == 99 || trainConfig == 101 || trainConfig == 103 || 
			trainConfig == 105 || trainConfig == 107 || trainConfig == 109 ){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_6080V0A",
																				 "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_6080V0A", "", "Pi0_Fit_Data_pPb_5023GeV_6080V0A",
																				 "Eta_Fit_Data_pPb_5023GeV_6080V0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_6080V0A",
																				 "Eta_Hijing_LHC13e7_pPb_5023GeV_6080V0A", "", "Pi0_Fit_Data_pPb_5023GeV_6080V0A", "Eta_Fit_Data_pPb_5023GeV_6080V0A");
				}
			}
		}   
		if ( trainConfig == 90 || trainConfig == 92 || trainConfig == 94 || trainConfig == 96 || trainConfig == 98 || trainConfig == 100 || trainConfig == 102 || trainConfig == 104 ||
			trainConfig == 106 || trainConfig == 108 || trainConfig == 108 ){
			if (doWeighting){
				analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_6080V0A",
																			 "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_6080V0A", "", "Pi0_Fit_Data_pPb_5023GeV_6080V0A",
																			 "Eta_Fit_Data_pPb_5023GeV_6080V0A");
			}
		}
		if ( trainConfig == 111 || trainConfig == 113 || trainConfig == 115 || trainConfig == 117 || trainConfig == 119 || trainConfig == 121 || trainConfig == 123 || trainConfig == 125 ||
			trainConfig == 127 || trainConfig == 129 || trainConfig == 131 ){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_60100V0A",
																				 "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_60100V0A", "", "Pi0_Fit_Data_pPb_5023GeV_60100V0A",
																				 "Eta_Fit_Data_pPb_5023GeV_60100V0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_60100V0A",
																				 "Eta_Hijing_LHC13e7_pPb_5023GeV_60100V0A", "","Pi0_Fit_Data_pPb_5023GeV_60100V0A","Eta_Fit_Data_pPb_5023GeV_60100V0A");
				}
			}
		}   
		if ( trainConfig == 112 || trainConfig == 114 || trainConfig == 116 || trainConfig == 118 || trainConfig == 120 || trainConfig == 122 || trainConfig == 124 || trainConfig == 126 || 
			trainConfig == 128 || trainConfig == 130 || trainConfig == 130){
			if (doWeighting){
				analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_60100V0A",
																			"Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_60100V0A", "", "Pi0_Fit_Data_pPb_5023GeV_60100V0A",
																			"Eta_Fit_Data_pPb_5023GeV_60100V0A");
			}
		}
		analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
		if (doEtaShiftIndCuts) {
			analysisEventCuts[i]->DoEtaShift(doEtaShiftIndCuts);
			analysisEventCuts[i]->SetEtaShift(stringShift);
		}
		if (trainConfig == 179 || trainConfig == 180) {
			analysisEventCuts[0]->DoEtaShift(kTRUE);
			analysisEventCuts[0]->SetEtaShift(-0.215);
		}
		
		EventCutList->Add(analysisEventCuts[i]);
		analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
		
		analysisCuts[i] = new AliConversionPhotonCuts();
		analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());	
		if ((trainConfig == 193 || trainConfig == 194 || trainConfig == 195 || trainConfig == 196 )&& i==3 ) {
		        analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);
		}
		analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
		ConvCutList->Add(analysisCuts[i]);
		analysisCuts[i]->SetFillCutHistograms("",kFALSE);
		
		analysisMesonCuts[i] = new AliConversionMesonCuts();
		if (trainConfig == 189 || trainConfig == 190 || trainConfig == 191 || trainConfig == 192) {
		        analysisMesonCuts[i]->SetOpeningAngleCut(0.0);
		}
		analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
		MesonCutList->Add(analysisMesonCuts[i]);
		analysisMesonCuts[i]->SetFillCutHistograms("");
		analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
	}
	if (trainConfig == 189 || trainConfig == 190 || trainConfig == 191 || trainConfig == 192) {
	  task->SetDoTHnSparse(kFALSE);
	}
	task->SetEventCutList(numberOfCuts,EventCutList);
	task->SetConversionCutList(numberOfCuts,ConvCutList);
	task->SetMesonCutList(numberOfCuts,MesonCutList);
	task->SetMoveParticleAccordingToVertex(kTRUE);
	task->SetDoMesonAnalysis(kTRUE);
	task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
	task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
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
