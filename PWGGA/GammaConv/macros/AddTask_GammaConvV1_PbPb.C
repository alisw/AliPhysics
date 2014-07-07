void AddTask_GammaConvV1_PbPb(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC 
                              Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvV1
                              Int_t enableQAPhotonTask = 0, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                              Int_t headerSelectionInt = 0,  // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
                              TString cutnumberAODBranch = "1000000060084000001500000",
                              TString periodName = "LHC13d2",  //name of the period for added signals and weighting
							  Bool_t doWeighting = kFALSE  //enable Weighting
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
	TString cutnumberPhoton = "000084001001500000000";
	TString cutnumberEvent = "1000000";
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
	//================================================
	AliAnalysisTaskGammaConvV1 *task=NULL;
	task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 5;

	TString *eventCutArray = new TString[numberOfCuts];
	TString *photonCutArray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	if (trainConfig == 1){ // Standard cuts
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "042092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "042092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "042092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "042092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
	} else if (trainConfig == 2) { // Standard cuts
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "042092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "042092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "042092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
		eventCutArray[ 4] = "5490001"; photonCutArray[ 4] = "042092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
	} else if (trainConfig == 3) { // Standard cuts only added signals
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "042092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "042092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
		eventCutArray[ 3] = "5120002"; photonCutArray[ 3] = "042092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
		eventCutArray[ 4] = "5020002"; photonCutArray[ 4] = "042092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
	} else if (trainConfig == 4) { // Standard cuts only added signals
		eventCutArray[ 0] = "5240002"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
		eventCutArray[ 1] = "5460002"; photonCutArray[ 1] = "042092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
		eventCutArray[ 2] = "5680002"; photonCutArray[ 2] = "042092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
		eventCutArray[ 3] = "5480002"; photonCutArray[ 3] = "042092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 20-40% 
		eventCutArray[ 4] = "5490002"; photonCutArray[ 4] = "042092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
	} else if (trainConfig == 5){ // R-minCut 7.5 cm
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "049092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "049092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "049092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "049092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "049092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
	} else if (trainConfig == 6) { // R-minCut 7.5 cm
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "049092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "049092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "049092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "049092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
		eventCutArray[ 4] = "5490001"; photonCutArray[ 4] = "049092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
	} else if (trainConfig == 7) {// R-minCut 7.5 cm
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "049092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "049092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "049092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
		eventCutArray[ 3] = "5120002"; photonCutArray[ 3] = "049092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
		eventCutArray[ 4] = "5020002"; photonCutArray[ 4] = "049092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
	} else if (trainConfig == 8) { // R-minCut 7.5 cm
		eventCutArray[ 0] = "5240002"; photonCutArray[ 0] = "040092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
		eventCutArray[ 1] = "5460002"; photonCutArray[ 1] = "049092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
		eventCutArray[ 2] = "5680002"; photonCutArray[ 2] = "049092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
		eventCutArray[ 3] = "5480002"; photonCutArray[ 3] = "049092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 20-40% 
		eventCutArray[ 4] = "5490002"; photonCutArray[ 4] = "049092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%   
	} else if (trainConfig == 9){ // R-minCut 12.5 cm
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "048092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "048092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "048092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "048092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "048092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
	} else if (trainConfig == 10) { // R-minCut 12.5 cm
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "048092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "048092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "048092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "048092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
		eventCutArray[ 4] = "5490001"; photonCutArray[ 4] = "048092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
	} else if (trainConfig == 11) {// R-minCut 12.5 cm
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "048092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "048092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "048092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
		eventCutArray[ 3] = "5120002"; photonCutArray[ 3] = "048092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
		eventCutArray[ 4] = "5020002"; photonCutArray[ 4] = "048092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
	} else if (trainConfig == 12) { // R-minCut 12.5 cm
		eventCutArray[ 0] = "5240002"; photonCutArray[ 0] = "040092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
		eventCutArray[ 1] = "5460002"; photonCutArray[ 1] = "048092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
		eventCutArray[ 2] = "5680002"; photonCutArray[ 2] = "048092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
		eventCutArray[ 3] = "5480002"; photonCutArray[ 3] = "048092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 20-40% 
		eventCutArray[ 4] = "5490002"; photonCutArray[ 4] = "048092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%      
	} else  if (trainConfig == 13){ // LHC10h standard, eta 0.65, y = 0.6 
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "032092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "032092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "032092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "032092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "032092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 14) {  // LHC10h standard, eta 0.65, y = 0.6 
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "032092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "032092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "032092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "032092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 40-80%
		eventCutArray[ 4] = "5490001"; photonCutArray[ 4] = "032092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else if (trainConfig == 15) { // LHC10h standard, eta 0.65, y = 0.6  added signals
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "032092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "032092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "032092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120002"; photonCutArray[ 3] = "032092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020002"; photonCutArray[ 4] = "032092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 16) { // LHC10h standard, eta 0.65, y = 0.6  added signals
		eventCutArray[ 0] = "5240002"; photonCutArray[ 0] = "032092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460002"; photonCutArray[ 1] = "032092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680002"; photonCutArray[ 2] = "032092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480002"; photonCutArray[ 3] = "032092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 20-40% 
		eventCutArray[ 4] = "5490002"; photonCutArray[ 4] = "032092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else  if (trainConfig == 17){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 1 
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "032092970023220020000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "032092970023220020000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "032092970023220020000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "032092970023220020000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "032092970023220020000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 18) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 1 
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "032092970023220020000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "032092970023220020000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "032092970023220020000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "032092970023220020000"; mesonCutArray[ 3] = "01523065009000"; // 40-80%
		eventCutArray[ 4] = "5490001"; photonCutArray[ 4] = "032092970023220020000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else if (trainConfig == 19) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 1  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "032092970023220020000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "032092970023220020000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "032092970023220020000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120002"; photonCutArray[ 3] = "032092970023220020000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020002"; photonCutArray[ 4] = "032092970023220020000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 20) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 1 added signal
		eventCutArray[ 0] = "5240002"; photonCutArray[ 0] = "032092970023220020000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460002"; photonCutArray[ 1] = "032092970023220020000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680002"; photonCutArray[ 2] = "032092970023220020000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480002"; photonCutArray[ 3] = "032092970023220020000"; mesonCutArray[ 3] = "01523065009000"; // 20-40% 
		eventCutArray[ 4] = "5490002"; photonCutArray[ 4] = "032092970023220020000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else  if (trainConfig == 21){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 2 
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "032092970023220030000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "032092970023220030000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "032092970023220030000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "032092970023220030000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "032092970023220030000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 22) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 2 
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "032092970023220030000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "032092970023220030000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "032092970023220030000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "032092970023220030000"; mesonCutArray[ 3] = "01523065009000"; // 40-80%
		eventCutArray[ 4] = "5490001"; photonCutArray[ 4] = "032092970023220030000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else if (trainConfig == 23) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 2  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "032092970023220030000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "032092970023220030000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "032092970023220030000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120002"; photonCutArray[ 3] = "032092970023220030000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020002"; photonCutArray[ 4] = "032092970023220030000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 24) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 2 added signal
		eventCutArray[ 0] = "5240002"; photonCutArray[ 0] = "032092970023220030000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460002"; photonCutArray[ 1] = "032092970023220030000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680002"; photonCutArray[ 2] = "032092970023220030000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480002"; photonCutArray[ 3] = "032092970023220030000"; mesonCutArray[ 3] = "01523065009000"; // 20-40% 
		eventCutArray[ 4] = "5490002"; photonCutArray[ 4] = "032092970023220030000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else  if (trainConfig == 25){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 3 
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "032092970023220040000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "032092970023220040000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "032092970023220040000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "032092970023220040000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "032092970023220040000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 26) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 3 
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "032092970023220040000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "032092970023220040000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "032092970023220040000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "032092970023220040000"; mesonCutArray[ 3] = "01523065009000"; // 40-80%
		eventCutArray[ 4] = "5490001"; photonCutArray[ 4] = "032092970023220040000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else if (trainConfig == 27) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 3  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "032092970023220040000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "032092970023220040000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "032092970023220040000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120002"; photonCutArray[ 3] = "032092970023220040000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020002"; photonCutArray[ 4] = "032092970023220040000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 28) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 3 added signal
		eventCutArray[ 0] = "5240002"; photonCutArray[ 0] = "032092970023220040000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460002"; photonCutArray[ 1] = "032092970023220040000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680002"; photonCutArray[ 2] = "032092970023220040000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480002"; photonCutArray[ 3] = "032092970023220040000"; mesonCutArray[ 3] = "01523065009000"; // 20-40% 
		eventCutArray[ 4] = "5490002"; photonCutArray[ 4] = "032092970023220040000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else  if (trainConfig == 29){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm 
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "037092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "037092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "037092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "037092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "037092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 30) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm 
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "037092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "037092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "037092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "037092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 40-80%
		eventCutArray[ 4] = "5490001"; photonCutArray[ 4] = "037092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else if (trainConfig == 31) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm  added signals
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "037092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "037092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "037092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120002"; photonCutArray[ 3] = "037092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020002"; photonCutArray[ 4] = "037092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 32) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm  added signals
		eventCutArray[ 0] = "5240002"; photonCutArray[ 0] = "037092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460002"; photonCutArray[ 1] = "037092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680002"; photonCutArray[ 2] = "037092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480002"; photonCutArray[ 3] = "037092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 20-40% 
		eventCutArray[ 4] = "5490002"; photonCutArray[ 4] = "037092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else  if (trainConfig == 33){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1 
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "037092970023220020000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "037092970023220020000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "037092970023220020000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "037092970023220020000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "037092970023220020000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 34) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1  
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "037092970023220020000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "037092970023220020000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "037092970023220020000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "037092970023220020000"; mesonCutArray[ 3] = "01523065009000"; // 40-80%
		eventCutArray[ 4] = "5490001"; photonCutArray[ 4] = "037092970023220020000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else if (trainConfig == 35) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1  added signals
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "037092970023220020000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "037092970023220020000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "037092970023220020000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120002"; photonCutArray[ 3] = "037092970023220020000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020002"; photonCutArray[ 4] = "037092970023220020000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 36) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1  added signals
		eventCutArray[ 0] = "5240002"; photonCutArray[ 0] = "037092970023220020000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460002"; photonCutArray[ 1] = "037092970023220020000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680002"; photonCutArray[ 2] = "037092970023220020000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480002"; photonCutArray[ 3] = "037092970023220020000"; mesonCutArray[ 3] = "01523065009000"; // 20-40% 
		eventCutArray[ 4] = "5490002"; photonCutArray[ 4] = "037092970023220020000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else  if (trainConfig == 37){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3 
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "037092970023220040000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "037092970023220040000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "037092970023220040000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "037092970023220040000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "037092970023220040000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 38) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3  
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "037092970023220040000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "037092970023220040000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "037092970023220040000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "037092970023220040000"; mesonCutArray[ 3] = "01523065009000"; // 40-80%
		eventCutArray[ 4] = "5490001"; photonCutArray[ 4] = "037092970023220040000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else if (trainConfig == 39) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3  added signals
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "037092970023220040000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "037092970023220040000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "037092970023220040000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
		eventCutArray[ 3] = "5120002"; photonCutArray[ 3] = "037092970023220040000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
		eventCutArray[ 4] = "5020002"; photonCutArray[ 4] = "037092970023220040000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
	} else if (trainConfig == 40) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3  added signals
		eventCutArray[ 0] = "5240002"; photonCutArray[ 0] = "037092970023220040000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
		eventCutArray[ 1] = "5460002"; photonCutArray[ 1] = "037092970023220040000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
		eventCutArray[ 2] = "5680002"; photonCutArray[ 2] = "037092970023220040000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
		eventCutArray[ 3] = "5480002"; photonCutArray[ 3] = "037092970023220040000"; mesonCutArray[ 3] = "01523065009000"; // 20-40% 
		eventCutArray[ 4] = "5490002"; photonCutArray[ 4] = "037092970023220040000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
	} else if (trainConfig == 41){ // Standard cuts, eta 0.9, only to be run on data :kSemiCentral
		eventCutArray[ 0] = "6016001"; photonCutArray[ 0] = "002092970023220000000"; mesonCutArray[ 0] = "01525045009000"; // 0-5%
		eventCutArray[ 1] = "6126001"; photonCutArray[ 1] = "002092970023220000000"; mesonCutArray[ 1] = "01525045009000"; // 5-10%
		eventCutArray[ 2] = "5016001"; photonCutArray[ 2] = "002092970023220000000"; mesonCutArray[ 2] = "01525045009000"; // 0-10%
		eventCutArray[ 3] = "6126001"; photonCutArray[ 3] = "002092970023220000000"; mesonCutArray[ 3] = "01525045009000"; // 10-20%
		eventCutArray[ 4] = "5026001"; photonCutArray[ 4] = "002092970023220000000"; mesonCutArray[ 4] = "01525045009000"; // 0-20%
	} else if (trainConfig == 42) { // Standard cuts, eta 0.9, only to be run on data :kSemiCentral
		eventCutArray[ 0] = "5236001"; photonCutArray[ 0] = "002092970023220000000"; mesonCutArray[ 0] = "01525065009000"; 
		eventCutArray[ 1] = "5346001"; photonCutArray[ 1] = "002092970023220000000"; mesonCutArray[ 1] = "01525065009000"; 
		eventCutArray[ 2] = "5456001"; photonCutArray[ 2] = "002092970023220000000"; mesonCutArray[ 2] = "01525065009000"; 
		eventCutArray[ 3] = "5566001"; photonCutArray[ 3] = "002092970023220000000"; mesonCutArray[ 3] = "01525065009000"; 
		eventCutArray[ 4] = "5676001"; photonCutArray[ 4] = "002092970023220000000"; mesonCutArray[ 4] = "01525065009000"; 
	} else if (trainConfig == 43){ // Standard cuts, eta 0.9, only to be run on data
		eventCutArray[ 0] = "5230001"; photonCutArray[ 0] = "002092970023220000000"; mesonCutArray[ 0] = "01525045009000"; 
		eventCutArray[ 1] = "5340001"; photonCutArray[ 1] = "002092970023220000000"; mesonCutArray[ 1] = "01525065009000"; 
		eventCutArray[ 2] = "5450001"; photonCutArray[ 2] = "002092970023220000000"; mesonCutArray[ 2] = "01525065009000"; 
		eventCutArray[ 3] = "5560001"; photonCutArray[ 3] = "002092970023220000000"; mesonCutArray[ 3] = "01525065009000"; 
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970023220000000"; mesonCutArray[ 4] = "01525065009000"; 
	} else if ( trainConfig == 44){ // qt elipse cut 0.05
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028220000000"; mesonCutArray[ 0] = "01525065009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028220000000"; mesonCutArray[ 1] = "01525065009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028220000000"; mesonCutArray[ 2] = "01525065009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "002092970028220000000"; mesonCutArray[ 3] = "01525065009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "002092970028220000000"; mesonCutArray[ 4] = "01525065009000"; // 0-20%
	} else if ( trainConfig == 45) { // qt elipse cut 0.05
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "002092970028220000000"; mesonCutArray[ 0] = "01525065009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "002092970028220000000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "002092970028220000000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "002092970028220000000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
		eventCutArray[ 4] = "5350001"; photonCutArray[ 4] = "002092970028220000000"; mesonCutArray[ 4] = "01525065009000"; // 30-50%  
	} else if ( trainConfig == 46){ // qt elipse cut 0.05
		eventCutArray[ 0] = "5230001"; photonCutArray[ 0] = "002092970028220000000"; mesonCutArray[ 0] = "01525065009000"; // 20-30% 
		eventCutArray[ 1] = "5340001"; photonCutArray[ 1] = "002092970028220000000"; mesonCutArray[ 1] = "01525065009000"; // 30-40% 
		eventCutArray[ 2] = "5450001"; photonCutArray[ 2] = "002092970028220000000"; mesonCutArray[ 2] = "01525065009000"; // 40-50% 
		eventCutArray[ 3] = "5560001"; photonCutArray[ 3] = "002092970028220000000"; mesonCutArray[ 3] = "01525065009000"; // 50-60% 
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028220000000"; mesonCutArray[ 4] = "01525065009000"; // 60-70% 
	} else if ( trainConfig == 47){ // cos(theta_point) cut 0.85
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970023220400000"; mesonCutArray[ 0] = "01525065009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970023220400000"; mesonCutArray[ 1] = "01525065009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970023220400000"; mesonCutArray[ 2] = "01525065009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "002092970023220400000"; mesonCutArray[ 3] = "01525065009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "002092970023220400000"; mesonCutArray[ 4] = "01525065009000"; // 0-20%
	} else if ( trainConfig == 48) { // cos(theta_point) cut 0.85
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "002092970023220400000"; mesonCutArray[ 0] = "01525065009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "002092970023220400000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "002092970023220400000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "002092970023220400000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
		eventCutArray[ 4] = "5350001"; photonCutArray[ 4] = "002092970023220400000"; mesonCutArray[ 4] = "01525065009000"; // 30-50%  
	} else if ( trainConfig == 49){ // cos(theta_point) cut 0.85
		eventCutArray[ 0] = "5230001"; photonCutArray[ 0] = "002092970023220400000"; mesonCutArray[ 0] = "01525065009000"; // 20-30% 
		eventCutArray[ 1] = "5340001"; photonCutArray[ 1] = "002092970023220400000"; mesonCutArray[ 1] = "01525065009000"; // 30-40% 
		eventCutArray[ 2] = "5450001"; photonCutArray[ 2] = "002092970023220400000"; mesonCutArray[ 2] = "01525065009000"; // 40-50% 
		eventCutArray[ 3] = "5560001"; photonCutArray[ 3] = "002092970023220400000"; mesonCutArray[ 3] = "01525065009000"; // 50-60% 
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970023220400000"; mesonCutArray[ 4] = "01525065009000"; // 60-70% 
	} else if ( trainConfig == 50){ // psi pair 2D 0.05
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970023260000000"; mesonCutArray[ 0] = "01525065009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970023260000000"; mesonCutArray[ 1] = "01525065009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970023260000000"; mesonCutArray[ 2] = "01525065009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "002092970023260000000"; mesonCutArray[ 3] = "01525065009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "002092970023260000000"; mesonCutArray[ 4] = "01525065009000"; // 0-20%
	} else if ( trainConfig == 51) { // psi pair 2D 0.05
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "002092970023260000000"; mesonCutArray[ 0] = "01525065009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "002092970023260000000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "002092970023260000000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "002092970023260000000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
		eventCutArray[ 4] = "5350001"; photonCutArray[ 4] = "002092970023260000000"; mesonCutArray[ 4] = "01525065009000"; // 30-50%  
	} else if ( trainConfig == 52){ // psi pair 2D 0.05
		eventCutArray[ 0] = "5230001"; photonCutArray[ 0] = "002092970023260000000"; mesonCutArray[ 0] = "01525065009000"; // 20-30% 
		eventCutArray[ 1] = "5340001"; photonCutArray[ 1] = "002092970023260000000"; mesonCutArray[ 1] = "01525065009000"; // 30-40% 
		eventCutArray[ 2] = "5450001"; photonCutArray[ 2] = "002092970023260000000"; mesonCutArray[ 2] = "01525065009000"; // 40-50% 
		eventCutArray[ 3] = "5560001"; photonCutArray[ 3] = "002092970023260000000"; mesonCutArray[ 3] = "01525065009000"; // 50-60% 
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970023260000000"; mesonCutArray[ 4] = "01525065009000"; // 60-70%       
	} else if ( trainConfig == 53){ // psi pair 2D 0.1
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970023250000000"; mesonCutArray[ 0] = "01525065009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970023250000000"; mesonCutArray[ 1] = "01525065009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970023250000000"; mesonCutArray[ 2] = "01525065009000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "002092970023250000000"; mesonCutArray[ 3] = "01525065009000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "002092970023250000000"; mesonCutArray[ 4] = "01525065009000"; // 0-20%
	} else if ( trainConfig == 54) { // psi pair 2D 0.1
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "002092970023250000000"; mesonCutArray[ 0] = "01525065009000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "002092970023250000000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "002092970023250000000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "002092970023250000000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
		eventCutArray[ 4] = "5350001"; photonCutArray[ 4] = "002092970023250000000"; mesonCutArray[ 4] = "01525065009000"; // 30-50%  
	} else if ( trainConfig == 55){ // psi pair 2D 0.1
		eventCutArray[ 0] = "5230001"; photonCutArray[ 0] = "002092970023250000000"; mesonCutArray[ 0] = "01525065009000"; // 20-30% 
		eventCutArray[ 1] = "5340001"; photonCutArray[ 1] = "002092970023250000000"; mesonCutArray[ 1] = "01525065009000"; // 30-40% 
		eventCutArray[ 2] = "5450001"; photonCutArray[ 2] = "002092970023250000000"; mesonCutArray[ 2] = "01525065009000"; // 40-50% 
		eventCutArray[ 3] = "5560001"; photonCutArray[ 3] = "002092970023250000000"; mesonCutArray[ 3] = "01525065009000"; // 50-60% 
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970023250000000"; mesonCutArray[ 4] = "01525065009000"; // 60-70%          
	} else if ( trainConfig == 56){ // cleaner cuts
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5120001"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5020001"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 57){ // cleaner cuts added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5120002"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5020002"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%   
	} else if ( trainConfig == 58) { // cleaner cuts
		eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 20-40%
		eventCutArray[ 1] = "5460001"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 40-60%
		eventCutArray[ 2] = "5680001"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 60-80%
		eventCutArray[ 3] = "5480001"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 40-80%
		eventCutArray[ 4] = "5350001"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 30-50%  
	} else if ( trainConfig == 59) { // cleaner cuts added signal
		eventCutArray[ 0] = "5240002"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 20-40%
		eventCutArray[ 1] = "5460002"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 40-60%
		eventCutArray[ 2] = "5680002"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 60-80%
		eventCutArray[ 3] = "5480002"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 40-80%
		eventCutArray[ 4] = "5350002"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 30-50%  		
	} else if ( trainConfig == 60){ // cleaner cuts
		eventCutArray[ 0] = "5230001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 20-30% 
		eventCutArray[ 1] = "5340001"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 30-40% 
		eventCutArray[ 2] = "5450001"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 40-50% 
		eventCutArray[ 3] = "5560001"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 50-60% 
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 60-70%                
	} else if ( trainConfig == 61){ // cleaner cuts added signal
		eventCutArray[ 0] = "5230002"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 20-30% 
		eventCutArray[ 1] = "5340002"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 30-40% 
		eventCutArray[ 2] = "5450002"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 40-50% 
		eventCutArray[ 3] = "5560002"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 50-60% 
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 60-70%                   
	} else if ( trainConfig == 62){ // cleaner cuts
		eventCutArray[ 0] = "6230001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6340001"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "6450001"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "6560001"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "6670001"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 63){ // cleaner cuts added signal
		eventCutArray[ 0] = "6230002"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6340002"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "6450002"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "6560002"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "6670002"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 64){ // cleaner cuts
		eventCutArray[ 0] = "6780001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6890001"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5670001"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5780001"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5890001"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 65){ // cleaner cuts added signal
		eventCutArray[ 0] = "6780002"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6890002"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5670002"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5780002"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5890002"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 66){ // cleaner cuts
		eventCutArray[ 0] = "7010001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "7120001"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "7230001"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "7340001"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "7450001"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 67){ // cleaner cuts added signal
		eventCutArray[ 0] = "7010002"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "7120002"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "7230002"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "7340002"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "7450002"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 68){ // cleaner cuts
		eventCutArray[ 0] = "7560001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "7670001"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "7780001"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "7890001"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 10-20%
		eventCutArray[ 4] = "7090001"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 0-20%	
	} else if ( trainConfig == 69){ // cleaner cuts added signal
		eventCutArray[ 0] = "7560002"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "7670002"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "7780002"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "7890002"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 10-20%
		eventCutArray[ 4] = "7090002"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 0-20%			
	} else if ( trainConfig == 70){ // variation eta  0.65
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "032092970028250400000"; mesonCutArray[ 0]= "01523065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "032092970028250400000"; mesonCutArray[ 1]= "01523065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "032092970028250400000"; mesonCutArray[ 2]= "01523065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "032092970028250400000"; mesonCutArray[ 3]= "01523065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "032092970028250400000"; mesonCutArray[ 4]= "01523065000000"; // 20-50% 
	} else if ( trainConfig == 71){ // variation eta  0.65 added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "032092970028250400000"; mesonCutArray[ 0]= "01523065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "032092970028250400000"; mesonCutArray[ 1]= "01523065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "032092970028250400000"; mesonCutArray[ 2]= "01523065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "032092970028250400000"; mesonCutArray[ 3]= "01523065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "032092970028250400000"; mesonCutArray[ 4]= "01523065000000"; // 20-50% 		
	} else if ( trainConfig == 72){ // variation eta  0.75
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "042092970028250400000"; mesonCutArray[ 0]= "01522065009000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "042092970028250400000"; mesonCutArray[ 1]= "01522065009000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "042092970028250400000"; mesonCutArray[ 2]= "01522065009000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "042092970028250400000"; mesonCutArray[ 3]= "01522065009000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "042092970028250400000"; mesonCutArray[ 4]= "01522065009000"; // 20-50% 
	} else if ( trainConfig == 73){ // variation eta  0.75 added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "042092970028250400000"; mesonCutArray[ 0]= "01522065009000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "042092970028250400000"; mesonCutArray[ 1]= "01522065009000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "042092970028250400000"; mesonCutArray[ 2]= "01522065009000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "042092970028250400000"; mesonCutArray[ 3]= "01522065009000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "042092970028250400000"; mesonCutArray[ 4]= "01522065009000"; // 20-50% 
	} else if ( trainConfig == 74){ // single pt 0.075
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002492970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002492970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002492970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002492970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002492970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 75){ // single pt 0.075 added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002492970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002492970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002492970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002492970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002492970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 76){ // single pt 0.1
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002192970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002192970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002192970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002192970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002192970028250400000"; mesonCutArray[ 4]= "01525065000000"; //20-50%
	} else if ( trainConfig == 77){ // single pt 0.1 added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002192970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002192970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002192970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002192970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002192970028250400000"; mesonCutArray[ 4]= "01525065000000"; //20-50%
	} else if ( trainConfig == 78){ // variation TPC cls 0.7
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002062970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002062970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002062970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002062970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002062970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 79){ // variation TPC cls 0.7 added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002062970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002062970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002062970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002062970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002062970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 80){ // variation TPC cls 0.35
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002082970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002082970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002082970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002082970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002082970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 81){ // variation TPC cls 0.35 added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002082970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002082970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002082970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002082970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002082970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 82){ // variation edEdx  -4,5
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002093970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002093970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002093970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002093970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002093970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 83){ // variation edEdx  -4,5  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002093970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002093970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002093970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002093970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002093970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 84){ // variation edEdx  -2.5,4
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002096970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002096970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002096970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002096970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002096970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 85){ // variation edEdx  -2.5,4  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002096970028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002096970028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002096970028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002096970028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002096970028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 86){ //variation pion p dEdx 0.3-5.
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092951028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092951028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092951028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092951028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092951028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 87){ //variation pion p dEdx 0.3-5.  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092951028250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092951028250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092951028250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092951028250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092951028250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 88){ // TOF el. PID -3,5
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970038250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970038250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970038250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970038250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970038250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 89){ // TOF el. PID -3,5  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970038250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970038250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970038250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970038250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970038250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%		
	} else if ( trainConfig == 90){ // TOF el. PID -2,3
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970048250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970048250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970048250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970048250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970048250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 91){ // TOF el. PID -2,3  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970048250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970048250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970048250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970048250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970048250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 92){ // qt 0.03
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970029250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970029250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970029250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970029250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970029250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 93){ // qt 0.03  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970029250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970029250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970029250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970029250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970029250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 94){ // qt 0.07 no2D
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970022250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970022250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970022250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970022250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970022250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig 	== 95){ // qt 0.07 no2D  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970022250400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970022250400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970022250400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970022250400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970022250400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 96){ // chi2  50.
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028150400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028150400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028150400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970028150400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028150400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 97){ // chi2  50.  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028150400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028150400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028150400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970028150400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028150400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%		
	} else if ( trainConfig == 98){ // chi2  20.
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028850400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028850400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028850400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970028850400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028850400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 99){ // chi2  20.  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028850400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028850400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028850400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970028850400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028850400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 100){ // psi pair 0.05
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028260400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028260400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028260400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970028260400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028260400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 101){ // psi pair 0.05  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028260400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028260400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028260400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970028260400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028260400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 102){ // cosPA -1
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028250000000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028250000000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028250000000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970028250000000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028250000000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 103){ // cosPA -1  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028250000000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028250000000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028250000000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970028250000000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028250000000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 104){ // variation alpha 0.75
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0]= "01525055000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1]= "01525055000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2]= "01525055000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3]= "01525055000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4]= "01525055000000"; // 20-50%
	} else if ( trainConfig == 105){ // variation alpha 0.75  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0]= "01525055000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1]= "01525055000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2]= "01525055000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3]= "01525055000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4]= "01525055000000"; // 20-50%	
	} else if ( trainConfig == 106){ // variation alpha 0.85
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0]= "01525075000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1]= "01525075000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2]= "01525075000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3]= "01525075000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4]= "01525075000000"; // 20-50%
	} else if ( trainConfig == 107){ // variation alpha 0.85  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0]= "01525075000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1]= "01525075000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2]= "01525075000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3]= "01525075000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4]= "01525075000000"; // 20-50%
	} else if ( trainConfig == 108){ // psi pair 0.2
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028280400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028280400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028280400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970028280400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028280400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 109){ // psi pair 0.2  added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028280400000"; mesonCutArray[ 0]= "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028280400000"; mesonCutArray[ 1]= "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028280400000"; mesonCutArray[ 2]= "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970028280400000"; mesonCutArray[ 3]= "01525065000000"; // 20-40%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028280400000"; mesonCutArray[ 4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 110){ // cleaner cuts
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 111){ // cleaner cuts added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%   		
	} else if ( trainConfig == 112){ // cleaner cuts photon Quality 1
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028250420000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028250420000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028250420000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970028250420000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028250420000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 113){ // cleaner cuts added signal photon Quality 1
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028250420000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028250420000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028250420000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970028250420000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028250420000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%   		
	} else if ( trainConfig == 114){ // cleaner cuts photon Quality 2
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028250430000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028250430000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028250430000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970028250430000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028250430000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 115){ // cleaner cuts added signal photon Quality 2
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028250430000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028250430000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028250430000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970028250430000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028250430000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%   		
	} else if ( trainConfig == 116){ // cleaner cuts photon Quality 3
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028250440000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "002092970028250440000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "002092970028250440000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "002092970028250440000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "002092970028250440000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 117){ // cleaner cuts added signal photon Quality 3
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "002092970028250440000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "002092970028250440000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "002092970028250440000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "002092970028250440000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "002092970028250440000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%   		
	} else if ( trainConfig == 118){ // cleaner cuts, min R = 35 cm
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "007092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "007092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "007092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "007092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "007092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 119){ // cleaner cuts, min R = 35 cm added signal
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "007092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "007092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "007092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "007092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "007092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%   		
	} else if ( trainConfig == 120){ // cleaner cuts, photon Quality 1, min R = 35 cm
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "007092970028250420000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "007092970028250420000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "007092970028250420000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "007092970028250420000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "007092970028250420000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 121){ // cleaner cuts added signal, photon Quality 1, min R = 35 cm
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "007092970028250420000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "007092970028250420000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "007092970028250420000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "007092970028250420000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "007092970028250420000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%   		
	} else if ( trainConfig == 122){ // cleaner cuts, photon Quality 3, min R = 35 cm
		eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "007092970028250440000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120001"; photonCutArray[ 1] = "007092970028250440000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010001"; photonCutArray[ 2] = "007092970028250440000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240001"; photonCutArray[ 3] = "007092970028250440000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250001"; photonCutArray[ 4] = "007092970028250440000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 123){ // cleaner cuts added signal, photon Quality 3, min R = 35 cm
		eventCutArray[ 0] = "6010002"; photonCutArray[ 0] = "007092970028250440000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		eventCutArray[ 1] = "6120002"; photonCutArray[ 1] = "007092970028250440000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		eventCutArray[ 2] = "5010002"; photonCutArray[ 2] = "007092970028250440000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		eventCutArray[ 3] = "5240002"; photonCutArray[ 3] = "007092970028250440000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
		eventCutArray[ 4] = "5250002"; photonCutArray[ 4] = "007092970028250440000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%   		
	} else {
		Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *EventCutList = new TList();
	TList *ConvCutList = new TList();
	TList *MesonCutList = new TList();

	TList *HeaderList = new TList();
	if (periodName.CompareTo("LHC13d2")==0){
		TObjString *Header1 = new TObjString("pi0_1");
		HeaderList->Add(Header1);
	//    TObjString *Header3 = new TObjString("eta_2");
	//    HeaderList->Add(Header3);

	} else if (periodName.CompareTo("LHC12a17x_fix")==0){
		TObjString *Header1 = new TObjString("PARAM");
		HeaderList->Add(Header1);
	} else if (periodName.CompareTo("LHC14a1a")==0){
		if (headerSelectionInt == 1){ 
			TObjString *Header1 = new TObjString("pi0_1");
			HeaderList->Add(Header1);
		} else if (headerSelectionInt == 2){
			TObjString *Header1 = new TObjString("eta_2");
			HeaderList->Add(Header1);
		} else {
			TObjString *Header1 = new TObjString("pi0_1");
			HeaderList->Add(Header1);
			TObjString *Header2 = new TObjString("eta_2");
			HeaderList->Add(Header2);
		}  
	} else if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0){
		TObjString *Header1 = new TObjString("BOX");
		HeaderList->Add(Header1);
	}	

	EventCutList->SetOwner(kTRUE);
	AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
	ConvCutList->SetOwner(kTRUE);
	AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

	for(Int_t i = 0; i<numberOfCuts; i++){
		analysisEventCuts[i] = new AliConvEventCuts();
		if (trainConfig == 1 ||trainConfig == 5 || trainConfig == 9 || trainConfig == 13 || trainConfig == 17 || trainConfig == 21 || trainConfig == 25 || trainConfig == 29 || trainConfig == 33 || trainConfig == 37){ // || trainConfig == 41 
			if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
			if (i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
			if (i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
			if (i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
			if (i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
		} else if (trainConfig == 2 ||trainConfig == 6 || trainConfig == 10 || trainConfig == 14 || trainConfig == 18 || trainConfig == 22 || trainConfig == 26 || trainConfig == 30  || trainConfig == 34 || trainConfig == 38){ // || trainConfig == 42
			if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
			if (i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
			if (i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
			if (i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
		} else if (trainConfig == 3 ||trainConfig == 7 || trainConfig == 11 || trainConfig == 15 || trainConfig == 19 || trainConfig == 23 || trainConfig == 27 || trainConfig == 31 || trainConfig == 35 || trainConfig == 39 ){ //|| trainConfig == 43 
			if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
			if (i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
			if (i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
			if (i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
			if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
		} else if (trainConfig == 4 ||trainConfig == 8 || trainConfig == 12 || trainConfig == 16 || trainConfig == 20 || trainConfig == 24 || trainConfig == 28 || trainConfig == 32 || trainConfig == 36 || trainConfig == 40){ // || trainConfig == 44 
			if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
			if (i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
			if (i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
			if (i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
		}

		if (trainConfig == 56 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
				if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
				if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
				if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_1020TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_1020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
				if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0020TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M","Eta_Fit_Data_PbPb_2760GeV_0020V0M");
			}	
		}	  
		if (trainConfig == 57 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
				if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
				if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
				if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
				if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0020TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M","Eta_Fit_Data_PbPb_2760GeV_0020V0M");
			}	
		}	  
		if (trainConfig == 58 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
				if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4060TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M","Eta_Fit_Data_PbPb_2760GeV_4060V0M");
				if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_6080TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_6080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M","Eta_Fit_Data_PbPb_2760GeV_6080V0M");
				if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4080TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
				if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_3050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
			}	
		}	  
		if (trainConfig == 59 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
				if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4060TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M","Eta_Fit_Data_PbPb_2760GeV_4060V0M");
				if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_6080TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_6080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M","Eta_Fit_Data_PbPb_2760GeV_6080V0M");
				if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4080TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
				if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
			}	
		}	  
	
		if (trainConfig == 60 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2030TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2030TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
				if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_3040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
				if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4050V0M","Eta_Fit_Data_PbPb_2760GeV_4050V0M");
				if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_5060TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_5060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_5060V0M","Eta_Fit_Data_PbPb_2760GeV_5060V0M");
				if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
			}	
		}	  
		if (trainConfig == 61 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
				if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
				if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4050V0M","Eta_Fit_Data_PbPb_2760GeV_4050V0M");
				if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_5060TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_5060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_5060V0M","Eta_Fit_Data_PbPb_2760GeV_5060V0M");
				if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
			}	
		}	  
		
		if ( trainConfig == 70 || trainConfig == 72  || trainConfig == 74  || trainConfig == 76  || trainConfig == 78  || trainConfig == 80  || trainConfig == 82  || trainConfig == 84 || trainConfig == 86  || trainConfig == 88  || trainConfig == 90 || trainConfig == 92 || trainConfig == 94 || trainConfig == 96  || trainConfig == 98  || trainConfig == 100 || trainConfig == 102  || trainConfig == 104 || trainConfig == 106 || trainConfig == 108 || trainConfig == 110 || trainConfig == 112 || trainConfig == 114 || trainConfig == 116 || trainConfig == 118 || trainConfig == 120 || trainConfig == 122){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
				if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
				if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
				if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
				if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
			}	
		} 
		if ( trainConfig == 71 || trainConfig == 73  || trainConfig == 75  || trainConfig == 77  || trainConfig == 79  || trainConfig == 81  || trainConfig == 83  || trainConfig == 85 || trainConfig == 87  || trainConfig == 89  || trainConfig == 91 || trainConfig == 93 || trainConfig == 95 || trainConfig == 97  || trainConfig == 99  || trainConfig == 101 || trainConfig == 103  || trainConfig == 105 || trainConfig == 107 || trainConfig == 109 || trainConfig == 111 || trainConfig == 113 || trainConfig == 115 || trainConfig == 117 || trainConfig == 119 || trainConfig == 121 || trainConfig == 123){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
				if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
				if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
				if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
				if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
			}	
		} 
		
		analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
		if (periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
			if (headerSelectionInt == 1) analysisEventCuts[i]->SetAddedSignalPDGCode(111);
			if (headerSelectionInt == 2) analysisEventCuts[i]->SetAddedSignalPDGCode(221);
		}
		EventCutList->Add(analysisEventCuts[i]);		
		if (trainConfig == 37 || trainConfig == 38){
			analysisEventCuts[i]->SelectSpecialTrigger(AliVEvent::kMB, "AliVEvent::kMB" );
		}
		if (trainConfig == 39 || trainConfig == 40){   
			analysisEventCuts[i]->SelectSpecialTrigger(AliVEvent::kCentral,"AliVEvent::kCentral" );
		}   
		if (trainConfig == 41 || trainConfig == 42){   
			analysisEventCuts[i]->SelectSpecialTrigger(AliVEvent::kSemiCentral,"AliVEvent::kSemiCentral" );
		}   
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
