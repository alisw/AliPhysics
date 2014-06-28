void AddTask_GammaConvCalo_PbPb(   	Int_t trainConfig = 1,  //change different set of cuts
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
		AddTaskPIDResponse(isMC,1,0,4,0,"",1,1,4);
	}
	
	//=========  Set Cutnumber for V0Reader ================================
	TString cutnumber = "1000000000084001001500000000"; 
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

		// Set AnalysisCut Number
		AliConversionCuts *fCuts=NULL;
		if(cutnumber!=""){
			fCuts= new AliConversionCuts(cutnumber.Data(),cutnumber.Data());
			fCuts->SetPreSelectionCutFlag(kTRUE);
			if(fCuts->InitializeCutsFromCutString(cutnumber.Data())){
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
	task->SetIsHeavyIon(1);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 5;

	TString *cutarray = new TString[numberOfCuts];
	TString *clustercutarray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// meson cuts
	// meson type (Dalitz or not), BG scheme, pool depth, rotation degrees, rapidity cut, radius cut, alpha, chi2, shared electrons, reject to close v0, MC smearing, dca, dca, dca
  
	if (trainConfig == 1){ 
		cutarray[ 0] = "6010001002092970028250400000"; clustercutarray[0] = "10000040022030000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002092970028250400000"; clustercutarray[1] = "10000040022030000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002092970028250400000"; clustercutarray[2] = "10000040022030000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002092970028250400000"; clustercutarray[3] = "10000040022030000"; mesonCutArray[ 3] = "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002092970028250400000"; clustercutarray[4] = "10000040022030000"; mesonCutArray[ 4] = "01525065000000"; // 20-50%
	} else {
		Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

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

	ConvCutList->SetOwner(kTRUE);
	AliConversionCuts **analysisCuts = new AliConversionCuts*[numberOfCuts];
   	ClusterCutList->SetOwner(kTRUE);
	AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

	for(Int_t i = 0; i<numberOfCuts; i++){
		analysisCuts[i] = new AliConversionCuts();
      
		if ( trainConfig == 1){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
			}	
		} 
		analysisCuts[i]->InitializeCutsFromCutString(cutarray[i].Data());
		if (periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
			if (headerSelectionInt == 1) analysisCuts[i]->SetAddedSignalPDGCode(111);
			if (headerSelectionInt == 2) analysisCuts[i]->SetAddedSignalPDGCode(221);
		}
		ConvCutList->Add(analysisCuts[i]);

		analysisClusterCuts[i] = new AliCaloPhotonCuts();
		analysisClusterCuts[i]->InitializeCutsFromCutString(clustercutarray[i].Data());
		ClusterCutList->Add(analysisClusterCuts[i]);
		analysisClusterCuts[i]->SetFillCutHistograms("");

		analysisCuts[i]->SetFillCutHistograms("",kFALSE);
		analysisMesonCuts[i] = new AliConversionMesonCuts();
		analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
		MesonCutList->Add(analysisMesonCuts[i]);
		analysisMesonCuts[i]->SetFillCutHistograms("");
		analysisCuts[i]->SetAcceptedHeader(HeaderList);
	}

	task->SetConversionCutList(numberOfCuts,ConvCutList);
	task->SetCaloCutList(numberOfCuts,ClusterCutList);
	task->SetMesonCutList(numberOfCuts,MesonCutList);
	task->SetMoveParticleAccordingToVertex(kTRUE);
	task->SetDoMesonAnalysis(kTRUE);
	task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
	task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA

	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer(Form("GammaConvCalo_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaConvCalo_%i.root",trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
