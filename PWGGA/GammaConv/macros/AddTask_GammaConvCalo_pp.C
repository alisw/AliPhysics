void AddTask_GammaConvCalo_pp(  Int_t trainConfig = 1,  //change different set of cuts
								Bool_t isMC = kFALSE, //run MC
								Int_t enableQAMesonTask = 1, //enable QA in AliAnalysisTaskGammaConvV1
								Int_t enableQAPhotonTask = 1, // enable additional QA task
								TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
								TString cutnumberAODBranch = "0000000060084001001500000" 
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
	TString cutnumber = "0000000002084000002200000000"; 
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
	task->SetIsHeavyIon(0);
	task->SetIsMC(isMC);
	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 3;

	TString *cutarray = new TString[numberOfCuts];
	TString *clustercutarray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];

	// meson cuts
	// meson type (Dalitz or not), BG scheme, pool depth, rotation degrees, rapidity cut, radius cut, alpha, chi2, shared electrons, reject to close v0, MC smearing, dca, dca, dca
	
	if (trainConfig == 1){ 
		cutarray[ 0] = "0000001002092970028250400000"; clustercutarray[0] = "10000040022030000"; mesonCutArray[0] = "01525065000000"; //standard cut LHC11h pp 2.76TeV, kMB // EMCAL clusters
		cutarray[ 1] = "0005101002092970028250400000"; clustercutarray[1] = "10000040022030000"; mesonCutArray[1] = "01525065000000"; //standard cut LHC11h pp 2.76TeV, kEMC1 // EMCAL clusters
		cutarray[ 2] = "0002001002092970028250400000"; clustercutarray[2] = "10000040022030000"; mesonCutArray[2] = "01525065000000"; //standard cut LHC11h pp 2.76TeV, SDD V0OR // EMCAL clusters
	} else if (trainConfig == 2){ 
		cutarray[ 0] = "0000001002092970028250400000"; clustercutarray[0] = "20000030022000000"; mesonCutArray[0] = "01525065000000"; //standard cut LHC11h pp 2.76TeV, kMB   // PHOS clusters
		cutarray[ 1] = "0006101002092970028250400000"; clustercutarray[1] = "20000030022000000"; mesonCutArray[1] = "01525065000000"; //standard cut LHC11h pp 2.76TeV, kPHI1 // PHOS clusters
		cutarray[ 2] = "0002001002092970028250400000"; clustercutarray[2] = "20000030022000000"; mesonCutArray[2] = "01525065000000"; //standard cut LHC11h pp 2.76TeV, SDD V0OR //PHOS clusters

		
	} else {
		Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}

	TList *ConvCutList = new TList();
	TList *ClusterCutList = new TList();
	TList *MesonCutList = new TList();

	TList *HeaderList = new TList();
	TObjString *Header1 = new TObjString("BOX");
	HeaderList->Add(Header1);

	ConvCutList->SetOwner(kTRUE);
	AliConversionCuts **analysisCuts = new AliConversionCuts*[numberOfCuts];
	ClusterCutList->SetOwner(kTRUE);
	AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

	for(Int_t i = 0; i<numberOfCuts; i++){
		analysisCuts[i] = new AliConversionCuts();
		analysisCuts[i]->InitializeCutsFromCutString(cutarray[i].Data());
		ConvCutList->Add(analysisCuts[i]);
		analysisCuts[i]->SetFillCutHistograms("",kFALSE);
		
		analysisClusterCuts[i] = new AliCaloPhotonCuts();
		analysisClusterCuts[i]->InitializeCutsFromCutString(clustercutarray[i].Data());
		ClusterCutList->Add(analysisClusterCuts[i]);
		analysisClusterCuts[i]->SetFillCutHistograms("");

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
