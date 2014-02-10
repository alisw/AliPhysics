void AddTask_GammaConvEtaPiPlPiMiGamma_pPb(    
										Int_t trainConfig = 1,
										Bool_t isMC       = kFALSE, //run MC 
										Bool_t enableQAMesonTask = kTRUE, //enable QA in AliAnalysisTaskEtaToPiPlPiMiGamma
										TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
										Bool_t doWeighting = kFALSE,  //enable Weighting
										TString generatorName = "HIJING",				
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
		Error(Form("AddTask_GammaConvEtaPiPlPiMiGamma_pPb_%i",trainConfig), "No analysis manager found.");
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
	TString ConvCutnumber = "8000000060084001001500000000";   //Online  V0 finder
	TString PionCuts      = "00000221";            //Electron Cuts
		
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

		// Set AnalysisCut Number
		AliConversionCuts *fCuts=NULL;
		if( ConvCutnumber !=""){
			fCuts= new AliConversionCuts(ConvCutnumber.Data(),ConvCutnumber.Data());
			fCuts->SetPreSelectionCutFlag(kTRUE);
			if(fCuts->InitializeCutsFromCutString(ConvCutnumber.Data())){
				fCuts->DoEtaShift(doEtaShift);
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

	
	
	AliAnalysisTaskEtaToPiPlPiMiGamma *task=NULL;

	task= new AliAnalysisTaskEtaToPiPlPiMiGamma(Form("GammaConvEtaPiPlPiMiGamma_%i",trainConfig));

	task->SetIsHeavyIon(2);
	task->SetIsMC(isMC);

	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 1;
	TString *ConvCutarray    = new TString[numberOfCuts];
	TString *PionCutarray    = new TString[numberOfCuts];

	TString *MesonCutarray   = new TString[numberOfCuts];

	Bool_t doEtaShiftIndCuts = kFALSE;
	TString stringShift = "";

	// Shifting in pPb direction

	doEtaShiftIndCuts = kTRUE;
	stringShift = "pPb";

	if( trainConfig == 1 ) {
		ConvCutarray[0] = "8000011002092170008260400000"; PionCutarray[0] = "00000362"; MesonCutarray[0] = "01039035009000"; //standard cut Pi0 PbPb 00-100			
	} 
	
	TList *ConvCutList  = new TList();
	TList *MesonCutList = new TList();
	TList *PionCutList  = new TList();

	TList *HeaderList = new TList();
	TObjString *Header1 = new TObjString("pi0_1");
	HeaderList->Add(Header1);
	TObjString *Header3 = new TObjString("eta_2");
	HeaderList->Add(Header3);
	
	ConvCutList->SetOwner(kTRUE);
	AliConversionCuts **analysisCuts             = new AliConversionCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
	PionCutList->SetOwner(kTRUE);
	AliPrimaryPionCuts **analysisPionCuts     = new AliPrimaryPionCuts*[numberOfCuts];

	for(Int_t i = 0; i<numberOfCuts; i++){

		analysisCuts[i] = new AliConversionCuts();
		if( ! analysisCuts[i]->InitializeCutsFromCutString(ConvCutarray[i].Data()) ) {
				cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
				return 0;
		} else {				
			ConvCutList->Add(analysisCuts[i]);
			analysisCuts[i]->SetFillCutHistograms("",kFALSE);
			analysisCuts[i]->SetAcceptedHeader(HeaderList);
		}

		analysisMesonCuts[i] = new AliConversionMesonCuts();
		
		if( ! analysisMesonCuts[i]->InitializeCutsFromCutString(MesonCutarray[i].Data()) ) {
			cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
			return 0;
		} else {
			MesonCutList->Add(analysisMesonCuts[i]);
			analysisMesonCuts[i]->SetFillCutHistograms("");
		}


		TString cutName( Form("%s_%s_%s",ConvCutarray[i].Data(),PionCutarray[i].Data(),MesonCutarray[i].Data() ) );
		analysisPionCuts[i] = new AliPrimaryPionCuts();
		if( !analysisPionCuts[i]->InitializeCutsFromCutString(PionCutarray[i].Data())) {
			cout<< "ERROR:  analysisPionCuts [ " <<i<<" ] "<<endl;
			return 0;
		} else { 
			PionCutList->Add(analysisPionCuts[i]);
			analysisPionCuts[i]->SetFillCutHistograms("",kFALSE,cutName); 
		}
		

	}


	task->SetConversionCutList(numberOfCuts,ConvCutList);
	task->SetMesonCutList(MesonCutList);
	task->SetPionCutList(PionCutList);

	task->SetMoveParticleAccordingToVertex(kTRUE);

	if(enableQAMesonTask) task->SetDoMesonQA(kTRUE);

	//connect containers
	AliAnalysisDataContainer *coutput =
	mgr->CreateContainer(Form("GammaConvEtaPiPlPiMiGamma_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaConvEtaPiPlPiMiGamma_%i.root",trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
