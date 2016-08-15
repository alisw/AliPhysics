AliAnalysisTask *AddTask_ConversionAODProduction( Int_t dataset                 = 0, 
                                                  Bool_t isMC                   = kFALSE, 
                                                  TString periodNameV0Reader    = "",
						  Bool_t addv0sInESDFilter 	= kTRUE
                                                ){

	// Before doing anything, we load the needed library
	gSystem->Load("libPWGGAGammaConv");
	// dataset 0: pp
	// dataset 1: PbPb
	// dataset 2: pPb

	//get the current analysis manager
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTask_V0ReaderV1", "No analysis manager found.");
		return 0;
	}

//========= Add PID Reponse to ANALYSIS manager ====
	if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
	AddTaskPIDResponse(isMC);
	}
	
	TString analysiscut;
	TString analysiscutEvent;
	TString analysiscutB;

	if(dataset == 1){
		analysiscutEvent = "10000003";
        analysiscut= "06000008400000001000000000";
		analysiscutB="16000008400000001000000000";
	} else if (dataset == 2){
		analysiscutEvent = "80000003";
		analysiscut= "06000008400000001000000000";
		analysiscutB="16000008400000001000000000";
	} else{
		analysiscutEvent = "00000003";
		analysiscut ="06000008400100001000000000";
		analysiscutB="16000008400100001000000000";
	}

	//========= Add V0 Reader to  ANALYSIS manager =====

	AliV0ReaderV1 *fV0Reader=new AliV0ReaderV1("ConvGammaAODProduction");
	if (periodNameV0Reader.CompareTo("") != 0) fV0Reader->SetPeriodName(periodNameV0Reader);
    fV0Reader->SetCreateAODs(kTRUE);
	fV0Reader->SetUseOwnXYZCalculation(kTRUE);
	fV0Reader->SetUseAODConversionPhoton(kTRUE);
	if (addv0sInESDFilter){fV0Reader->SetAddv0sInESDFilter(kTRUE);}
//     fV0Reader->CheckAODConsistency();

	AliV0ReaderV1 *fV0ReaderB=new AliV0ReaderV1("ConvGammaAODProductionB");
    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderB->SetPeriodName(periodNameV0Reader);
	fV0ReaderB->SetCreateAODs(kTRUE);
	fV0ReaderB->SetUseOwnXYZCalculation(kTRUE);
	fV0ReaderB->SetUseAODConversionPhoton(kTRUE);
	if (addv0sInESDFilter){fV0ReaderB->SetAddv0sInESDFilter(kTRUE);}
//     fV0ReaderB->CheckAODConsistency();

	AliConvEventCuts *fEventCutsA=NULL;
	AliConvEventCuts *fEventCutsB=NULL;
	if(analysiscutEvent!=""){
		fEventCutsA= new AliConvEventCuts(analysiscutEvent.Data(),analysiscutEvent.Data());
		fEventCutsA->SetPreSelectionCutFlag(kTRUE);
		fEventCutsA->SetV0ReaderName("ConvGammaAODProduction");
		if(fEventCutsA->InitializeCutsFromCutString(analysiscutEvent.Data())){
			fV0Reader->SetEventCuts(fEventCutsA);
		}
		fEventCutsB= new AliConvEventCuts(analysiscutEvent.Data(),analysiscutEvent.Data());
		fEventCutsB->SetPreSelectionCutFlag(kTRUE);
		fEventCutsB->SetV0ReaderName("ConvGammaAODProductionB");
		if(fEventCutsB->InitializeCutsFromCutString(analysiscutEvent.Data())){
			fV0ReaderB->SetEventCuts(fEventCutsB);
		}
	}
	
	// Set AnalysisCut Number
	AliConversionPhotonCuts *fCuts= new AliConversionPhotonCuts(analysiscut.Data(),analysiscut.Data());
	AliConversionPhotonCuts *fCutsB= new AliConversionPhotonCuts(analysiscutB.Data(),analysiscutB.Data());
	if(fCuts->InitializeCutsFromCutString(analysiscut.Data())){
		fCuts->SetIsHeavyIon(dataset);
		fV0Reader->SetConversionCuts(fCuts);
	}
	
	fV0Reader->Init();

	if(fCutsB->InitializeCutsFromCutString(analysiscutB.Data())){
		fCutsB->SetIsHeavyIon(dataset);
		fV0ReaderB->SetConversionCuts(fCutsB);
	}
	fV0ReaderB->Init();

	AliLog::SetGlobalLogLevel(AliLog::kInfo);

	//================================================
	//              data containers
	//================================================
	//            find input container
	//below the trunk version
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutputPCMv0sA = mgr->CreateContainer("PCM offlineV0Finder container", TBits::Class(),AliAnalysisManager::kExchangeContainer);
	AliAnalysisDataContainer *coutputPCMv0sB = mgr->CreateContainer("PCM onflyV0Finder container", TBits::Class(),AliAnalysisManager::kExchangeContainer);
	mgr->ConnectOutput(fV0Reader,1,coutputPCMv0sA);
	mgr->ConnectOutput(fV0ReaderB,1,coutputPCMv0sB);

	// connect input V0Reader
	mgr->AddTask(fV0Reader);
	mgr->ConnectInput (fV0Reader,0,cinput);

	mgr->AddTask(fV0ReaderB);
	mgr->ConnectInput (fV0ReaderB,0,cinput);

	return fV0Reader;
}
