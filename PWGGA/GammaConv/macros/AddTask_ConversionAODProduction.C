AliAnalysisTask *AddTask_ConversionAODProduction(Int_t dataset=0, Bool_t isMC = kFALSE){

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
		analysiscutEvent = "1000000";
		analysiscut= "060000084000001500000000";
		analysiscutB="160000084000001500000000";
	} else if (dataset == 2){
		analysiscutEvent = "8000000";
		analysiscut= "060000084000001500000000";
		analysiscutB="160000084000001500000000";
	} else{
		analysiscutEvent = "0000000";
		analysiscut ="060000084001001500000000";
		analysiscutB="160000084001001500000000";
	}

	//========= Add V0 Reader to  ANALYSIS manager =====

	AliV0ReaderV1 *fV0Reader=new AliV0ReaderV1("ConvGammaAODProduction");
	fV0Reader->SetCreateAODs(kTRUE);
	fV0Reader->SetUseOwnXYZCalculation(kTRUE);
	fV0Reader->SetUseAODConversionPhoton(kTRUE);
//     fV0Reader->CheckAODConsistency();

	AliV0ReaderV1 *fV0ReaderB=new AliV0ReaderV1("ConvGammaAODProductionB");
	fV0ReaderB->SetCreateAODs(kTRUE);
	fV0ReaderB->SetUseOwnXYZCalculation(kTRUE);
	fV0ReaderB->SetUseAODConversionPhoton(kTRUE);
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

	// connect input V0Reader
	mgr->AddTask(fV0Reader);
	mgr->ConnectInput (fV0Reader,0,cinput);

	mgr->AddTask(fV0ReaderB);
	mgr->ConnectInput (fV0ReaderB,0,cinput);

	return fV0Reader;
}
