void AddTask_mmarquar_mptPbPb_lowmult()
{
	/*
	   CheckLoadLibrary("libPWG0base");
	   CheckLoadLibrary("libPWG0dep");
	   CheckLoadLibrary("libPWG0selectors");
	 */

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	if (!mgr) {
		Error("AddTask_mmarquar_mptPbPb_lowmult", "No analysis manager found.");
		return 0;
	}

	// Switch off all AliInfo (too much output!!!)
	AliLog::SetGlobalLogLevel(AliLog::kError);
	mgr->SetDebugLevel(0);

	//
	// Create physics trigger selection class
	//


	// AliPhysicsSelection *physTrigSel =  new AliPhysicsSelection();

	//
	// Create event cuts
	//

	Float_t zvWindow = 30. ;

	AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
	evtCuts->SetZvRange(-zvWindow,zvWindow);
	evtCuts->SetMeanXYZv(0.0,0.0,0.0);
	evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
	evtCuts->SetTriggerRequired(kTRUE);
	evtCuts->SetEventSelectedRequired(kFALSE);

	//
	// Create geom. acceptance cuts
	//
	Float_t etaWindow = 0.3;
	Float_t ptMin = 0.15 ;

	AlidNdPtAcceptanceCuts *accCuts = new AlidNdPtAcceptanceCuts("AlidNdPtAcceptanceCuts","Geom. acceptance cuts");
	accCuts->SetEtaRange(-etaWindow,etaWindow);
	accCuts->SetPtRange(ptMin,1.e10);

	//
	// Create standard esd track cuts
	//
	Int_t cutMode = 200;

	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C");
	//gROOT->LoadMacro("./CreatedNdPtTrackCuts.C");
	AliESDtrackCuts* esdTrackCuts = CreatedNdPtTrackCuts(cutMode);
	if (!esdTrackCuts) {
		printf("ERROR: esdTrackCuts could not be created\n");
		return;
	} else {
		//esdTrackCuts->SetHistogramsOn(kTRUE);
		esdTrackCuts->SetHistogramsOn(kFALSE);
	}
	esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);
	esdTrackCuts->SetMaxChi2PerClusterITS(36.);

	AliESDtrackCuts* esdmultTrackCuts = CreatedNdPtTrackCuts(201);

	esdmultTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);
	esdmultTrackCuts->SetMaxChi2PerClusterITS(36.);



	Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

	//
	// Create task
	//
	AlidNdPtTask *task = new AlidNdPtTask("AlidNdPtTask_TPCITS");
	task->SetUseMCInfo(hasMC);

	// trigger
	// task->SelectCollisionCandidates(AliVEvent::kMB); 

	//
	// set analysis options from the Helper here !!!
	//

	//AliTriggerAnalysis::Trigger trigger = AliTriggerAnalysis::kMB1;
	AlidNdPtHelper::AnalysisMode analysisMode = AlidNdPtHelper::kTPCITS;
	AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kAllPart ;

	//
	// Create analysis object
	//

	AlimPtAnalysis *fmPtAnalysis = new AlimPtAnalysis("mPtAnalysis_TPCITS","dN/dPt Analysis with TPC-ITS tracking");
	fmPtAnalysis->SetmultTrackCuts(esdmultTrackCuts);
	fmPtAnalysis->SetEventCuts(evtCuts);
	fmPtAnalysis->SetAcceptanceCuts(accCuts);
	fmPtAnalysis->SetTrackCuts(esdTrackCuts);
	//fmPtAnalysis->SetBackgroundCuts(backCuts);
	fmPtAnalysis->SetAnalysisMode(analysisMode); 
	fmPtAnalysis->SetParticleMode(particleMode);
	//fmPtAnalysis->SetTrigger(trigger);
	fmPtAnalysis->SetTriggerMask(AliVEvent::kMB);
	//fmPtAnalysis->SetTriggerMask(AliVEvent::kEMC1);
	fmPtAnalysis->SetCentSelection(kFALSE);
	//fmPtAnalysis->SetCentLimit(8);
	if(hasMC) 
	{
		//physTrigSel->SetAnalyzeMC();
		//fmPtAnalysis->SetPhysicsTriggerSelection(physTrigSel); 

		fmPtAnalysis->SetUseMCInfo(kTRUE);
		fmPtAnalysis->SetHistogramsOn(kTRUE);
		//fmPtAnalysis->SetHistogramsOn(kFALSE);
	}
	else { // online trigger
		//     physTrigSel->SetUseBXNumbers();
		//     physTrigSel->SetComputeBG();
		//     fmPtAnalysis->SetPhysicsTriggerSelection(physTrigSel); 
	}

	// change binning
	Double_t binsMult[402];
	Int_t multNbins = 401;
	for (int i=0; i<=multNbins; i++) { binsMult[i] = -0.5 + i; }

	const Int_t ptNbins = 73;
	Double_t bins[74] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 };
	Double_t* binsPt = new Double_t [74];
	for (int i=0; i<74; i++) {binsPt[i] = bins[i];}
	fmPtAnalysis->SetBinsPt(ptNbins, binsPt);
	fmPtAnalysis->SetBinsPtCorr(ptNbins, binsPt);  
	fmPtAnalysis->SetBinsMult(multNbins, binsMult);


	// Add analysis object
	task->AddAnalysisObject( fmPtAnalysis );

	// Add task
	mgr->AddTask(task);

	// Create containers for input
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	mgr->ConnectInput(task, 0, cinput);

	AliAnalysisDataContainer *coutput = mgr->CreateContainer("mmarquar_mptPbPb_lowmult", TList::Class(), AliAnalysisManager::kOutputContainer, "mmarquar_mptPbPb_lowmult.root");
	mgr->ConnectOutput(task, 1, coutput);

}

