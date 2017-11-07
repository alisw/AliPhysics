//Authors: Sergio Iga ,sergio.iga@correo.nucleares.unam.mx


AliAnalysisTask* AddTaskPPvsMultINEL0(
		Bool_t AnalysisMC = kTRUE,
		TString  type ="ESD",
		UInt_t kTriggerInt = AliVEvent::kINT7, //kINT7 is the good trigger
		Bool_t ispileuprej = kTRUE)
{
	// Get the pointer to the existing analysis manager via the static
	//access method
	//=========================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskPPvsMultINEL0", "No analysis manager to connect to.");
		return NULL;
	}  

	// Check the analysis type using the event handlers connected to the
	// analysis manager The availability of MC handler can also be
	// checked here.
	// =========================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTaskPPvsMultINEL0", "This task requires an input event handler");
		return NULL;
	}  

        gROOT->LoadMacro("$(ALICE_PHYSICS)/PWGJE/macros/CreateTrackCutsPWGJE.C");
	
	/*
	AliAnalysisFilter* trackFilterTPC = new AliAnalysisFilter("trackFilterTPC");
	AliESDtrackCuts* esdTrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	trackFilterTPC->AddCuts(esdTrackCutsTPC);

	AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilterGolden");
	AliESDtrackCuts* esdTrackCutsGolden = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);	
	trackFilterGolden->AddCuts(esdTrackCutsGolden);
	*/
	
	

	AliAnalysisTaskPPvsMultINEL0 * taskPPvsMultINEL0 = new AliAnalysisTaskPPvsMultINEL0("taskPPvsMultINEL0");


	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	taskPPvsMultINEL0->SetAnalysisMC(AnalysisMC);
	taskPPvsMultINEL0->SetAnalysisType(type);
	taskPPvsMultINEL0->SetDebugLevel(0);
	taskPPvsMultINEL0->SetEtaCut(0.8);
	taskPPvsMultINEL0->SetVtxCut(10.0);
	taskPPvsMultINEL0->SetTrigger(kTriggerInt);
	taskPPvsMultINEL0->SetPileUpRej(ispileuprej);		

	mgr->AddTask(taskPPvsMultINEL0);

	// Create ONLY the output containers for the data produced by the
	// task.  Get and connect other common input/output containers via
	// the manager as below
	//=======================================================================
	TString outputFileName = AliAnalysisManager::GetCommonFileName();

	AliAnalysisDataContainer *cout_histPPvsMultINEL0;
	cout_histPPvsMultINEL0=0;
	cout_histPPvsMultINEL0 = mgr->CreateContainer("outputPPvsMultINEL0", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
	mgr->ConnectInput (taskPPvsMultINEL0, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(taskPPvsMultINEL0, 1, cout_histPPvsMultINEL0);

	
	// Return task pointer at the end
	return taskPPvsMultINEL0;



}
