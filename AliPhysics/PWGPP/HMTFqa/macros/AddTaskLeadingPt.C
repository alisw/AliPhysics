AliAnalysisTask* AddTaskLeadingPt(
		Bool_t AnalysisMC = kFALSE,
		Float_t ptleadingCut = 0.5,
		Int_t typerun =0, // 0 for pp and 1 for Pb-Pb or pPb
		TString  type ="ESD",
		UInt_t kTriggerInt = AliVEvent::kMB, //for pPb kINT7, for pp or PbPb kMB
		//UInt_t kTriggerInt = AliVEvent::kINT7, //LHC11c
		Bool_t ispileuprej = kTRUE,
		Bool_t ispileuprejMV = kTRUE,
		Int_t nContributors = 5
		)
{
	// Creates a pid task and adds it to the analysis manager

	// Get the pointer to the existing analysis manager via the static
	//access method
	//=========================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskLeadingPt", "No analysis manager to connect to.");
		return NULL;
	}  

	// Check the analysis type using the event handlers connected to the
	// analysis manager The availability of MC handler can also be
	// checked here.
	// =========================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTaskLeadingPt", "This task requires an input event handler");
		return NULL;
	}  


	AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilter");
	//AliESDtrackCuts* esdTrackCutsGolden = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
	AliESDtrackCuts* esdTrackCutsGolden = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
	trackFilterGolden->AddCuts(esdTrackCutsGolden);

	AliAnalysisFilter* trackFilterTPC = new AliAnalysisFilter("trackFilterTPC");
	AliESDtrackCuts* esdTrackCutsTPC = 
		AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	trackFilterTPC->AddCuts(esdTrackCutsTPC);

	AliAnalysisTaskLeadingPt * taskLeadingPt = new AliAnalysisTaskLeadingPt("taskLeadingPt");
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	taskLeadingPt->SetAnalysisType(type);
	taskLeadingPt->SetPtLeadingCut(ptleadingCut);
	taskLeadingPt->SetDebugLevel(0);
	taskLeadingPt->SetEtaCut(0.8);
	taskLeadingPt->SetVtxCut(10.0);
	taskLeadingPt->SetTrigger(kTriggerInt);
	taskLeadingPt->SetPileUpRej(ispileuprej);
	taskLeadingPt->SetPileUpRejMV(ispileuprejMV);
	taskLeadingPt->SetNcontributors(nContributors);
	//Set Filters
	taskLeadingPt->SetTrackFilterGolden(trackFilterGolden);
	taskLeadingPt->SetTrackFilterTPC(trackFilterTPC);
	mgr->AddTask(taskLeadingPt);

	// Create ONLY the output containers for the data produced by the
	// task.  Get and connect other common input/output containers via
	// the manager as below
	//=======================================================================
	TString outputFileName = AliAnalysisManager::GetCommonFileName();

	AliAnalysisDataContainer *cout_histLeadingPt;
	cout_histLeadingPt=0;
	cout_histLeadingPt = mgr->CreateContainer("outputLeadingPt", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
	mgr->ConnectInput (taskLeadingPt, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(taskLeadingPt, 1, cout_histLeadingPt);

	// Return task pointer at the end
	return taskLeadingPt;



}
