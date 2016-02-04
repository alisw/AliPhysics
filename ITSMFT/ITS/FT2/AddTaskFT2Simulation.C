AliAnalysisTask* AddTaskFT2Simulation(Bool_t ft2OutputTree = 0, Bool_t tuneOnDataOrMC=0,Bool_t simMat=0,TString tpcParameterizationFile="", TString crossSectionFile="",Bool_t usePID=0,Double_t maxStepTGeo=0,Double_t dNdY=0,Bool_t useKalman=0,Bool_t allowDecay=0, Bool_t allowAbsorbtion=0,Bool_t allowConversions=0,TString outputFile="",Int_t runnumber=0,Bool_t standaloneTune=kFALSE, Int_t debugLevel=0,Int_t level=0,TString taskname=""){
	
	// Get the pointer to the existing analysis manager via the static access method.
	//==============================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		::Error("AddTaskFT2Simulation", "No analysis manager to connect to.");
		return NULL;
	}
	//
	// Create event cuts
	//
	Float_t zvWindow = 30. ;
	AliFilteredTreeEventCuts *evtCuts = new AliFilteredTreeEventCuts("AliFilteredTreeEventCuts","Event cuts");
	evtCuts->SetZvRange(-zvWindow,zvWindow);
	evtCuts->SetMeanXYZv(0.0,0.0,0.0);
	evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
	evtCuts->SetTriggerRequired(kFALSE);
	//
	// Create geom. acceptance cuts
	//
	Float_t etaWindow = 1.5 ;
	Float_t ptMin = 0.1 ;
	//
	// Create standard esd track cuts
	//
	AliFilteredTreeAcceptanceCuts *accCuts = new AliFilteredTreeAcceptanceCuts("AliFilteredTreeAcceptanceCuts","Geom. acceptance cuts");
	accCuts->SetEtaRange(-etaWindow,etaWindow);
	accCuts->SetPtRange(ptMin,1.e10);
	//
	// Aanalysis task
	//
	AliAnalysisTaskSEFT2Simulation *TPResidualTask = new AliAnalysisTaskSEFT2Simulation();
	TPResidualTask->SetDebugLevel(debugLevel);
	TPResidualTask->SetEventCuts(evtCuts);
	TPResidualTask->SetAcceptanceCuts(accCuts);
	TPResidualTask->SetMonitorFT2Ouput(ft2OutputTree);
	// FT2 processing flags
	TPResidualTask->SetTuneOnDataOrMC(tuneOnDataOrMC);
	TPResidualTask->SetSimMat(simMat);
	TPResidualTask->SetUsePid(usePID);
	TPResidualTask->SetUseKalman(useKalman);
	TPResidualTask->SetAllowDecay(allowDecay);
	TPResidualTask->SetAllowAbsorbtion(allowAbsorbtion);
	TPResidualTask->SetUseConversionExtension(allowConversions);
	TPResidualTask->SetdNdY(dNdY);
	TPResidualTask->SetMaxStepTGeo(maxStepTGeo);
	TPResidualTask->SetTPCParaFile(tpcParameterizationFile.Data());
	TPResidualTask->SetXSectionFile(crossSectionFile.Data());
	TPResidualTask->SetRunNumber(runnumber);
	TPResidualTask->SetStandaloneTune(standaloneTune);
	TPResidualTask->SetStreamLevel(level);
	mgr->AddTask(TPResidualTask);
	
	//
	// Create containers for input/output
	AliAnalysisDataContainer *coutputmassD08 = mgr->CreateContainer("Simu_hMonitor",TH1F::Class(),AliAnalysisManager::kOutputContainer,outputFile.Data());
	AliAnalysisDataContainer *coutputmassD19 = mgr->CreateContainer("Simu_hCpuTimeEvent",TH1F::Class(),AliAnalysisManager::kOutputContainer,outputFile.Data());
	AliAnalysisDataContainer *coutputmassD20 = mgr->CreateContainer("Simu_hRealTimeEvent",TH1F::Class(),AliAnalysisManager::kOutputContainer,outputFile.Data());
	AliAnalysisDataContainer *coutputmassD21 = mgr->CreateContainer("Simu_hCpuTimeSingleTrack",TH1F::Class(),AliAnalysisManager::kOutputContainer,outputFile.Data());
	AliAnalysisDataContainer *coutputmassD22 = mgr->CreateContainer("Simu_hRealTimeSingleTrack",TH1F::Class(),AliAnalysisManager::kOutputContainer,outputFile.Data());
	AliAnalysisDataContainer *coutputmassD25 = mgr->CreateContainer("esdTree",TTree::Class(),AliAnalysisManager::kOutputContainer,"jstiller_AliESDs.root");
	AliAnalysisDataContainer *coutputmassD26 = mgr->CreateContainer("filteredTree", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());

	mgr->ConnectInput(TPResidualTask,0,mgr->GetCommonInputContainer());
	mgr->ConnectOutput(TPResidualTask,1,coutputmassD08);
	mgr->ConnectOutput(TPResidualTask,2,coutputmassD19);
	mgr->ConnectOutput(TPResidualTask,3,coutputmassD20);
	mgr->ConnectOutput(TPResidualTask,4,coutputmassD21);
	mgr->ConnectOutput(TPResidualTask,5,coutputmassD22);
	mgr->ConnectOutput(TPResidualTask,6,coutputmassD25);
	mgr->ConnectOutput(TPResidualTask,7,coutputmassD26);

	return TPResidualTask;
}