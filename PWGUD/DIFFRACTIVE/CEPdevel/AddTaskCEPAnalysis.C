//______________________________________________________________________________

AliAnalysisTaskCEPAnalysis* AddAnalysisTaskCEP(
  Bool_t  isMC         = kFALSE,
  Bool_t  isSaveAll    = kFALSE,
  Int_t   numTracksMax = 100000)
 
  )
{

	// get the manager and task
	AliAnalysisManager *aam = AliAnalysisManager::GetAnalysisManager();

  // check for MonteCarlo
  if (isMC) {
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // create the analysis task
  UInt_t taskConfig  = AliCEPBase::kBitConfigurationSet;
  if (isSaveAll) taskConfig |= AliCEPBase::kBitSaveAllEvents; 
	taskConfig |= AliCEPBase::kBitConfigurationVersion;

  TString name = TString("CEPAnalysis");
	AliAnalysisTaskCEP *task = new AliAnalysisTaskCEP (
    name.Data(),taskConfig, numTracksMax);

	// get input and output managers
	AliAnalysisDataContainer *aadci = aam->GetCommonInputContainer();
	AliAnalysisDataContainer *aadco1 = aam->CreateContainer
	(
		Form("CEPHist"),
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:CEP", AliAnalysisManager::GetCommonFileName())
	);

	AliAnalysisDataContainer *aadco2 = aam->CreateContainer
	(
		Form("CEPTree"),
		TTree::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:CEP", AliAnalysisManager::GetCommonFileName())
	);
	
	// add task and connect input and output managers
	aam->AddTask(task);
	aam->ConnectInput (task,0,aadci);
	aam->ConnectOutput(task,1,aadco1);
	aam->ConnectOutput(task,2,aadco2);

	// return pointer to Task
	return task;

}

//______________________________________________________________________________

