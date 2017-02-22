//______________________________________________________________________________
// creates an AliAnalysisTaskCEP and the output containers
//
AliAnalysisTaskSE* AddTaskCEPAnalysis (
  TString taskname        = TString("CEPAnalysis"),
  Bool_t  isMC            = kFALSE,
  Bool_t  isSaveAll       = kFALSE,
  Int_t   numTracksMax    = 6,
  UInt_t  ETmask          = AliCEPBase::kBitBaseLine,
  UInt_t  ETpattern       = AliCEPBase::kBitBaseLine,
  UInt_t  TTmask          = AliCEPBase::kTTBaseLine,
  UInt_t  TTpattern       = AliCEPBase::kTTBaseLine
  )
 
{

	// get the manager and task
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // check for MonteCarlo
  if (isMC) {
    AliMCEventHandler* MChandler = new AliMCEventHandler;
    MChandler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(MChandler);
  }

  // create the analysis task
  UInt_t taskConfig  = AliCEPBase::kBitConfigurationSet;
  if (isSaveAll)       taskConfig |= AliCEPBase::kBitSaveAllEvents; 
	taskConfig |= AliCEPBase::kBitConfigurationVersion;

	AliAnalysisTaskCEP *task = new AliAnalysisTaskCEP (
    taskname.Data(),taskConfig,
    numTracksMax,ETmask,ETpattern,TTmask,TTpattern );

	// get input and output managers
	AliAnalysisDataContainer *aadci  = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *aadco1 = mgr->CreateContainer
	(
		Form("CEPHist"),
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s", AliAnalysisManager::GetCommonFileName())
	);

	AliAnalysisDataContainer *aadco2 = mgr->CreateContainer
	(
		Form("CEPTree"),
		TTree::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s", AliAnalysisManager::GetCommonFileName())
	);
	
	// add task and connect input and output managers
	mgr->AddTask(task);
	mgr->ConnectInput (task,0,aadci);
	mgr->ConnectOutput(task,1,aadco1);
	mgr->ConnectOutput(task,2,aadco2);

	// return pointer to Task
	return task;

}

//______________________________________________________________________________

