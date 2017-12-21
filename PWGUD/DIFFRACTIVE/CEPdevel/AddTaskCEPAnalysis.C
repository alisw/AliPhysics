//______________________________________________________________________________
// creates an AliAnalysisTaskCEP and the output containers
//
AliAnalysisTaskSE* AddTaskCEPAnalysis (
  TString taskname        = TString("CEPAnalysis"),
  UInt_t  taskconfig      = AliCEPBase::kBitConfigurationSet,
  Int_t rnummin,
  Int_t rnummax,
  Int_t   numTracksMax    = 6,
  Double_t fracDG         = 1.0,
  Double_t fracNDG        = 0.0,
  UInt_t  ETmaskDG        = AliCEPBase::kETBaseLine,
  UInt_t  ETpatternDG     = AliCEPBase::kETBaseLine,
  UInt_t  ETmaskNDG       = AliCEPBase::kETBaseLine,
  UInt_t  ETpatternNDG    = AliCEPBase::kETBaseLine,
  UInt_t  TTmask          = AliCEPBase::kTTBaseLine,
  UInt_t  TTpattern       = AliCEPBase::kTTBaseLine
  )
{

	// get the manager and task
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // check for MonteCarlo
  Bool_t isMC = (taskconfig&AliCEPBase::kBitisMC)==AliCEPBase::kBitisMC;
  if (isMC) {
    AliMCEventHandler* MChandler = new AliMCEventHandler;
    MChandler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(MChandler);
  }

  // create the analysis task
	AliAnalysisTaskCEP *task = new AliAnalysisTaskCEP (
    taskname.Data(),taskconfig,
    rnummin,rnummax,
    numTracksMax,
    fracDG,fracNDG,
    ETmaskDG,ETpatternDG,ETmaskNDG,ETpatternNDG,TTmask,TTpattern
  );

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

