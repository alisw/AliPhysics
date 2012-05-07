

AliAnalysisTask *AddTaskCombinedHadronSpectra(Int_t identifier = 0, Bool_t isMC = kFALSE, Bool_t isTPConly = kFALSE, Bool_t writeOwnFile = kFALSE, Bool_t saveMotherPDG = kFALSE, Bool_t smallTHnSparse = kFALSE, Double_t nSigmaTPCLow= -3., Double_t nSigmaTPCHigh = 3., Double_t rapidityLow = -0.2, Double_t rapidityHigh = 0.2, Bool_t setTrackCuts = kFALSE, AliESDtrackCuts *ESDtrackCuts = 0){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_janielsk_CombinedHadron", "No analysis manager found.");
    return 0;
  }
  //============= Set Task Name ===================
  TString taskName=("AliAnalysisCombinedHadronSpectra.cxx+g");
  //===============================================
  //            Load the task
  gROOT->LoadMacro(taskName.Data());


  
  //========= Add task to the ANALYSIS manager =====

  //normal tracks
  AliAnalysisCombinedHadronSpectra *task = new AliAnalysisCombinedHadronSpectra("janielskTaskCombinedHadron");
  task->SelectCollisionCandidates(AliVEvent::kMB);

  //switches
  if (isMC) task->SetIsMCtrue(isMC);
  if (isTPConly)task->SetUseTPConlyTracks(isTPConly);
  if (saveMotherPDG) task->SetSaveMotherPDG(saveMotherPDG);
  if (smallTHnSparse){
    task->SetSmallTHnSparse(kTRUE);
    task->SetTPCnSigmaCuts(-3.,3.);
    task->SetRapidityCuts(-0.2,0.2);
  }

  //initialize task
  task->Initialize();

  //esd cuts need to be set after initialize or cuts will be replaced by standard cuts in initialize
  if (setTrackCuts) task->SetESDtrackCuts(ESDtrackCuts);

  //add task to manager
  mgr->AddTask(task);


  

  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

/*
  //dumm output container
  AliAnalysisDataContainer *coutput0 =
      mgr->CreateContainer(Form("akalweit_tree%i",identifier),
                           TTree::Class(),
                           AliAnalysisManager::kExchangeContainer,
                           Form("akalweit_default%i",identifier));


  //define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
      mgr->CreateContainer(Form("akalweit_CombinedHadron%i",identifier), TList::Class(),
                           AliAnalysisManager::kOutputContainer,Form("akalweit_CombinedHadron%i.root",identifier));
*/
  if (!writeOwnFile) {
  	AliAnalysisDataContainer *coutput1 =  mgr->CreateContainer(Form("janielsk_CombinedHadron%i",identifier), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:janielsk_CombinedHadron", AliAnalysisManager::GetCommonFileName())); 
	}
  else {
	AliAnalysisDataContainer *coutput1 =  mgr->CreateContainer(Form("janielsk_CombinedHadron%i",identifier), TList::Class(), AliAnalysisManager::kOutputContainer, Form("janielsk_CombinedHadron.root"));
	}



  //connect containers

  //
  mgr->ConnectInput  (task,  0, cinput );
  //mgr->ConnectOutput (task,  0, coutput0);
  mgr->ConnectOutput (task,  1, coutput1);

  return task;
}
