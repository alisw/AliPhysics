AliAnalysisTask *AddTaskCombinedHadronSpectra(Bool_t isMC=kFALSE, Bool_t tpcOnly = kFALSE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_akalweit_CombinedHadron", "No analysis manager found.");
    return 0;
  }
  //============= Set Task Name ===================
  TString taskName=("AliAnalysisCombinedHadronSpectra.cxx+g");
  //===============================================
  //            Load the task
  //gROOT->LoadMacro(taskName.Data());


  
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisCombinedHadronSpectra *task = new AliAnalysisCombinedHadronSpectra("akalweitTaskCombinedHadron");
  task->SelectCollisionCandidates(AliVEvent::kMB);

  
  if (isMC)  task->SetIsMCtrue();
  if (tpcOnly) {
    task->SetUseTPConlyTracks(kTRUE);
    task->Initialize();
  }

  mgr->AddTask(task);


  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  //dumm output container
  AliAnalysisDataContainer *coutput0 =
      mgr->CreateContainer("akalweit_tree",
                           TTree::Class(),
                           AliAnalysisManager::kExchangeContainer,
                           "akalweit_default");

  //define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
      mgr->CreateContainer("akalweit_CombinedHadron", TList::Class(),
                           AliAnalysisManager::kOutputContainer,"akalweit_CombinedHadron.root");

  //connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  0, coutput0);
  mgr->ConnectOutput (task,  1, coutput1);

  return task;
}
