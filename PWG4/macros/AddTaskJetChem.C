AliAnalysisTask *AddTaskJetChem(){



  //cout<<" OB : add JetTasks inlcude path ! "<<endl;
  //gSystem->AddIncludePath("-I$ALICE_ROOT/PWG4/JetTasks");

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_obusch_jets", "No analysis manager found.");
    return 0;
  }

  // physics event selection task - not for AOD analysis

  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(); 
  //AliPhysicsSelection* physSel = physSelTask->GetPhysicsSelection();
  //physSel->AddBackgroundIdentification(new AliBackgroundSelection()); 
  //physSel->SetAnalyzeMC();

  /*
  //============= Set Task Name ===================
  TString taskName=("AliAnalysisTaskJetChem.cxx+");
  //===============================================
  //            Load the task
  gROOT->LoadMacro(taskName.Data());
  if (gProof){
    TString taskSO=gSystem->pwd();
    taskSO+="/";
    taskSO+=taskName(0,taskName.First('.'))+"_cxx.so";
    gProof->Exec(Form("gSystem->Load(\"%s\")",taskSO.Data()),kTRUE);
  }
  */
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskJetChem *task = new AliAnalysisTaskJetChem;
  
  // configure task 


  task->SetUseLOConeJets();
  task->SetUseLOConeMCJets();
  
  //task->SetUsePythiaJets();
  task->SetConeRadius(0.4);
  task->SetTrackPtCutJF(0.150); //
  task->SetFilterBitJF(0x01);   // official PWG4 high pt filter bit 0x10, but not all AliAnalysisTaskESDFilter had configured this ESDTrackCut
  task->SetRequireITSRefitJF(); // 0x01 + ITS refit = 0x10 
  task->SetRejectK0TracksJF();  // uncomment for K0 analysis running jet finder in task - modifies jet spectrum

  task->SetJetPtCut(2.0);
  //task->SetJetPtCut(0.150); // lower pt cut: for plot of diffractive contribution to jet spectrum (goes up to 2 GeV ...) 
  task->SetJetEtaCut(0.5);

  task->SetFilterBit((UInt_t) 0X01); // std AOD track cuts
  task->SetTrackPtCut(0.150); 
  task->SetTrackEtaCut(0.9);

  task->SetUseOnFlyV0s(); 
  task->SetCutnSigdEdx(2); 
 
  //task->ReadDeltaAOD(); // uncomment for DeltaAODs
  //task->SelectDeltaAODBranch("bla"); 
  //task->SelectAODBranch("jetsAOD_FASTKT04");
  //task->SelectAODBranch("jetsAOD_UA107");
  //task->SelectAODBranch("jets");

  task->SelectCollisionCandidates(); // either here or in userExec of task - but not for AODs ...

  //  AliLog::SetGlobalLogLevel(AliLog::kInfo);    // kInfo // kDebug // kFatal
  //  task->SetDebugLevel(10);
  //  mgr->SetDebugLevel(10);

  mgr->AddTask(task);


  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("PWG4_JetChem", 
							   TList::Class(), AliAnalysisManager::kOutputContainer,
							   Form("%s:PWG4_JetChem",AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput(task,0,cinput );
  mgr->ConnectOutput(task,1,coutput);
 
  return task;
}
