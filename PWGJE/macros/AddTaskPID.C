AliAnalysisTask *AddTaskPID(TString nameSuffix, Bool_t writeOutputToSeparateFiles = kTRUE,
                            Bool_t useConvolutedGauss = kTRUE, TString centralityEstimator = "V0A",
                            Bool_t considerJets = kTRUE, Bool_t overrideStoreCentralityPercentile = kFALSE,
                            Bool_t overrideStoreCentralityPercentileValue = kFALSE,
                            TString listOfFiles = "")
{
  // Macro to set up and add PID task with default settings.
  //
  // Typical parameters to run on 11a1* (MC_pp@7TeV):
  // "PWGJE_taskPID_Jets", kTRUE, kTRUE, "V0A", kTRUE
  // and as a second task
  // "PWGJE_taskPID_Jets_Inclusive", kTRUE, kTRUE, "V0A", kTRUE
  
  TString taskName = Form("PWGJE_taskPID%s", nameSuffix.Data());
  
  
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPID", "No analysis manager found.");
    return 0x0;
  }
  
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskPID *task = new AliAnalysisTaskPID(taskName.Data());
  task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kINT7);
  
  
  
  printf("\nSetting up task %s:\n", taskName.Data());
  
  if (!considerJets) {
    //
    // Add track filters
    //
    
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    AliESDtrackCuts* esdTrackCutsL = 0x0;
    
    if (listOfFiles.Contains("LHC11") || listOfFiles.Contains("LHC12") || listOfFiles.Contains("LHC13")) {
      esdTrackCutsL = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
      printf("Using standard ITS-TPC track cuts 2011\n");
      trackFilter->SetTitle("Standard ITS-TPC track cuts 2011");
    }
    else if (listOfFiles.Contains("LHC10")) {
      esdTrackCutsL = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
      printf("Using standard ITS-TPC track cuts 2010\n");
      trackFilter->SetTitle("Standard ITS-TPC track cuts 2010");
    }
    else  {
      esdTrackCutsL = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
      printf("WARNING: Cuts not configured for this period!!! Using standard ITS-TPC track cuts 2011\n");
      trackFilter->SetTitle("Standard ITS-TPC track cuts 2011");
    }
    /*
    esdTrackCutsL->SetMinNCrossedRowsTPC(120);
    esdTrackCutsL->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCutsL->SetMaxChi2PerClusterITS(36);
    esdTrackCutsL->SetMaxFractionSharedTPCClusters(0.4);
    esdTrackCutsL->SetMaxChi2TPCConstrainedGlobal(36);
    */
    
    trackFilter->AddCuts(esdTrackCutsL);
    task->SetTrackFilter(trackFilter);
  }
  
  // Test whether we have pPb or Pbp
  if (listOfFiles.Contains("pPb") || listOfFiles.Contains("Pbp")) {
    task->SetIsPbpOrpPb(kTRUE);
    printf("pPb/Pbp detected -> Adapting vertex cuts!\n");
    task->SetCentralityEstimator("V0A");
  }
  else  {
    task->SetIsPbpOrpPb(kFALSE);
    printf("Collision type different from pPb/Pbp detected -> Using standard vertex cuts!\n");
  }
  
  // Do not store centrality percentile for pp (will be set to -1 for all events) - or use special centrality estimator and store
  if (listOfFiles.Contains("pp")) {
    //task->SetStoreCentralityPercentile(kFALSE);
    
    task->SetStoreCentralityPercentile(kTRUE);
    task->SetCentralityEstimator("ITSTPCtracklets");
  }
  else
    task->SetStoreCentralityPercentile(kTRUE);
  
  if (overrideStoreCentralityPercentile)
    task->SetStoreCentralityPercentile(overrideStoreCentralityPercentileValue);
  
  task->SetEtaAbsCutRange(0.0, 0.9);
  task->SetUsePhiCut(kFALSE);
  
  // Do not set the covolution parameters if they are not used (saves some cpu time for the initialisation)
  if (useConvolutedGauss) {
    if ((listOfFiles.Contains("pPb") || listOfFiles.Contains("Pbp")) && listOfFiles.Contains("LHC13")) {
      printf("\n13* pPb @ 5.023 ATeV data detected -> Setting corresponding convolution parameters...\n");
      task->SetConvolutedGaussLambdaParameter(3.0);
    }
    
    if (listOfFiles.Contains("LHC11a_without_SDD")) {
      printf("\n11a (pp 2.76 TeV) detected -> Setting corresponding convolution parameters...\n");
      task->SetConvolutedGaussLambdaParameter(3.0);
    }
    else {
      printf("\nUsing default convolution parameters...\n");
      task->SetConvolutedGaussLambdaParameter(2.0);
    }
  }
  
  task->SetCentralityEstimator(centralityEstimator.Data());
  
  task->SetUseMCidForGeneration(kTRUE);
  task->SetUseITS(kTRUE);
  task->SetUseTOF(kTRUE);
  task->SetUsePriors(kTRUE);
  task->SetUseTPCDefaultPriors(kTRUE);
  task->SetUseConvolutedGaus(useConvolutedGauss);
  task->SetTakeIntoAccountMuons(kTRUE);
  
  task->SetInputFromOtherTask(considerJets);
  task->SetStoreAdditionalJetInformation(considerJets);
  
  task->PrintSettings();
  
  mgr->AddTask(task);


  //================================================
  //              data containers
  //================================================

  //define output containers
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer(Form("%s", taskName.Data()),
                         TObjArray::Class(),
                         AliAnalysisManager::kOutputContainer,
                         writeOutputToSeparateFiles
                          ? Form("%s.root", taskName.Data())
                          : Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskName.Data()));
  
  if (task->GetDoEfficiency()) {
    AliAnalysisDataContainer *coutput2 = 
        mgr->CreateContainer(Form("%s_efficiency", taskName.Data()),
                            AliCFContainer::Class(),
                            AliAnalysisManager::kOutputContainer,
                            writeOutputToSeparateFiles
                              ? Form("%s_efficiency.root", taskName.Data())
                              : Form("%s:%s_efficiency", AliAnalysisManager::GetCommonFileName(), taskName.Data()));
    mgr->ConnectOutput (task,  2, coutput2);
  }
  
  if (task->GetDoDeDxCheck() || task->GetDoPtResolution()) {
    AliAnalysisDataContainer *coutput3 = 
        mgr->CreateContainer(Form("%s_PtResolution", taskName.Data()),
                            TObjArray::Class(),
                            AliAnalysisManager::kOutputContainer,
                            writeOutputToSeparateFiles
                              ? Form("%s_PtResolution.root", taskName.Data())
                              : Form("%s:%s_PtResolution", AliAnalysisManager::GetCommonFileName(), taskName.Data()));
    mgr->ConnectOutput (task,  3, coutput3);
  }
  //connect containers
  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  0, mgr->GetCommonOutputContainer()); // comment to run local
  mgr->ConnectOutput (task,  1, coutput1);

  return task;
}
