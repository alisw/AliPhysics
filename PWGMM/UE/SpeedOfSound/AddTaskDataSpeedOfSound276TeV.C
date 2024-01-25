AliAnalysisTask* AddTaskDataSpeedOfSound276TeV(
    Bool_t AnalysisMC = kFALSE,
    Int_t typerun = 1,  // 0 for pp and 1 for Pb-Pb or pPb
    TString type = "ESD",
    UInt_t kTriggerInt = AliVEvent::kMB,  // for pPb kINT7, for pp or PbPb kMB
    Float_t minCent = 0., Float_t maxCent = 80.,
    const char* centralityEstimator = "V0M",  // for pPb V0A for PbPb V0M
    Bool_t ispileuprej = kTRUE, const char* suffix = "") {
  // Creates a pid task and adds it to the analysis manager

  // Get the pointer to the existing analysis manager via the static
  // access method
  //=========================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHighPtDeDx", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskHighPtDeDx", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts* esdTrackCutsGolden =
      AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 1);
  trackFilterGolden->AddCuts(esdTrackCutsGolden);

  AliAnalysisFilter* trackFilterTPC = new AliAnalysisFilter("trackFilterTPC");
  AliESDtrackCuts* esdTrackCutsTPC =
      AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  trackFilterTPC->AddCuts(esdTrackCutsTPC);

  AliAnalysisTaskDataSpeedOfSound276TeV* mytask =
      new AliAnalysisTaskDataSpeedOfSound276TeV("SpeedOfSound");
  mytask->SetAnalysisType(type);
  mytask->SetAnalysisMC(AnalysisMC);
  if (typerun == 1) {
    mytask->SetAnalysisPbPb(kTRUE);
    mytask->SetMinCent(minCent);
    mytask->SetMaxCent(maxCent);
    mytask->SetCentralityEstimator(centralityEstimator);
  } else {
    mytask->SetAnalysisPbPb(kFALSE);
  }
  mytask->SetDebugLevel(0);
  mytask->SetEtaCut(0.8);
  mytask->SetPtMin(0.15);
  mytask->SetEtaGappT(0.4);
  mytask->SetEtaGapNch(0.7, 1.4);
  mytask->SetVtxCut(10.0);
  mytask->SetTrigger(kTriggerInt);
  mytask->SetPileUpRej(ispileuprej);
  mytask->SetStoreMcIn(AnalysisMC);  // def: kFALSE
  mgr->AddTask(mytask);

  const char* taskname = "SpeedOfSound";
  mgr->ConnectInput(mytask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(
      mytask, 1,
      mgr->CreateContainer(
          Form("cList%s_%s", taskname, suffix), TList::Class(),
          AliAnalysisManager::kOutputContainer,
          Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname)));

  return mytask;
}
