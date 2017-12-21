// ### Settings that make sense when using the Alien plugin
//==============================================================================
Int_t       iCollision         = 0;       // 0=pp, 1=Pb-Pb
//==============================================================================
Bool_t      usePhysicsSelection = kFALSE; // use physics selection
Bool_t      useCentrality       = kFALSE; // centrality
Bool_t      useDBG              = kFALSE; // activate debugging
Bool_t      useMC               = kTRUE;  // use MC info
Bool_t      useKFILTER          = kTRUE;  // use Kinematics filter
Bool_t      useTR               = kTRUE;  // use track references
Bool_t      useSysInfo          = kFALSE; // use sys info

// ### Analysis modules to be included. Some may not be yet fully implemented.
//==============================================================================
Int_t       iAODhandler        = 1;      // Analysis produces an AOD or dAOD's
Int_t       iMUONCDBConnect    = 0;      // Task to load MUON OCDB objects                       (> 1 = use parfile)
Int_t       iESDMCLabelAddition= 1;      // Recompute MC labels for MUON                         (> 1 = use parfile)
Int_t       iESDfilter         = 1;      // ESD to AOD filter (barrel + muon tracks)
Int_t       iMUONcopyAOD       = 1;      // Task that copies only muon events in a separate AOD (PWG3)
Int_t       iMUONRefit         = 0;      // Refit ESD muon tracks before producing AODs          (> 1 = use parfile)
Int_t       iMUONRefitVtx      = 0;      // Refit ESD muon tracks at vtx before producing AODs   (> 1 = use parfile)
Int_t       iMUONQA            = 1;      // run muon QA task on ESDs                             (> 1 = use parfile)
Int_t       iMUONPerformance   = 1;      // Task to study the muon performances in MC simulation (> 1 = use parfile)
Int_t       iMUONEfficiency    = 1;      // Task to measure the muon efficiency                  (> 1 = use parfile)

// ### Configuration macros used for each module and OCDB settings
//==============================================================================
TString defaultStorage = VAR_OCDB_PATH;
TString alignStorage = VAR_REC_ALIGNDATA;
Int_t alignVersion = -1;
Int_t alignSubVersion = -1;
TString recoParamStorage = "";

// Temporaries.
void AODmerge();
void AddAnalysisTasks(Int_t);
Bool_t LoadAnalysisLibraries();
TChain *CreateChain();

//______________________________________________________________________________
void AODtrainsim(Int_t merge=0)
{
  // Main analysis train macro.
  // merge = 0: production
  // merge = 1: intermediate merging
  // merge = 2: final merging + terminate
  // merge = 3: terminate only

  if (merge) {
    TGrid::Connect("alien://");
    if (!gGrid || !gGrid->IsConnected()) {
      ::Error("AODtrainsim.C", "No grid connection");
      return;
    }
  }
  // Set temporary merging directory to current one
  gSystem->Setenv("TMPDIR", gSystem->pwd());
  // Set temporary compilation directory to current one
  gSystem->SetBuildDir(gSystem->pwd(), kTRUE);
  printf("==================================================================\n");
  printf("===========    RUNNING FILTERING TRAIN   ==========\n");
  printf("==================================================================\n");
  printf("=  Configuring analysis train for:                               =\n");
  if (usePhysicsSelection)   printf("=  Physics selection                                                =\n");
  if (iESDfilter)   printf("=  ESD filter                                                    =\n");
  if (iMUONcopyAOD) printf("=  MUON copy AOD                                                 =\n");

  // Load common libraries and set include path
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  // Make the analysis manager and connect event handlers
  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Production train");
  if (useSysInfo) mgr->SetNSysInfo(100);

  // Load ParFiles
  if (!LoadAnalysisLibraries()) {
    ::Error("AODtrainsim.C", "Could not load analysis libraries");
    return;
  }
  
  // Create input handler (input container created automatically)
  // ESD input handler
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdHandler);
  // Monte Carlo handler
  if (useMC) {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
    mcHandler->SetReadTR(useTR);
  }
  // AOD output container, created automatically when setting an AOD handler
  if (iAODhandler) {
    // AOD output handler
    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName("AliAOD.root");
    mgr->SetOutputEventHandler(aodHandler);
  }
  // Debugging if needed
  if (useDBG) mgr->SetDebugLevel(3);

  AddAnalysisTasks(merge);
  if (merge) {
    if (merge < 3) AODmerge();
    if (merge > 1) {
      mgr->InitAnalysis();
      mgr->SetGridHandler(new AliAnalysisAlien());
      mgr->StartAnalysis("grid terminate",0);
    }
    return;
  }
  // Run the analysis
  //
  TChain *chain = CreateChain();
  if (!chain) return;

  TStopwatch timer;
  timer.Start();
  mgr->SetSkipTerminate(kTRUE);
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  timer.Print();
}

//______________________________________________________________________________
void AddAnalysisTasks(Int_t merge){
  // Add all analysis task wagons to the train
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  //
  // Tender and supplies. Needs to be called for every event.
  //

  if ( VAR_OCDB_SNAPSHOT )
  {
      AliCDBManager::Instance()->SetDefaultStorage("local:///donotexist"); 
      // explicitely wrong default storage
      // so we take _EVERYTHING_ from the snapshot.
      AliCDBManager::Instance()->SetRun(0);
      AliCDBManager::Instance()->SetSnapshotMode("OCDB_rec.root");
  }

  if ( VAR_MAKE_COMPACT_ESD )
  {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/muondep/AddTaskCompactTreeMaker.C");
      AddTaskCompactTreeMaker();
  }
  
  if (iMUONCDBConnect) {
    if (iMUONCDBConnect > 1) gROOT->LoadMacro("AddTaskMuonCDBConnect.C");
    else gROOT->LoadMacro("$ALICE_PHYSICS/PWG/muondep/AddTaskMuonCDBConnect.C");
    AliAnalysisTaskMuonCDBConnect *cdbConnect = AddTaskMuonCDBConnect();
    if (!defaultStorage.IsNull()) cdbConnect->SetDefaultStorage(defaultStorage.Data());
    if (!alignStorage.IsNull() || alignVersion >= 0 || alignSubVersion >= 0)
      cdbConnect->SetAlignStorage(alignStorage.Data(), alignVersion, alignSubVersion);
    if (!recoParamStorage.IsNull()) cdbConnect->SetRecoParamStorage(recoParamStorage.Data());
    cdbConnect->LoadMagField();
    cdbConnect->LoadGeometry();
    cdbConnect->LoadMapping();
  }
  
  UInt_t offlineTriggerMask = 0;
  if (usePhysicsSelection) {
    // Physics selection task
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    mgr->RegisterExtraFile("event_stat.root");
    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(kTRUE);
    offlineTriggerMask = AliVEvent::kAny;
  }
  
  // Centrality (only Pb-Pb)
  if (iCollision && useCentrality) {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask *mult = AddTaskMultSelection(kFALSE);
  }
  
  // track selection
  AliMuonTrackCuts *trackCuts = 0x0;
  if (iMUONEfficiency) {
    trackCuts = new AliMuonTrackCuts("stdCuts", "stdCuts");
    trackCuts->SetAllowDefaultParams();
    trackCuts->SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
    trackCuts->SetIsMC(kTRUE);
  }
  
  if (iMUONRefit) {
    if (iMUONRefit > 1) gROOT->LoadMacro("AddTaskMuonRefit.C");
    else gROOT->LoadMacro("$ALICE_PHYSICS/PWG/muondep/AddTaskMuonRefit.C");
    AliAnalysisTaskMuonRefit* refit = AddTaskMuonRefit(-1., -1., kTRUE, -1., -1.);
    if (!defaultStorage.IsNull()) refit->SetDefaultStorage(defaultStorage.Data());
    if (!alignStorage.IsNull()) refit->SetAlignStorage(alignStorage.Data());
    refit->RemoveMonoCathodClusters(kTRUE, kFALSE);
  }
  
  if (iMUONRefitVtx) {
    if (iMUONRefitVtx > 1) gROOT->LoadMacro("AddTaskMuonRefitVtx.C");
    else gROOT->LoadMacro("$ALICE_PHYSICS/PWG/muondep/AddTaskMuonRefitVtx.C");
    AliAnalysisTaskMuonRefitVtx* refitVtx = AddTaskMuonRefitVtx(kFALSE, kFALSE, kTRUE);
    if (!defaultStorage.IsNull()) refitVtx->SetDefaultStorage(defaultStorage.Data());
  }
  
  if(iESDMCLabelAddition) {
    if(iESDMCLabelAddition > 1) gROOT->LoadMacro("AddTaskESDMCLabelAddition.C");
    else gROOT->LoadMacro("$ALICE_PHYSICS/PWG/muondep/AddTaskESDMCLabelAddition.C");
    AliAnalysisTaskESDMCLabelAddition *esdmclabel = AddTaskESDMCLabelAddition();
    if (!defaultStorage.IsNull()) esdmclabel->SetDefaultStorage(defaultStorage.Data());
    if (!alignStorage.IsNull()) esdmclabel->SetAlignStorage(alignStorage.Data());
    if (!recoParamStorage.IsNull()) esdmclabel->SetRecoParamStorage(recoParamStorage.Data());
    esdmclabel->DecayAsFake(kTRUE);
  }
  
  if (iMUONQA) {
    if (iMUONQA > 1) gROOT->LoadMacro("AddTaskMuonQA.C");
    else gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskMuonQA.C");
    AliAnalysisTaskMuonQA* muonQA = AddTaskMuonQA(kFALSE);
    if (usePhysicsSelection) muonQA->SelectCollisionCandidates(offlineTriggerMask);
    muonQA->SetTrackCuts(trackCuts);
  }
  
  if (useMC && useTR && iMUONPerformance) {
    if (iMUONPerformance > 1) gROOT->LoadMacro("AddTaskMuonPerformance.C");
    else gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/MUON/dep/AddTaskMuonPerformance.C");
    AliAnalysisTaskMuonPerformance* muonPerformance = AddTaskMuonPerformance();
    if (usePhysicsSelection) muonPerformance->SelectCollisionCandidates(offlineTriggerMask);
    if (!defaultStorage.IsNull()) muonPerformance->SetDefaultStorage(defaultStorage.Data());
    if (!alignStorage.IsNull()) muonPerformance->SetAlignStorage(alignStorage.Data());
    if (!recoParamStorage.IsNull()) muonPerformance->SetRecoParamStorage(recoParamStorage.Data());
    muonPerformance->UseMCKinematics(kTRUE);
    muonPerformance->SetMCTrigLevelFromMatchTrk(kTRUE);
  }
  
  if (iMUONEfficiency) {
    if (iMUONEfficiency > 1) gROOT->LoadMacro("AddTaskMUONTrackingEfficiency.C");
    else gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/MUON/dep/AddTaskMUONTrackingEfficiency.C");
    AliAnalysisTaskMuonTrackingEff* muonEfficiency = AddTaskMUONTrackingEfficiency(kFALSE,kTRUE,"");
    if (usePhysicsSelection) muonEfficiency->SelectCollisionCandidates(offlineTriggerMask);
    if (!defaultStorage.IsNull()) muonEfficiency->SetDefaultStorage(defaultStorage.Data());
    if (!alignStorage.IsNull()) muonEfficiency->SetAlignStorage(alignStorage.Data());
    if (!recoParamStorage.IsNull()) muonEfficiency->SetRecoParamStorage(recoParamStorage.Data());
    muonEfficiency->SetMuonTrackCuts(*trackCuts);
    muonEfficiency->SetMuonPtCut(VAR_EFFTASK_PTMIN);
    muonEfficiency->UseMCLabel(kTRUE);
    muonEfficiency->EnableDisplay(kFALSE);
  }
  
  TString addExtraTasks = VAR_EXTRATASKS_CONFIGMACRO;
  if (!addExtraTasks.IsNull()) gROOT->ProcessLineSync(TString::Format(".x %s",addExtraTasks.Data()));

  if (iESDfilter) {
    //  ESD filter task configuration.
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C");
    if (iMUONcopyAOD) {
      printf("Registering delta AOD file\n");
      mgr->RegisterExtraFile("AliAOD.Muons.root");
    }
    AliAnalysisTaskESDfilter *taskesdfilter = AddTaskESDFilter(useKFILTER, iMUONcopyAOD, kFALSE, kFALSE /*usePhysicsSelection*/,kFALSE,kTRUE,kTRUE,kTRUE,1500,VAR_MUONMCMODE,kFALSE,kFALSE,kTRUE,kFALSE); // others
    taskesdfilter->DisablePmdClusters();
    taskesdfilter->DisableCaloClusters();
    taskesdfilter->DisableCells();
    taskesdfilter->DisableCaloTrigger("PHOS");
    taskesdfilter->DisableCaloTrigger("EMCAL");
    taskesdfilter->SetPropagateTrackToEMCal(kFALSE);

    if ( 0 && VAR_USE_ITS_RECO ) /* 0 for the moment to get this macro running also with AliRoot <= .... */
    {
      AliAnalysisTaskESDMuonFilter* muFilter = mgr->GetTask("ESD Muon Filter");
      if ( !muFilter )
      {
        std::cout << "ERROR : got a NULL muFilter ! so I cannot ask to keep SPD tracklets !" << std::endl;
      }
      else
      {
        muFilter->SetWithSPDtracklets(kTRUE);
      }
    }
  }
}

//______________________________________________________________________________
Bool_t LoadAnalysisLibraries()
{
  // Load ParFiles.
  if (iMUONQA > 1) {
    if (!AliAnalysisAlien::SetupPar("PWGPPMUONlite.par")) return kFALSE;
  }
  if ((iMUONRefit > 1) || (iMUONRefitVtx > 1) || (iESDMCLabelAddition > 1) || (iMUONCDBConnect > 1)) {
    if (!AliAnalysisAlien::SetupPar("PWGmuondep.par")) return kFALSE;
  }
  if ((useMC && useTR && (iMUONPerformance > 1)) || (iMUONEfficiency > 1)) {
    if (!AliAnalysisAlien::SetupPar("PWGPPMUONdep.par")) return kFALSE;
  }
  ::Info("AODtrainsim.C::LoadAnalysisLibraries", "Load other libraries:   SUCCESS");
  return kTRUE;
}

//______________________________________________________________________________
TChain *CreateChain()
{
  // Create the input chain
  chain = new TChain("esdTree");
  if (gSystem->AccessPathName("AliESDs.root"))
    ::Error("AODtrainsim.C::CreateChain", "File: AliESDs.root not in ./data dir");
  else
    chain->Add("AliESDs.root");
  if (chain->GetNtrees()) return chain;
  return NULL;
}

//______________________________________________________________________________
void AODmerge()
{
  // Merging method. No staging and no terminate phase.
  TStopwatch timer;
  timer.Start();
  TString outputDir = "wn.xml";
  TString outputFiles = VAR_AOD_MERGE_FILES;
  TString mergeExcludes = "";
  TObjArray *list = outputFiles.Tokenize(",");
  TIter *iter = new TIter(list);
  TObjString *str;
  TString outputFile;
  Bool_t merged = kTRUE;
  while((str=(TObjString*)iter->Next())) {
    outputFile = str->GetString();
    // Skip already merged outputs
    if (!gSystem->AccessPathName(outputFile)) {
      printf("Output file <%s> found. Not merging again.",outputFile.Data());
      continue;
    }
    if (mergeExcludes.Contains(outputFile.Data())) continue;
    merged = AliAnalysisAlien::MergeOutput(outputFile, outputDir, 10, 0);
    if (!merged) {
      printf("ERROR: Cannot merge %s\n", outputFile.Data());
      return;
    }
  }
  // all outputs merged, validate
  ofstream out;
  out.open("outputs_valid_merge", ios::out);
  out.close();
  timer.Print();
}
