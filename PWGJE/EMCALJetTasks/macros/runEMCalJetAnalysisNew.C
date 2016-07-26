// runEMCalJetAnalysisNew.C

class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisManager;
class AliPhysicsSelectionTask;
class AliCentralitySelectionTask;
class AliEmcalSetupTask;
class AliAnalysisGrid;

void LoadMacros();
void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode);
AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers,
    const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, Int_t workerTTL, Bool_t isMC);

//______________________________________________________________________________
AliAnalysisManager* runEMCalJetAnalysisNew(
    const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
    const char   *cLocalFiles    = "fileLists/files_LHC11h_2_AOD145.txt",   // set the local list file
    UInt_t        iNumFiles      = 100,                                     // number of files analyzed locally
    UInt_t        iNumEvents     = 5000,                                    // number of events to be analyzed
    const char   *cRunPeriod     = "LHC11h",                                // set the run period
    UInt_t        kPhysSel       = AliVEvent::kAnyINT,                      // physics selection
    const char   *cTaskName      = "JetAna",                                // sets name of analysis manager
    const Bool_t  bDoChargedJets = kTRUE,
    const Bool_t  bDoFullJets    = kTRUE,
    const Bool_t  bDoMCJets      = kFALSE,
    const char   *cOCDBpath      = "uselocal",                              // change to "raw://" if running on the grid
    // 0 = only prepare the analysis manager but do not start the analysis
    // 1 = prepare the analysis manager and start the analysis
    // 2 = launch a grid analysis
    Int_t         iStartAnalysis = 1,
    const char   *cGridMode      = "test"
)
{
  TString sRunPeriod(cRunPeriod);
  sRunPeriod.ToLower();

  AliAnalysisTaskEmcal::BeamType iBeamType = AliAnalysisTaskEmcal::kpp;

  Bool_t bIsRun2 = kFALSE;
  if (sRunPeriod.Length() == 6 && sRunPeriod.BeginsWith("lhc15")) bIsRun2 = kTRUE;

  if (sRunPeriod == "lhc10h" || sRunPeriod == "lhc11h" || sRunPeriod == "lhc15o") {
    iBeamType = AliAnalysisTaskEmcal::kAA;
  }
  else if (sRunPeriod == "lhc12g" || sRunPeriod == "lhc13b" || sRunPeriod == "lhc13c" ||
      sRunPeriod == "lhc13d" || sRunPeriod == "lhc13e" || sRunPeriod == "lhc13f") {
    iBeamType = AliAnalysisTaskEmcal::kpA;
  }

  Double_t kGhostArea = 0.01;
  if (iBeamType != AliAnalysisTaskEmcal::kpp) kGhostArea = 0.005;

  Bool_t bSeparateEMCalDCal = bIsRun2;

  AliTrackContainer::SetDefTrackCutsPeriod(sRunPeriod);

  Printf("Default track cut period set to: %s", AliTrackContainer::GetDefTrackCutsPeriod().Data());

  const Bool_t   bDoTender            = kTRUE;
  const Bool_t   bDoHadCorr           = kTRUE;
  const Double_t kJetRadius           = 0.4;
  const Double_t kClusPtCut           = 0.30;
  const Double_t kTrackPtCut          = 0.15;
  const Double_t kJetPtCut            = 1.;
  const Double_t kHadCorrF            = 2.;
  
  const AliAnalysisTaskEmcalJetSpectraQA::EHistoType_t kHistoType = AliAnalysisTaskEmcalJetSpectraQA::kTHnSparse;

  enum eDataType { kAod, kEsd };

  eDataType iDataType;
  if (!strcmp(cDataType, "ESD")) {
    iDataType = kEsd;
  }
  else if (!strcmp(cDataType, "AOD")) {
    iDataType = kAod;
  }
  else {
    Printf("Incorrect data type option, check third argument of run macro.");
    Printf("datatype = AOD or ESD");
    return 0;
  }

  Printf("%s analysis chosen.", cDataType);

  TString sLocalFiles(cLocalFiles);
  if (sLocalFiles == "") {
    Printf("You need to provide the list of local files!");
    return 0;
  }
  Printf("Setting local analysis for %d files from list %s, max events = %d", iNumFiles, sLocalFiles.Data(), iNumEvents);

  LoadMacros();

  // Analysis manager
  AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);

  if (iDataType == kAod) {
    AliAODInputHandler* pAODHandler = AddAODHandler();
  }
  else {  
    AliESDInputHandler* pESDHandler = AddESDHandler();
  }

  // Physics selection task
  if (iDataType == kEsd) {
    AliPhysicsSelectionTask *pPhysSelTask = AddTaskPhysicsSelection();
  }

  // Centrality task
  if (iDataType == kEsd && iBeamType != AliAnalysisTaskEmcal::kpp) {
    AliCentralitySelectionTask *pCentralityTask = AddTaskCentrality(kTRUE);
    pCentralityTask->SelectCollisionCandidates(AliVEvent::kAny);
  }

  // Setup task
  if (bDoFullJets || iDataType == kEsd) {
    AliEmcalSetupTask *pSetupTask = AddTaskEmcalSetup();
    pSetupTask->SelectCollisionCandidates(AliVEvent::kAny);
    pSetupTask->SetOcdbPath(cOCDBpath);
  }

  if (bDoMCJets) {
    AliAnalysisTaskEmcalJetQA* pQATaskMC = AddTaskEmcalJetQA("mcparticles", "", "");
    pQATaskMC->SelectCollisionCandidates(kPhysSel);

    if (bDoChargedJets) {
      AliEmcalJetTask *pChJetTaskMC = AddTaskEmcalJet("mcparticles", "", AliJetContainer::antikt_algorithm, kJetRadius, AliJetContainer::kChargedJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 1., kFALSE, kFALSE);
      pChJetTaskMC->SelectCollisionCandidates(kPhysSel);
    }

    if (bDoFullJets) {
      AliEmcalJetTask *pChJetTaskMC = AddTaskEmcalJet("mcparticles", "", AliJetContainer::antikt_algorithm, kJetRadius, AliJetContainer::kFullJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 1., kFALSE, kFALSE);
      pChJetTaskMC->SelectCollisionCandidates(kPhysSel);
    }
  }

  if (bDoFullJets && bDoTender) {
    // QA task
    AliAnalysisTaskEmcalJetQA* pQATaskBefore = AddTaskEmcalJetQA("", "usedefault", "usedefault", "BeforeTender");
    pQATaskBefore->GetClusterContainer(0)->SetClusECut(0.15);
    pQATaskBefore->GetClusterContainer(0)->SetClusPtCut(0.);
    pQATaskBefore->SetSeparateEMCalDCal(bSeparateEMCalDCal);
    pQATaskBefore->SetHistoBins(300, 0, 150);
    pQATaskBefore->SelectCollisionCandidates(kPhysSel);

    const char *cPass        = 0;
    Bool_t   bDistBC         = kFALSE; //switch for recalculation cluster position from bad channel
    Bool_t   bRecalibClus    = kFALSE;
    Bool_t   bRecalcClusPos  = kFALSE;
    Bool_t   bNonLinearCorr  = kFALSE;
    Bool_t   bRemExoticCell  = kFALSE;
    Bool_t   bRemExoticClus  = kFALSE;
    Bool_t   bFidRegion      = kFALSE;
    Bool_t   bCalibEnergy    = kTRUE;
    Bool_t   bCalibTime      = kTRUE;
    Bool_t   bRemBC          = kTRUE;
    UInt_t   iNonLinFunct    = AliEMCALRecoUtils::kNoCorrection;
    Bool_t   bReclusterize   = kFALSE;
    Float_t  fSeedThresh     = 0.1;      // 100 MeV
    Float_t  fCellThresh     = 0.05;     // 50 MeV
    UInt_t   iClusterizer    = AliEMCALRecParam::kClusterizerv2;
    Bool_t   bTrackMatch     = kFALSE;
    Bool_t   bUpdateCellOnly = kTRUE;
    Float_t  fEMCtimeMin     = -50e-6;
    Float_t  fEMCtimeMax     =  50e-6;
    Float_t  fEMCtimeCut     =  1e6;
    if (sRunPeriod == "lhc11h") {
      fEMCtimeMin = -50e-9;
      fEMCtimeMax = 100e-9;
    }

    AliAnalysisTaskSE *pTenderTask = AddTaskEMCALTender(bDistBC, bRecalibClus, bRecalcClusPos, bNonLinearCorr, bRemExoticCell, bRemExoticClus,
        bFidRegion, bCalibEnergy, bCalibTime, bRemBC, iNonLinFunct, bReclusterize, fSeedThresh,
        fCellThresh, iClusterizer, bTrackMatch, bUpdateCellOnly, fEMCtimeMin, fEMCtimeMax, fEMCtimeCut, cPass);
    pTenderTask->SelectCollisionCandidates(kPhysSel);

    AliAnalysisTaskEMCALClusterizeFast *pClusterizerTask = AddTaskClusterizerFast("ClusterizerFast", "", "", iClusterizer,
        0.05, 0.1, fEMCtimeMin, fEMCtimeMax, fEMCtimeCut,
        kFALSE, kFALSE, AliAnalysisTaskEMCALClusterizeFast::kFEEData);
    pClusterizerTask->SelectCollisionCandidates(kPhysSel);

    bRemExoticClus  = kTRUE;
    iNonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrected;

    AliEmcalClusterMaker *pClusterMakerTask = AddTaskEmcalClusterMaker(iNonLinFunct, bRemExoticClus, "usedefault", "", 0., kFALSE);
    pClusterMakerTask->GetClusterContainer(0)->SetClusPtCut(0.);
    pClusterMakerTask->GetClusterContainer(0)->SetClusECut(0.);
    pClusterMakerTask->SelectCollisionCandidates(kPhysSel);
  }

  if (bDoFullJets && bDoHadCorr) {
    // Cluster-track matcher task
    AliEmcalClusTrackMatcherTask *pMatcherTask = AddTaskEmcalClusTrackMatcher("usedefault", "usedefault", 0.1, kFALSE, kTRUE, kTRUE, kFALSE);
    pMatcherTask->SelectCollisionCandidates(kPhysSel);
    pMatcherTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
    pMatcherTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.15);
    pMatcherTask->GetClusterContainer(0)->SetClusECut(0.);
    pMatcherTask->GetClusterContainer(0)->SetClusPtCut(0.);

    if (iDataType == kEsd) {
      pMatcherTask->SetDoPropagation(kTRUE);
    }

    // Hadronic correction task
    AliHadCorrTask *pHadCorrTask = AddTaskHadCorr("usedefault", "usedefault", "",
        kHadCorrF, 0.15, 0.030, 0.015, 0, kTRUE, kFALSE);
    pHadCorrTask->SelectCollisionCandidates(kPhysSel);
    pHadCorrTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.15);
    pHadCorrTask->GetClusterContainer(0)->SetClusECut(0);
    pHadCorrTask->GetClusterContainer(0)->SetClusPtCut(0.);
  }

  // QA task
  AliAnalysisTaskEmcalJetQA *pQATask = 0;
  if (bDoFullJets) {
    pQATask = AddTaskEmcalJetQA("usedefault", "usedefault", "usedefault");
    pQATask->GetClusterContainer(0)->SetClusECut(0.);
    pQATask->GetClusterContainer(0)->SetClusPtCut(0.);
    pQATask->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
    pQATask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    pQATask->SetSeparateEMCalDCal(bSeparateEMCalDCal);
  }
  else {
    pQATask = AddTaskEmcalJetQA("usedefault", "", "");
  }
  pQATask->GetParticleContainer(0)->SetParticlePtCut(0.15);
  pQATask->SetHistoBins(300, 0, 150);
  pQATask->SelectCollisionCandidates(kPhysSel);

  // Background
  TString sRhoChName;
  TString sRhoFuName;
  if (iBeamType != AliAnalysisTaskEmcal::kpp) {
    sRhoChName = "Rho";
    sRhoFuName = "Rho_Scaled";

    AliEmcalJetTask *pKtChJetTask = AddTaskEmcalJet("usedefault", "", AliJetContainer::kt_algorithm, kJetRadius, AliJetContainer::kChargedJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 0., kFALSE, kFALSE);
    pKtChJetTask->SelectCollisionCandidates(kPhysSel);

    AliAnalysisTaskRho* pRhoTask = AddTaskRhoNew("usedefault", "usedefault", sRhoChName, kJetRadius);
    pRhoTask->SetExcludeLeadJets(2);
    pRhoTask->SelectCollisionCandidates(kPhysSel);

    if (bDoFullJets) {
      TString sFuncPath = "alien:///alice/cern.ch/user/s/saiola/LHC11h_ScaleFactorFunctions.root";
      TString sFuncName = "LHC11h_HadCorr20_ClustersV2";
      pRhoTask->LoadRhoFunction(sFuncPath, sFuncName);
    }
  }

  // Charged jet analysis
  if (bDoChargedJets) {
    AliEmcalJetTask *pChJetTask = AddTaskEmcalJet("usedefault", "", AliJetContainer::antikt_algorithm, kJetRadius, AliJetContainer::kChargedJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 1., kFALSE, kFALSE);
    pChJetTask->SelectCollisionCandidates(kPhysSel);
  }

  // Full jet analysis
  if (bDoFullJets) {
    AliEmcalJetTask *pFuJetTask = AddTaskEmcalJet("usedefault", "usedefault", AliJetContainer::antikt_algorithm, kJetRadius, AliJetContainer::kFullJet, 0.15, 0.30, kGhostArea, AliJetContainer::pt_scheme, "Jet", 1., kFALSE, kFALSE);
    pFuJetTask->SelectCollisionCandidates(kPhysSel);
  }

  AliAnalysisTaskEmcalJetSpectraQA *pSpectraTask = AddTaskEmcalJetSpectraQA("usedefault", "usedefault");

  if (bDoFullJets) {
    AliJetContainer* jetCont = pSpectraTask->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, kJetRadius, AliJetContainer::kEMCALfid, "Jet");
    if (iBeamType != AliAnalysisTaskEmcal::kpp) {
      jetCont->SetRhoName(sRhoFuName);
      jetCont->SetPercAreaCut(0.6);
    }
  }

  if (bDoChargedJets) {
    AliJetContainer* jetCont = pSpectraTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, kJetRadius, AliJetContainer::kTPCfid, "Jet");
    if (iBeamType != AliAnalysisTaskEmcal::kpp) {
      jetCont->SetRhoName(sRhoChName);
      jetCont->SetPercAreaCut(0.6);
    }
  }

  if (bDoMCJets) {
    AliMCParticleContainer* MCpartCont = pSpectraTask->AddMCParticleContainer("mcparticles");

    if (bDoChargedJets) {
      pSpectraTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, kJetRadius, AliJetContainer::kTPCfid, MCpartCont, 0);
    }

    if (bDoFullJets) {
      pSpectraTask->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, kJetRadius, AliJetContainer::kTPCfid, MCpartCont, 0);
    }
  }

  pSpectraTask->SetNLeadingJets(1);
  pSpectraTask->SelectCollisionCandidates(kPhysSel);
  pSpectraTask->SetHistoType(kHistoType);

  TObjArray *pTopTasks = pMgr->GetTasks();
  for (Int_t i = 0; i < pTopTasks->GetEntries(); ++i) {
    AliAnalysisTaskSE *pTask = dynamic_cast<AliAnalysisTaskSE*>(pTopTasks->At(i));
    if (!pTask) continue;
    if (pTask->InheritsFrom("AliAnalysisTaskEmcal")) {
      AliAnalysisTaskEmcal *pTaskEmcal = static_cast<AliAnalysisTaskEmcal*>(pTask);
      Printf("Setting beam type %d for task %s", iBeamType, pTaskEmcal->GetName());
      pTaskEmcal->SetForceBeamType(iBeamType);
    }
  }

  if (!pMgr->InitAnalysis()) return 0;
  pMgr->PrintStatus();
    
  pMgr->SetUseProgressBar(kTRUE, 250);
  
  TFile *pOutFile = new TFile("train.root","RECREATE");
  pOutFile->cd();
  pMgr->Write();
  pOutFile->Close();
  delete pOutFile;

  if (iStartAnalysis == 1) { // start local analysis
    TChain* pChain = 0;
    if (iDataType == kAod) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
      pChain = CreateAODChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
    }
    else {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      pChain = CreateESDChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
    }

    // start analysis
    Printf("Starting Analysis...");
    pMgr->StartAnalysis("local", pChain, iNumEvents);
  }
  else if (iStartAnalysis == 2) {  // start grid analysis
    StartGridAnalysis(pMgr, cTaskName, cGridMode);
  }

  return pMgr;
}

void LoadMacros()
{
  // Aliroot macros
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALTender.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskClusterizerFast.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C"); 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetSpectraQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C");
}

void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode)
{
  Int_t maxFilesPerWorker = 4;
  Int_t workerTTL = 7200;
  const char* runNumbers = "180720";
  const char* pattern = "pass2/AOD/*/AliAOD.root";
  const char* gridDir = "/alice/data/2012/LHC12c";
  const char* additionalCXXs = "";
  const char* additionalHs = "";

  AliAnalysisGrid *plugin = CreateAlienHandler(uniqueName, gridDir, cGridMode, runNumbers, pattern, additionalCXXs, additionalHs, maxFilesPerWorker, workerTTL, kFALSE);
  pMgr->SetGridHandler(plugin);

  // start analysis
   Printf("Starting GRID Analysis...");
   pMgr->SetDebugLevel(0);
   pMgr->StartAnalysis("grid");
}

AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers,
    const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, Int_t workerTTL, Bool_t isMC)
{
  TDatime currentTime;
  TString tmpName(uniqueName);

  // Only add current date and time when not in terminate mode! In this case the exact name has to be supplied by the user
  if (strcmp(gridMode, "terminate")) {
    tmpName += "_";
    tmpName += currentTime.GetDate();
    tmpName += "_";
    tmpName += currentTime.GetTime();
  }

  TString macroName("");
  TString execName("");
  TString jdlName("");
  macroName = Form("%s.C", tmpName.Data());
  execName = Form("%s.sh", tmpName.Data());
  jdlName = Form("%s.jdl", tmpName.Data());

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  plugin->SetRunMode(gridMode);

  // Here you can set the (Ali)PHYSICS version you want to use
  plugin->SetAliPhysicsVersion("vAN-20160203-1");

  plugin->SetGridDataDir(gridDir); // e.g. "/alice/sim/LHC10a6"
  plugin->SetDataPattern(pattern); //dir structure in run directory

  if (!isMC) plugin->SetRunPrefix("000");

  plugin->AddRunList(runNumbers);

  plugin->SetGridWorkingDir(Form("work/%s",tmpName.Data()));
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output

  plugin->SetAnalysisSource(additionalCode.Data());

  plugin->SetDefaultOutputs(kTRUE);
  plugin->SetAnalysisMacro(macroName.Data());
  plugin->SetSplitMaxInputFileNumber(maxFilesPerWorker);
  plugin->SetExecutable(execName.Data());
  plugin->SetTTL(workerTTL);
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName(jdlName.Data());
  plugin->SetPrice(1);
  plugin->SetSplitMode("se");

  // merging via jdl
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(2);

  return plugin;
}
