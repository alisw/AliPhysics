// runEMCalJetAnalysisOld.C

class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisManager;
class AliEmcalPhysicsSelectionTask;
class AliCentralitySelectionTask;
class AliEmcalSetupTask;

void LoadLibs();
void LoadMacros();

//______________________________________________________________________________
void runEMCalJetAnalysisOld(
    const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
    const char   *cLocalFiles    = "fileLists/files_LHC11h_2_AOD145.txt",   // set the local list file
    UInt_t        iNumFiles      = 100,                                     // number of files analyzed locally
    UInt_t        iNumEvents     = 5000,                                    // number of events to be analyzed
    const char   *cRunPeriod     = "LHC11h",                                // set the run period
    const char   *cTaskName      = "JetAna"                                 // sets name of analysis manager
)
{
  const Bool_t   bDoEmcal             = kTRUE;
  const Bool_t   bUseNewCode          = kTRUE;
  const Bool_t   bDoChargedJets       = kTRUE;
  const Bool_t   bDoFullJets          = kTRUE;
  const Bool_t   bDoTender            = kTRUE;
  const Bool_t   bDoReclusterize      = kTRUE;
  const Bool_t   bDoHadCorr           = kTRUE;
  const Double_t kJetRadius           = 0.4;
  const Double_t kClusPtCut           = 0.30;
  const Double_t kTrackPtCut          = 0.15;
  const Double_t kJetPtCut            = 1.;
  const Double_t kGhostArea           = 0.005;
  const Double_t kHadCorrF            = 2.;
  const Int_t    kHistoType           = 1;
  const UInt_t   kClusterizerType     = AliEMCALRecParam::kClusterizerv2;
  const Bool_t   bForcePP             = kFALSE;

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
    return;
  }

  Printf("%s analysis chosen.", cDataType);

  TString sLocalFiles(cLocalFiles);
  if (sLocalFiles == "") {
    Printf("You need to provide the list of local files!");
    return;
  }
  Printf("Setting local analysis for %d files from list %s, max events = %d", iNumFiles, sLocalFiles.Data(), iNumEvents);

  LoadLibs();
  LoadMacros();

  TString sTracksName("PicoTracks");
  TString sClusName("EmcCaloClusters");

  TString sCellName;
  TString sOrigClusName;
  TString sCorrClusName;

  TString sChJetsName;
  TString sFuJetsName;

  if (iDataType == kAod) {
    sCellName = "emcalCells";
    sOrigClusName = "caloClusters";
  }
  else {
    sCellName = "EMCALCells";
    sOrigClusName = "CaloClusters";
  }

  // AliEmcalPhysicsSelection::kEmcalOk, AliEmcalPhysicsSelection::kEmcalH,
  // AliVEvent::kINT7, AliVEvent::kMB, AliVEvent::kCentral, AliVEvent::kSemiCentral,
  // AliVEvent::kEMCEGA, AliVEvent::kEMCEJE
  UInt_t kPrePhysSel = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  UInt_t kPhysSel = AliEmcalPhysicsSelection::kEmcalOk; 

  // Analysis manager
  AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);

  if (iDataType == kAod) {
    AliAODInputHandler* pAODHandler = AddAODHandler();
  }
  else {   //ESD or skimmed ESD
    AliESDInputHandler* pESDHandler = AddESDHandler();
  }

  // Physics selection task
  AliEmcalPhysicsSelectionTask *pPhysSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE,
      kPrePhysSel,
      5, 5, 10, kTRUE);

  // Centrality task
  if (iDataType == kEsd) {
    AliCentralitySelectionTask *pCentralityTask = AddTaskCentrality(kTRUE);
    pCentralityTask->SelectCollisionCandidates(kPhysSel);
  }

  // Setup task
  if (bDoEmcal) {
    AliEmcalSetupTask *pSetupTask = AddTaskEmcalSetup();
  }

  if (iDataType == kEsd) {
    // Hybrid tracks maker for ESD
    TString trackCuts("Hybrid_");
    trackCuts += cRunPeriod;
    AliEmcalEsdTrackFilterTask *pHybTask = AddTaskEmcalEsdTrackFilter("HybridTracks", trackCuts);
    pHybTask->SelectCollisionCandidates(kPhysSel);
    pHybTask->SetDoPropagation(bDoHadCorr);
  }
  else if (iDataType == kAod) {
    // Hybrid tracks maker for AOD
    AliEmcalAodTrackFilterTask *pHybTask = AddTaskEmcalAodTrackFilter("HybridTracks", "tracks", cRunPeriod);
    pHybTask->SelectCollisionCandidates(kPhysSel);
  }

  AliEmcalPicoTrackMaker *pPicoTrackTask = AddTaskEmcalPicoTrackMaker(sTracksName, "HybridTracks");
  pPicoTrackTask->SelectCollisionCandidates(kPhysSel);

  if (bDoTender || bDoReclusterize) {
    // QA task
    AliAnalysisTaskSAQA *pQATaskBefore = AddTaskSAQA("", sOrigClusName, sCellName, "", "",
        0, 0, 0, 0., 0., "TPC", "BeforeTender");
    pQATaskBefore->GetClusterContainer(0)->SetClusECut(0.15);
    pQATaskBefore->GetClusterContainer(0)->SetClusPtCut(0.);
    pQATaskBefore->SetHistoBins(200, 0, 30);
    pQATaskBefore->SelectCollisionCandidates(kPhysSel);
  }

  // Tender Supplies
  if (bDoTender) {
    Bool_t  bDistBC         = kFALSE; //switch for recalculation cluster position from bad channel 
    Bool_t  bRecalibClus    = kFALSE;
    Bool_t  bRecalcClusPos  = kFALSE;
    Bool_t  bNonLinearCorr  = kFALSE;
    Bool_t  bRemExoticCell  = kFALSE;
    Bool_t  bRemExoticClus  = kFALSE;
    Bool_t  bFidRegion      = kFALSE;
    Bool_t  bCalibEnergy    = kTRUE;
    Bool_t  bCalibTime      = kTRUE;
    Bool_t  bRemBC          = kTRUE;
    UInt_t  iNonLinFunct    = 0;
    Bool_t  bReclusterize   = kFALSE;
    Float_t fSeedThresh     = 0.1;      // 100 MeV
    Float_t fCellThresh     = 0.05;     // 50 MeV 
    UInt_t  iClusterizer    = 0;
    Bool_t  bTrackMatch     = kFALSE;
    Bool_t  bUpdateCellOnly = kFALSE;
    Float_t fTimeMin        = -50e-6;   // minimum time of physical signal in a cell/digit
    Float_t fTimeMax        =  50e-6;   // maximum time of physical signal in a cell/digit
    Float_t fTimeCut        =  25e-6;
    const char *cPass       = 0;

    AliAnalysisTaskSE *pTenderTask = AddTaskEMCALTender(bDistBC, bRecalibClus, bRecalcClusPos, bNonLinearCorr, bRemExoticCell, bRemExoticClus,
        bFidRegion, bCalibEnergy, bCalibTime, bRemBC, iNonLinFunct, bReclusterize, fSeedThresh,
        fCellThresh, iClusterizer, bTrackMatch, bUpdateCellOnly, fTimeMin, fTimeMax, fTimeCut, cPass);
    pTenderTask->SelectCollisionCandidates(kPhysSel);
  }

  if (bDoReclusterize) {

    const Double_t kEMCtimeMin          = -50e-6;
    const Double_t kEMCtimeMax          = 100e-6;
    const Double_t kEMCtimeCut          =  75e-6;

    AliAnalysisTaskEMCALClusterizeFast *pClusterizerTask = AddTaskClusterizerFast("ClusterizerFast", "", "", kClusterizerType, 
        0.05, 0.1, kEMCtimeMin, kEMCtimeMax, kEMCtimeCut,
        kFALSE, kFALSE, AliAnalysisTaskEMCALClusterizeFast::kFEEData);

    pClusterizerTask->SelectCollisionCandidates(kPhysSel);

    AliEmcalClusterMaker *pClusterMakerTask = AddTaskEmcalClusterMaker(AliEMCALRecoUtils::kBeamTestCorrected, kTRUE, 0, sClusName, 0., kTRUE);
    pClusterMakerTask->GetClusterContainer(0)->SetClusPtCut(0.);
    pClusterMakerTask->GetClusterContainer(0)->SetClusECut(0.);
    pClusterMakerTask->SelectCollisionCandidates(kPhysSel);
  }

  if (bDoTender || bDoReclusterize) {
    // QA task

    AliAnalysisTaskSAQA *pQATaskAfter = AddTaskSAQA("", sOrigClusName, sCellName, "", "",
        0, 0, 0, 0., 0., "TPC", "AfterTender");
    pQATaskAfter->GetClusterContainer(0)->SetClusECut(0.);
    pQATaskAfter->GetClusterContainer(0)->SetClusPtCut(0.);
    pQATaskAfter->GetClusterContainer(0)->SetExoticCut(kFALSE);
    pQATaskAfter->SetHistoBins(200, 0, 30);
    pQATaskAfter->SelectCollisionCandidates(kPhysSel);

    AliAnalysisTaskSAQA *pQATaskAfterMaker = AddTaskSAQA("", sClusName, "", "", "",
        0, 0, 0, 0., 0., "TPC", "AfterClusterMaker");
    pQATaskAfterMaker->GetClusterContainer(0)->SetClusECut(0.15);
    pQATaskAfterMaker->GetClusterContainer(0)->SetClusPtCut(0.);
    pQATaskAfterMaker->SetHistoBins(200, 0, 30);
    pQATaskAfterMaker->SelectCollisionCandidates(kPhysSel);

  }

  if (bDoHadCorr) {    
    sCorrClusName = "CaloClustersCorr";

    TString sEmcalTracksName("EmcalTracks_");
    TString sEmcalClusName("EmcalClusters_");

    sEmcalTracksName += sTracksName;
    sEmcalClusName += sClusName;

    // EmcalParticles maker
    if (!bUseNewCode) {
      AliEmcalParticleMaker *pPartMakerTask = AddTaskEmcalParticleMaker(sTracksName, sClusName, sEmcalTracksName, sEmcalClusName, "AliEmcalParticleMaker");
      pPartMakerTask->GetClusterContainer(0)->SetClusECut(0.15);
      pPartMakerTask->GetClusterContainer(0)->SetClusPtCut(0.);
      pPartMakerTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
      pPartMakerTask->SelectCollisionCandidates(kPhysSel);
    }

    // Cluster-track matcher task
    if (bUseNewCode) {
      AliEmcalClusTrackMatcherTask *pMatcherTask = AddTaskEmcalClusTrackMatcher(sTracksName, sClusName, 0.1, kTRUE, kTRUE, kTRUE, kTRUE);
      pMatcherTask->SelectCollisionCandidates(kPhysSel);
      pMatcherTask->GetClusterContainer(0)->SetClusECut(0.15);
      pMatcherTask->GetClusterContainer(0)->SetClusPtCut(0.);
      pMatcherTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
    }
    else {
      AliEmcalClusTrackMatcherTask *pMatcherTask = AddTaskEmcalClusTrackMatcher(sEmcalTracksName, sEmcalClusName, 0.1, kFALSE, kTRUE);
      pMatcherTask->SelectCollisionCandidates(kPhysSel);
      pMatcherTask->GetParticleContainer(0)->SetParticlePtCut(0.);
      pMatcherTask->GetParticleContainer(1)->SetParticlePtCut(0.);
      pMatcherTask->SetModifyObjs(kTRUE);
    }

    // Hadronic correction task
    if (bUseNewCode) {
      AliHadCorrTask *pHadCorrTask = AddTaskHadCorr(sTracksName, sClusName, sCorrClusName, 
          kHadCorrF, 0.15, 0.030, 0.015, 0, kTRUE, kTRUE);
      pHadCorrTask->SetHistoBins(150,0,150);
      pHadCorrTask->SelectCollisionCandidates(kPhysSel);
    }
    else {
      AliHadCorrTask *pHadCorrTask = AddTaskHadCorr(sEmcalTracksName, sEmcalClusName, sCorrClusName, 
          kHadCorrF, 0.15, 0.030, 0.015, 0, kTRUE, kTRUE);
      pHadCorrTask->SetHistoBins(150,0,150);
      pHadCorrTask->SelectCollisionCandidates(kPhysSel); 
    }
  }
  else {
    sCorrClusName = sClusName;
  }

  // QA task
  AliAnalysisTaskSAQA *pQATask = AddTaskSAQA(sTracksName, sCorrClusName, sCellName, "", "", 0., 0, 0, 0., 0., "TPC");
  pQATask->GetClusterContainer(0)->SetClusECut(0.30);
  pQATask->GetClusterContainer(0)->SetClusPtCut(0.);
  pQATask->GetParticleContainer(0)->SetParticlePtCut(0.15);
  pQATask->GetParticleContainer(0)->SetClassName("AliPicoTrack");
  pQATask->SelectCollisionCandidates(kPhysSel);
  pQATask->SetHistoBins(200, 0, 30);

  // Charged jet analysis
  if (bDoChargedJets) {
    AliEmcalJetTask *pChJetTask = AddTaskEmcalJet(sTracksName, "", 1, kJetRadius, 1, kTrackPtCut, kClusPtCut, kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, kFALSE);
    pChJetTask->SelectCollisionCandidates(kPhysSel);
    sChJetsName = pChJetTask->GetName();

    AliAnalysisTaskSAJF *pSpectraChTask = AddTaskSAJF(sTracksName, "", sChJetsName, "",  kJetRadius, kJetPtCut, 0., "TPC");
    pSpectraChTask->SetNLeadingJets(1);
    pSpectraChTask->SelectCollisionCandidates(kPhysSel);
    pSpectraChTask->SetHistoType(kHistoType);
  }

  // Full jet analysis
  if (bDoFullJets) {
    AliEmcalJetTask *pFuJetTask = AddTaskEmcalJet(sTracksName, sCorrClusName, 1, kJetRadius, 0, kTrackPtCut, kClusPtCut, kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, kFALSE);
    pFuJetTask->SelectCollisionCandidates(kPhysSel);   
    sFuJetsName = pFuJetTask->GetName();

    AliAnalysisTaskSAJF *pSpectraFuTask = AddTaskSAJF(sTracksName, sCorrClusName, sFuJetsName, "", kJetRadius, kJetPtCut, 0., "EMCAL"); 
    pSpectraFuTask->SetNLeadingJets(1);
    pSpectraFuTask->SelectCollisionCandidates(kPhysSel);
    pSpectraFuTask->SetHistoType(kHistoType);
  }

  TObjArray *pTopTasks = pMgr->GetTasks();
  for (Int_t i = 0; i < pTopTasks->GetEntries(); ++i) {
    AliAnalysisTaskSE *pTask = dynamic_cast<AliAnalysisTaskSE*>(pTopTasks->At(i));
    if (!pTask) continue;
    if (pTask->InheritsFrom("AliAnalysisTaskEmcal")) {
      AliAnalysisTaskEmcal *pTaskEmcal = static_cast<AliAnalysisTaskEmcal*>(pTask);
      if (bForcePP) {
        Printf("Setting beam type for task %s", pTaskEmcal->GetName());
        pTaskEmcal->SetForceBeamType(0);
      }
    }
  }

  if (!pMgr->InitAnalysis()) return;
  pMgr->PrintStatus();

  TChain* pChain = 0;
  if (iDataType == kAod) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
    pChain = CreateAODChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE, "");
  }
  else { 
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
    pChain = CreateESDChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
  }

  // start analysis
  Printf("Starting Analysis...");
  pMgr->SetUseProgressBar(1, 250);
  //pMgr->SetDebugLevel(1);

  // To have more debug info
  //pMgr->AddClassDebug("AliEmcalClusTrackMatcherTask", AliLog::kDebug+100);

  TFile *pOutFile = new TFile("train.root","RECREATE");
  pOutFile->cd();
  pMgr->Write();
  pOutFile->Close();
  delete pOutFile;

  pMgr->StartAnalysis("local", pChain, iNumEvents);
}

//______________________________________________________________________________
void LoadLibs()
{
  // load fastjet libraries 3.x
  gSystem->Load("libCGAL");
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libsiscone_spherical");
  gSystem->Load("$FASTJET/lib/libfastjetplugins");
  gSystem->Load("$FASTJET/lib/libfastjetcontribfragile");
}

void LoadMacros()
{
  // Aliroot macros
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalSetup.C");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskTrackingQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALTender.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskClusterizerFast.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C"); 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAJF.C");
}
