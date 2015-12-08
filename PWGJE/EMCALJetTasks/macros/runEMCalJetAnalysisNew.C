// runEMCalJetAnalysisNew.C

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
void runEMCalJetAnalysisNew(
    const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
    const char   *cLocalFiles    = "fileLists/files_LHC11h_2_AOD145.txt",   // set the local list file
    UInt_t        iNumFiles      = 100,                                     // number of files analyzed locally
    UInt_t        iNumEvents     = 5000,                                    // number of events to be analyzed
    const char   *cRunPeriod     = "LHC11h",                                // set the run period
    const char   *cTaskName      = "JetAna"                                 // sets name of analysis manager
)
{
  const Bool_t   bDoChargedJets       = kTRUE;
  const Bool_t   bDoFullJets          = kTRUE;
  const Bool_t   bDoTender            = kTRUE;
  const Bool_t   bDoHadCorr           = kTRUE;
  const Double_t kJetRadius           = 0.4;
  const Double_t kClusPtCut           = 0.30;
  const Double_t kTrackPtCut          = 0.15;
  const Double_t kJetPtCut            = 1.;
  const Double_t kGhostArea           = 0.005;
  const Double_t kHadCorrF            = 2.;
  const Int_t    kHistoType           = 1;
  
  enum eDataType { kAod, kEsd };

  eDataType iDataType;
  if (!strcmp(cDataType, "ESD")) {
    iDataType = kEsd;
    Printf("For the moment, this macro is only available for AOD analysis!");
    return;
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
  else {  
    AliESDInputHandler* pESDHandler = AddESDHandler();
  }

  if (0) {
    AliAODHandler* pAODOutHandler = AddAODOutputHandler();
  }

  // Physics selection task
  if (1) {
    AliEmcalPhysicsSelectionTask *pPhysSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE,
                                                                              kPrePhysSel,
                                                                              5, 5, 10, kTRUE);
  }

  // Centrality task
  if (iDataType == kEsd) {
    AliCentralitySelectionTask *pCentralityTask = AddTaskCentrality(kTRUE);
    pCentralityTask->SelectCollisionCandidates(kPhysSel);
  }
    
  // Setup task
  if (bDoFullJets || bDoTender || bDoHadCorr) {
    AliEmcalSetupTask *pSetupTask = AddTaskEmcalSetup();
  }

  TString sTracksName;
  TString sClusName;
  TString sCellName;
  TString sChJetsName;
  TString sFuJetsName;

  if (iDataType == kAod) {
    sCellName = "emcalCells";
    sTracksName = "tracks";
    sClusName = "caloClusters";
  }
  else {
    sCellName = "EMCALCells";
    sTracksName = "Tracks";
    sClusName = "CaloClusters";
  }

  AliEmcalTrackPropagatorTask* pTrackProp = AddTaskEmcalTrackPropagator();
  pTrackProp->SelectCollisionCandidates(kPhysSel);

  if (bDoTender) {
    // QA task
    if (1) {
      AliAnalysisTaskSAQA *pQATaskBefore = AddTaskSAQA("", sClusName, sCellName, "", "",
                                                       0, 0, 0, 0., 0., "TPC", "BeforeTender");
      pQATaskBefore->GetClusterContainer(0)->SetClusECut(0.15);
      pQATaskBefore->GetClusterContainer(0)->SetClusPtCut(0.);
      pQATaskBefore->SetHistoBins(200, 0, 30);
      pQATaskBefore->SelectCollisionCandidates(kPhysSel);
    }

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
    Bool_t   bUpdateCellOnly = kFALSE;
    Float_t  fEMCtimeMin     = -50e-6;
    Float_t  fEMCtimeMax     =  50e-6;
    Float_t  fEMCtimeCut     =  1e6;
    if (strcmp(cRunPeriod, "LHC11h") == 0) {
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

    AliEmcalClusterMaker *pClusterMakerTask = AddTaskEmcalClusterMaker(iNonLinFunct, bRemExoticClus, 0, "", 0., kTRUE);
    pClusterMakerTask->GetClusterContainer(0)->SetClusPtCut(0.);
    pClusterMakerTask->GetClusterContainer(0)->SetClusECut(0.);
    pClusterMakerTask->SelectCollisionCandidates(kPhysSel);
  }
  
  if (bDoHadCorr) {    
    // Cluster-track matcher task
    AliEmcalClusTrackMatcherTask *pMatcherTask = AddTaskEmcalClusTrackMatcher(sTracksName, sClusName, 0.1, kFALSE, kTRUE, kTRUE, kTRUE);
    pMatcherTask->SelectCollisionCandidates(kPhysSel);
    pMatcherTask->GetParticleContainer(0)->SetClassName("AliAODTrack");
    pMatcherTask->GetParticleContainer(0)->SetFilterHybridTracks(kTRUE);
    pMatcherTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
    pMatcherTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.15);
    pMatcherTask->GetClusterContainer(0)->SetClusECut(0.);
    pMatcherTask->GetClusterContainer(0)->SetClusPtCut(0.);

    if (iDataType == kEsd) {
      pMatcherTask->SetDoPropagation(kTRUE);
    }

    // Hadronic correction task
    AliHadCorrTask *pHadCorrTask = AddTaskHadCorr(sTracksName, sClusName, "", 
                                                  kHadCorrF, 0.15, 0.030, 0.015, 0, kTRUE, kTRUE);
    pHadCorrTask->SelectCollisionCandidates(kPhysSel);
    pHadCorrTask->GetParticleContainer(0)->SetClassName("AliAODTrack");
    pHadCorrTask->GetParticleContainer(0)->SetFilterHybridTracks(kTRUE);
    pHadCorrTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.15);
    pHadCorrTask->GetClusterContainer(0)->SetClusECut(0);
    pHadCorrTask->GetClusterContainer(0)->SetClusPtCut(0.);
    pHadCorrTask->SetHistoBins(200, 0, 30);
  }

  if (1) {
    AliAnalysisTaskSAQA *pQATaskAfter = AddTaskSAQA("", sClusName, sCellName, "", "", 0, 0, 0, 0., 0., "TPC", "AliAnalysisTaskSAQA_AfterTender");
    pQATaskAfter->GetClusterContainer(0)->SetClusECut(0.15);
    pQATaskAfter->GetClusterContainer(0)->SetClusPtCut(0.);
    pQATaskAfter->GetClusterContainer(0)->SetExoticCut(kFALSE);
    pQATaskAfter->SetHistoBins(200, 0, 30);
    pQATaskAfter->SelectCollisionCandidates(kPhysSel);
  }

  if (1) {
    AliAnalysisTaskSAQA *pQATaskAfterMaker = AddTaskSAQA("", sClusName, "", "", "", 0, 0, 0, 0., 0., "TPC", "AliAnalysisTaskSAQA_AfterClusterMaker");
    pQATaskAfterMaker->GetClusterContainer(0)->SetClusECut(0.);
    pQATaskAfterMaker->GetClusterContainer(0)->SetClusPtCut(0.);
    pQATaskAfterMaker->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.15);
    pQATaskAfterMaker->SelectCollisionCandidates(kPhysSel);
    pQATaskAfterMaker->SetHistoBins(200, 0, 30);
    pQATaskAfterMaker->SetDefaultClusterEnergy(AliVCluster::kNonLinCorr);
  }
  
  // QA task
  AliAnalysisTaskSAQA *pQATask = AddTaskSAQA(sTracksName, sClusName, sCellName, "", "", 0, 0, 0, 0., 0., "TPC");
  pQATask->GetClusterContainer(0)->SetClusECut(0.);
  pQATask->GetClusterContainer(0)->SetClusPtCut(0.);
  pQATask->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
  pQATask->GetParticleContainer(0)->SetClassName("AliAODTrack");
  pQATask->GetParticleContainer(0)->SetFilterHybridTracks(kTRUE);
  pQATask->GetParticleContainer(0)->SetParticlePtCut(0.15);
  pQATask->SetAODfilterBits(256, 512);
  pQATask->SelectCollisionCandidates(kPhysSel);
  pQATask->SetHistoBins(200, 0, 30);
  pQATask->SetDefaultClusterEnergy(AliVCluster::kHadCorr);

  // Charged jet analysis
  if (bDoChargedJets) {
    AliEmcalJetTask *pChJetTask = AddTaskEmcalJet(sTracksName, "", 1, kJetRadius, 1, kTrackPtCut, kClusPtCut, kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, kFALSE);
    pChJetTask->SelectCollisionCandidates(kPhysSel);
    pChJetTask->SetFilterHybridTracks(kTRUE);
    sChJetsName = pChJetTask->GetName();

    AliAnalysisTaskSAJF *pSpectraChTask = AddTaskSAJF(sTracksName, "", sChJetsName, "",  kJetRadius, kJetPtCut, 0., "TPC");
    pSpectraChTask->SetNLeadingJets(1);
    pSpectraChTask->SelectCollisionCandidates(kPhysSel);
    pSpectraChTask->SetHistoType(kHistoType);
  }

  // Full jet analysis
  if (bDoFullJets) {
    AliEmcalJetTask *pFuJetTask = AddTaskEmcalJet(sTracksName, sClusName, 1, kJetRadius, 0, kTrackPtCut, kClusPtCut, kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, kFALSE);
    pFuJetTask->SelectCollisionCandidates(kPhysSel);
    pFuJetTask->SetFilterHybridTracks(kTRUE);
    pFuJetTask->SetClusterEnergyType(AliVCluster::kHadCorr);
    sFuJetsName = pFuJetTask->GetName();

    AliAnalysisTaskSAJF *pSpectraFuTask = AddTaskSAJF(sTracksName, sClusName, sFuJetsName, "", kJetRadius, kJetPtCut, 0., "EMCAL"); 
    pSpectraFuTask->SetNLeadingJets(1);
    pSpectraFuTask->SelectCollisionCandidates(kPhysSel);
    pSpectraFuTask->SetHistoType(kHistoType);
  }

  if (!pMgr->InitAnalysis()) return;
  pMgr->PrintStatus();

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
  pMgr->SetUseProgressBar(1, 250);
  //pMgr->SetDebugLevel(2);

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
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalTrackPropagator.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALTender.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskClusterizerFast.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C"); 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAJF.C");
}
