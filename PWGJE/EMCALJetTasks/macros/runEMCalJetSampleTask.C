/// \file runEMCalJetSampleTask.C
/// \brief Example macro to run a EMCal-Jet-Framework-based analysis
///
/// \ingroup EMCALJETFW
/// This macros illustrates how to run a jet analysis. It can run locally or on
/// grid (test/full/terminate modes).
/// The script runEMCalJetSampleTask.sh in the same folder allows to easily set
/// the input values
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Apr 27, 2016


class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisManager;
class AliPhysicsSelectionTask;
class AliCentralitySelectionTask;
class AliEmcalCorrectionTask;
class AliEmcalJetTask;
class AliAnalysisGrid;

void LoadMacros();
void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode);
AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers,
    const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, Int_t workerTTL, Bool_t isMC);

//______________________________________________________________________________
AliAnalysisManager* runEMCalJetSampleTask(
    const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
    const char   *cRunPeriod     = "LHC11h",                                // set the run period
    const char   *cLocalFiles    = "fileLists/files_LHC11h_2_AOD145.txt",   // set the local list file
    const UInt_t  iNumEvents     = 5000,                                    // number of events to be analyzed
    const UInt_t  kPhysSel       = AliVEvent::kAnyINT |
    AliVEvent::kCentral | AliVEvent::kSemiCentral,                          // physics selection
    const char   *cTaskName      = "EMCalJetAna",                           // sets name of analysis manager
    const Bool_t  bDoChargedJets = kTRUE,
    const Bool_t  bDoFullJets    = kTRUE,
    const char   *obsolete       = "",                                      // Previous handled the ocdb settings, but obsolete due to CDBconnect task
    // 0 = only prepare the analysis manager but do not start the analysis
    // 1 = prepare the analysis manager and start the analysis
    // 2 = launch a grid analysis
    Int_t         iStartAnalysis = 2,
    const UInt_t  iNumFiles      = 100,                                     // number of files analyzed locally
    const char   *cGridMode      = "full"
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
      sRunPeriod == "lhc13d" || sRunPeriod == "lhc13e" || sRunPeriod == "lhc13f" ||
      sRunPeriod == "lhc16q" || sRunPeriod == "lhc16r" || sRunPeriod == "lhc16s" ||
      sRunPeriod == "lhc16t") {
    iBeamType = AliAnalysisTaskEmcal::kpA;
  }

  Double_t kGhostArea = 0.01;
  if (iBeamType != AliAnalysisTaskEmcal::kpp) kGhostArea = 0.005;

  AliTrackContainer::SetDefTrackCutsPeriod(sRunPeriod);

  Printf("Default track cut period set to: %s", AliTrackContainer::GetDefTrackCutsPeriod().Data());

  const Bool_t   bDoEmcalCorrections  = kTRUE;
  const Double_t kClusECut            = 0.30;
  const Double_t kTrackPtCut          = 0.15;
  const Double_t kJetPtCut            = 1.;

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
  if (iStartAnalysis == 1) {
    if (sLocalFiles == "") {
      Printf("You need to provide the list of local files!");
      return 0;
    }
    Printf("Setting local analysis for %d files from list %s, max events = %d", iNumFiles, sLocalFiles.Data(), iNumEvents);
  }

  LoadMacros();

  // Analysis manager
  AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);

  if (iDataType == kAod) {
    AliAODInputHandler* pAODHandler = AliAnalysisTaskEmcal::AddAODHandler();
  }
  else {  
    AliESDInputHandler* pESDHandler = AliAnalysisTaskEmcal::AddESDHandler();
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

  // CDBconnect task
  if (bDoFullJets || iDataType == kEsd) {
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
    taskCDB->SetFallBackToRaw(kTRUE);
  }

  if (bDoEmcalCorrections) {
    // Configuration of the Correction Task is handled via a YAML file, which is setup below
    // NOTE: Calling this function is equivalent to the AddTask macro, just the function is defined in the class.
    //       The AddTask macro still exists for use on the LEGO train
    AliEmcalCorrectionTask * correctionTask = AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask();
    correctionTask->SelectCollisionCandidates(kPhysSel);
    correctionTask->SetForceBeamType(iBeamType);

    // Configure and initialize
    correctionTask->SetUserConfigurationFilename("$ALICE_PHYSICS/PWG/EMCAL/config/PWGJESampleConfig.yaml");
    correctionTask->Initialize();
  }

  // Background
  TString sRhoChName;
  TString sRhoFuName;
  if (iBeamType != AliAnalysisTaskEmcal::kpp) {
    sRhoChName = "Rho";
    sRhoFuName = "Rho_Scaled";

    AliEmcalJetTask *pKtChJetTask = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "", AliJetContainer::kt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 0., kFALSE, kFALSE);
    pKtChJetTask->SelectCollisionCandidates(kPhysSel);

    AliAnalysisTaskRho* pRhoTask = AddTaskRhoNew("usedefault", "usedefault", sRhoChName, 0.4);
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
    AliEmcalJetTask *pChJet02Task = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "", AliJetContainer::antikt_algorithm, 0.2, AliJetContainer::kChargedJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 1., kFALSE, kFALSE);
    pChJet02Task->SelectCollisionCandidates(kPhysSel);

    AliEmcalJetTask *pChJet04Task = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 1., kFALSE, kFALSE);
    pChJet04Task->SelectCollisionCandidates(kPhysSel);
  }

  // Full jet analysis
  if (bDoFullJets) {
    AliEmcalJetTask *pFuJet02Task = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "usedefault", AliJetContainer::antikt_algorithm, 0.2, AliJetContainer::kFullJet, 0.15, 0.30, kGhostArea, AliJetContainer::pt_scheme, "Jet", 1., kFALSE, kFALSE);
    pFuJet02Task->SelectCollisionCandidates(kPhysSel);
    pFuJet02Task->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);

    AliEmcalJetTask *pFuJet04Task = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "usedefault", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kFullJet, 0.15, 0.30, kGhostArea, AliJetContainer::pt_scheme, "Jet", 1., kFALSE, kFALSE);
    pFuJet04Task->SelectCollisionCandidates(kPhysSel);
    pFuJet04Task->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  // Sample task
  AliAnalysisTaskEmcalJetSample *sampleTask = 0;
  sampleTask = AddTaskEmcalJetSample("usedefault", "usedefault", "usedefault", "new");
  sampleTask->GetClusterContainer(0)->SetClusECut(0.);
  sampleTask->GetClusterContainer(0)->SetClusPtCut(0.);
  sampleTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
  sampleTask->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
  sampleTask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  sampleTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
  sampleTask->SetHistoBins(600, 0, 300);
  sampleTask->SelectCollisionCandidates(kPhysSel);

  if (bDoFullJets) {
    AliJetContainer* jetCont02 = sampleTask->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.2, AliEmcalJet::kEMCALfid, "Jet");
    AliJetContainer* jetCont04 = sampleTask->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.4, AliEmcalJet::kEMCALfid, "Jet");

    if (iBeamType != AliAnalysisTaskEmcal::kpp) {
      jetCont02->SetRhoName(sRhoFuName);
      jetCont02->SetPercAreaCut(0.6);
      jetCont04->SetRhoName(sRhoFuName);
      jetCont04->SetPercAreaCut(0.6);
    }
  }

  if (bDoChargedJets) {
    AliJetContainer* jetCont02 = sampleTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.2, AliEmcalJet::kTPCfid, "Jet");
    AliJetContainer* jetCont04 = sampleTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.4, AliEmcalJet::kTPCfid, "Jet");
    if (iBeamType != AliAnalysisTaskEmcal::kpp) {
      jetCont02->SetRhoName(sRhoChName);
      jetCont02->SetPercAreaCut(0.6);
      jetCont04->SetRhoName(sRhoChName);
      jetCont04->SetPercAreaCut(0.6);
    }
  }

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
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetSample.C");
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
  plugin->SetAliPhysicsVersion("vAN-20170628-1");

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
