// $Id$

void runJetAna(
  const char     *datatype     = "aod",            // aod, esd, sesd
  const char     *runtype      = "local",          // local or grid (when local gridmode specifies input txt file)
  const char     *gridmode     = "aod_files.txt",  // grid mode (can be "full", "test", "offline", "submit" or "terminate")
  const char     *taskname     = "JetAna")         // name of grid generated macros
{

  enum eDataType { kAod, kEsd, kSesd };
  enum eRunType  { kLocal, kGrid };

  eRunType rType = kLocal;
  if (strcmp(runtype, "grid")==0) 
    rType = kGrid;
  eDataType dType = kAod;
  if (strcmp(datatype, "esd")==0)
    dType = kEsd;
  else if (strcmp(datatype, "sesd")==0)
    dType = kSesd;

  // load the libraries
  LoadLibs();
   
  // analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(taskname);

  if (dType == kAod) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
    AliAODInputHandler* inH = AddAODHandler();
  } else {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
    AliESDInputHandler* inH = AddESDHandler();
  }

  if (0) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODOutputHandler.C");
    AliAODHandler* aodoutHandler = AddAODOutputHandler();
  }

  // PSel task
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE, 
                                                                      AliVEvent::kAnyINT /*| AliVEvent::kCentral| AliVEvent::kSemiCentral*/,
                                                                      10,10);

  // Setup task
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
  AliEmcalSetupTask *setupTask = AddTaskEmcalSetup();

  // Compatibility task (for skimmed ESD)
  if (dType == kSesd) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalCompat.C");
    AliEmcalCompatTask *comptask = AddTaskEmcalCompat();
  }

  // Centrality task
  if (dType == kEsd) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *centralityTask = AddTaskCentrality();
  }

  TString inputTracks("tracks");
  if (dType == kEsd) {
    inputTracks = "HybridTracks";

    // Hybrid tracks maker for ESD
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTpcTrack.C");
    AliEmcalEsdTpcTrackTask *hybTask = AddTaskEmcalEsdTpcTrack(inputTracks);

    // Track propagator
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalTrackPropagator.C");
    AliEmcalTrackPropagatorTask *propTask = AddTaskEmcalTrackPropagator(inputTracks);
  }
  else if (dType == kSesd) {
    inputTracks = "Tracks";
  }

  // PicoTracks maker
  TString tracksName("PicoTracks");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(tracksName, inputTracks, "LHC11h");

  TString clusName("CaloClusters");
  if (dType == kAod)
    clusName = "caloClusters";

  // EmcalParticle maker
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  AliEmcalParticleMaker *emcalParts = 
    AddTaskEmcalParticleMaker(tracksName,clusName,"EmcalTracks","EmcalClusters");

  // Cluster-track matcher
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  AliEmcalClusTrackMatcherTask *matcherTask = AddTaskEmcalClusTrackMatcher("EmcalTracks","EmcalClusters");

  // Hadronic correction task
  TString clusNameCorr(Form("%sCorr",clusName.Data()));
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskHadCorr.C");
  AliHadCorrTask *hcorr = AddTaskHadCorr("EmcalTracks","EmcalClusters",clusNameCorr);

  if (0) {
  // Embedding task
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbedding.C");
  AliJetEmbeddingTask* jemb = AddTaskJetEmbedding(tracksName, clusNameCorr, "JetEmbeddingTask", 10, 10, -0.9, 0.9);

  // Jet finder
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  AliEmcalJetTask *jetTask = AddTaskEmcalJet(tracksName, clusNameCorr);

  // Scale task
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskScale.C");
  AliAnalysisTaskScale *scaleTask = AddTaskScale(tracksName, clusNameCorr);
  }
  if (1) {
    UInt_t val = AliVEvent::kAny;
    //val = AliVEvent::kAnyINT | AliVEvent::kCentral| AliVEvent::kSemiCentral;
    //val = AliVEvent::kEMCEGA;
    //val = AliVEvent::kEMCEJE;
    val = AliEmcalPhysicsSelection::kEmcalHT;

    TObjArray *toptasks = mgr->GetTasks();
    for (Int_t i=0; i<toptasks->GetEntries(); ++i) {
      AliAnalysisTaskSE *task = dynamic_cast<AliAnalysisTaskSE*>(toptasks->At(i));
      if (!task)
        continue;
      TString name(task->ClassName());
      if (name.Contains("PhysicsSelection"))
        continue;
      ::Info("setPSel", "Set physics selection for %s (%s)", task->GetName(), task->ClassName());
      task->SelectCollisionCandidates(val);
    }
  }

  mgr->SetDebugLevel(0);
  mgr->SetUseProgressBar(1, 25);

  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  
  if (rType == kGrid) {
    mgr->StartAnalysis(gridmode);
  } else {
    const char *txtfile = gridmode;
    if (dType == kAod) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/CreateAODChain.C");
      chain = CreateAODChain(txtfile, 5);
    } else {
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/CreateESDChain.C");
      TChain* chain = CreateESDChain(txtfile, 5);
    }
    mgr->StartAnalysis("local", chain);
  }

  return;
}

void LoadLibs()
{
  // load root libraries
  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libGui");
  gSystem->Load("libXMLParser");
  gSystem->Load("libMinuit");
  gSystem->Load("libMinuit2");
  gSystem->Load("libProof");
  gSystem->Load("libGeom");
  gSystem->Load("libPhysics");

  // load aliroot libraries
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libSTEER");
  gSystem->Load("libCORRFW");
  gSystem->Load("libEMCALUtils");
  gSystem->Load("libPHOSUtils");
  gSystem->Load("libPWGCaloTrackCorrBase");
  gSystem->Load("libPWGGACaloTrackCorrelations");
  gSystem->Load("libPWGGACaloTasks");
  gSystem->Load("libEMCALraw");
  gSystem->Load("libEMCALbase");
  gSystem->Load("libEMCALrec");
  gSystem->Load("libTRDbase");
  gSystem->Load("libVZERObase");
  gSystem->Load("libVZEROrec");
  gSystem->Load("libTENDER");   
  gSystem->Load("libTENDERSupplies"); 
 
  // load fastjet libraries
  gSystem->Load("libJETAN");
  gSystem->Load("/usr/lib64/libCGAL.so");
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libSISConePlugin");
  gSystem->Load("libFASTJETAN");
  gSystem->Load("libPWGEMCAL");
  gSystem->Load("libPWGGAEMCALTasks");
  gSystem->Load("libPWGJEEMCALJetTasks");
}

AliAnalysisGrid* CreateAlienHandler(const char *taskname, 
                                    const char *gridmode, 
                                    const char *proofcluster, 
                                    const char *proofdataset)
{
  // TODO
  return 0;
}

