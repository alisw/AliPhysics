// $Id$

void runJetAna(const char     *datatype     = "aod",        // aod, esd, sesd
               const char     *runtype      = "local",      // local or grid
               const char     *gridmode     = "test",       // run mode (can be "full", "test", "offline", "submit" or "terminate")
               const char     *txtfile      = "ifiles.txt", // text file with input files for local mode
               const char     *taskname     = "JetAna")     // name of grid generated macros
{

  enum eDataType { kAod, kEsd, kSesd };
  enum eRunType  { kLocal, kGrid };

  eRunType rType = kLocal;
  if (!strcmp(runtype, "grid")) 
    rType = kGrid;
  eDataType dType = kAod;
  if (!strcmp(runtype, "esd"))
    dType = kEsd;
  else if (!strcmp(runtype, "sesd"))
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
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEmcalPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskEmcalPhysicsSelelection(kTRUE);

  // Centrality task
  if (dType == kEsd) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *centralityTask = AddTaskCentrality();
  }

  // Setup task
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEmcalSetup.C");
  AliEmcalSetupTask *setupTask = AddTaskEmcalSetup();

  // Track maker
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEmcalPicoTrackMaker.C");
  AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker("PicoTracks", "tracks", "LHC11h");

  // Cluster-track matcher
  if (0) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEmcalClusTrackMatcher.C");
    AliEmcalClusTrackMatcherTask *matcherTask = AddTaskEmcalClusTrackMatcher("PicoTracks", "caloClusters");
  }

  // Hadronic correction task
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALJetTasks/macros/AddTaskHadCorr.C");
  AliHadCorrTask *hcorr = AddTaskHadCorr("PicoTracks", "caloClusters", "caloClustersCorr");
  
  // Jet finder
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  AliEmcalJetTask *jetTask = AddTaskEmcalJet("PicoTracks", "caloClustersCorr");

  if (1) {
    UInt_t val = AliVEvent::kAny;
    val = AliVEvent::kAnyINT | AliVEvent::kCentral| AliVEvent::kSemiCentral;
    //val = AliVEvent::kEMCEGA;
    //val = AliVEvent::kEMCEJE;

    TObjArray *toptasks = mgr->GetTasks();
    for (Int_t i=0; i<toptasks->GetEntries(); ++i) {
      AliAnalysisTaskSE *task = dynamic_cast<AliAnalysisTaskSE*>(toptasks->At(i));
      if (!task)
        continue;
      TString name(task->GetName());
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
    if (dType == kAod) {
      gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALJetTasks/macros/CreateAODChain.C");
      chain = CreateAODChain("files_aod95.txt", 50);
    } else {
      gROOT->LoadMacro("$ALICE_ROOT/PWGUD/macros/CreateESDChain.C");
      TChain* chain = CreateESDChain(ifiles.txt, 50);
    }
    mgr->StartAnalysis("local", chain);
  }

  return;
}

/*
  //Tracks maker
  //gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEmcalAodTrackFilter.C");
  //AliEmcalAodTrackFilterTask *eTask = AddTaskEmcalAodTrackFilter(tracksName, "tracks", "LHC11h");

  //Tender Supplies
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEmcalAodTender.C");
  //geometry EMCAL_COMPLETEV1 or EMCAL_FIRSTYEARV1; data pp or PbPb
  AliEmcalTenderTask *tender = AddTaskEmcalAodTender("EMCAL_COMPLETEV1", "PbPb");
  /*
  if (runtype == "grid") {
    tender->SetDefaultCDBStorage("raw://"); //uncomment if you work on grid
  }
  else if (runtype == "local") {
    tender->SetDefaultCDBStorage("local://$ALICE_ROOT/OCDB"); //uncomment if you work local
  }
  */
  /*
  //V1unfold Clusterizer
  TString AODbranchName;
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEMCALClusterize.C");
  AliAnalysisTaskEMCALClusterize *v1UnfoldClusTask = AddTaskEMCALClusterize(AODbranchName);
  
  //L0-L1 Clusterizer
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskClusterizerFW.C");
  //AliAnalysisTaskEMCALClusterizeFast *L0ClusTask = AddTaskClusterizerFW("L0");
  AliAnalysisTaskEMCALClusterizeFast *L0ClusTask = AddTaskClusterizerFW("L1GAMMA");
  //Hadronic correction task
  //gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALJetTasks/macros/AddTaskHadCorr.C");
  //AliHadCorrTask *hcorr = AddTaskHadCorr(tracksName, clustersName, corrClusName);
  
  //Cluster Track matcher
  //gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEmcalClusTrackMatcher.C");
  //AliEmcalClusTrackMatcherTask *matcherTask = AddTaskEmcalClusTrackMatcher(tracksName, corrClusName, 1, 1);

  //Jet finder
  //gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  //AliEmcalJetTask *jetTask = AddTaskEmcalJet(tracksName, corrClusName);
	
  // create task
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEmcalIsolatedPhotons.C");
  //AliEmcalIsolatedPhotonsTask *task = AddTaskEmcalIsolatedPhotons(tracksName, corrClusName, jetsName);
  //task->SelectCollisionCandidates(AliVEvent::kAnyINT);  // Any MB trigger
  //task->SelectCollisionCandidates(AliVEvent::kEMCEGA);  // Gamma trigger
  //task->SelectCollisionCandidates(AliVEvent::kEMCEJE);  // Jet trigger
  */

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
  gSystem->Load("libPWGGAEMCALTasks");
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
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libSISConePlugin");
  gSystem->Load("libFASTJETAN");
  gSystem->Load("libPWGGAEMCALJetTasks");
}

AliAnalysisGrid* CreateAlienHandler(const char *taskname, 
                                    const char *gridmode, 
                                    const char *proofcluster, 
                                    const char *proofdataset)
{
  // TODO
  return 0;
}

