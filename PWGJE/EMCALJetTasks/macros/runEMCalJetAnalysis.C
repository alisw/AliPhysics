// runEMCalJetAnalysis.C
// =====================
// This macro can be used to run a jet analysis within the EMCal Jet Framework.
//
// Examples:
// -> Analyze ESDs from the pA pilot run on the AliEn grid with your task in AnaClass.cxx/.h
//     dataType = "ESD", useGrid = kTRUE, pattern = "*ESDs/pass2/*ESDs.root", addCXXs = "AnaClass.cxx", 
//     addHs = "AnaClass.h", gridDir = "/alice/data/2012/LHC12g", gridMode = "full", runNumbers = "188359 188362"
//     
// -> Analyze AODs (up to 96 files) locally given in files_aod.txt
//     dataType = "AOD", useGrid = kFALSE, numLocalFiles = 96
//
// MERGING ON ALIEN
// ++++++++++++++++
// If you run on the grid, you can monitor the jobs with alimonitor.cern.ch. When enough of them are in DONE state,
// you have to merge the output. This can be done automatically, if you just change the gridMode to "terminate" and
// give the EXACT name of the task whose output should be merged in uniqueName.
// 
//
// Authors: R. Haake, S. Aiola

#include <ctime>
#include "TGrid.h"

AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers, 
                                     const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, 
                                     Int_t workerTTL, Bool_t isMC);
                                    
//______________________________________________________________________________
void runEMCalJetAnalysis(
         Bool_t         useGrid             = kTRUE,                      // local or grid
         const char*    gridMode            = "test",                      // set the grid run mode (can be "full", "test", "offline", "submit" or "terminate")
         const char*    dataType            = "AOD",                       // set the analysis type, AOD, ESD or sESD
         const char*    pattern             = "*ESDs/pass2/AOD145/*AOD.root",    // file pattern (here one can specify subdirs like passX etc.) (used on grid)
         const char*    gridDir             = "/alice/data/2011/LHC11h_2",   // dir on alien, where the files live (used on grid)
         const char*    runNumbers          = "167903 167915",             // considered run numbers (used on grid)
         UInt_t         numLocalFiles       = 5,                           // number of files analyzed locally
         const char*    runPeriod           = "LHC11h",                    // set the run period (used on grid)
         const char*    uniqueName          = "EMCalJF_LEGOTrainTest",     // sets base string for the name of the train on the grid
         UInt_t         pSel                = AliVEvent::kAny,             // used event selection for every task
         Bool_t         useTender           = kTRUE,                       // trigger, if tender, track and cluster selection should be used (always)
         Bool_t         isMC                = kFALSE,                      // trigger, if MC handler should be used
         Bool_t         doBkg               = kTRUE,
         // Here you have to specify additional code files you want to use but that are not in aliroot
         const char*    addCXXs             = "",
         const char*    addHs               = "",

         // These two settings depend on the dataset and your quotas on the AliEN services
         Int_t          maxFilesPerWorker   = 4,
         Int_t          workerTTL           = 7200

         )
{

  // Some pre-settings and constants
  enum AlgoType {kKT, kANTIKT};
  enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};
  gSystem->SetFPEMask();
  gSystem->Setenv("ETRAIN_ROOT", ".");
  gSystem->Setenv("ETRAIN_PERIOD", runPeriod);
  // change this objects to strings
  TString usedData(dataType);
  TString additionalCXXs(addCXXs);
  TString additionalHs(addHs);
  cout << dataType << " analysis chosen" << endl;
  if (useGrid)  
  {
    cout << "-- using AliEn grid.\n";
    if (usedData == "sESD") 
    {
      cout << "Skimmed ESD analysis not available on the grid!" << endl;
      return;
    }
  }
  else
    cout << "-- using local analysis.\n";
  

  // Load necessary libraries
  LoadLibs();

  // Create analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(uniqueName);

  // Check type of input and create handler for it
  TString localFiles("-1");
  if(usedData == "AOD")
  {
    localFiles = "files_aod.txt";
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
    AliAODInputHandler* aodH = AddAODHandler();
  }
  else if((usedData == "ESD") || (usedData == "sESD"))
  {
    if (usedData == "ESD")
      localFiles = "files_esd.txt";
    else
      localFiles = "files_sesd.txt";
    
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
    AliESDInputHandler* esdH = AddESDHandler();
  }
  else
  {
    cout << "Data type not recognized! You have to specify ESD, AOD, or sESD!\n";
  }

  if(!useGrid)
    cout << "Using " << localFiles.Data() << " as input file list.\n";

  // Create MC handler, if MC is demanded
  if (isMC && (usedData != "AOD"))
  {
    AliMCEventHandler* mcH = new AliMCEventHandler();
    mcH->SetPreReadMode(AliMCEventHandler::kLmPreRead);
    mcH->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(mcH); 
  }
  
  // ################# Now: Add some basic tasks

  // Physics selection task
  if(!isMC) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
    AliPhysicsSelectionTask *physSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE, pSel, 5, 5, 10, kTRUE, -1, -1, -1, -1);
    if (!physSelTask) {
      cout << "no physSelTask but running on data" << endl; 
      return; 
    }
  }

  // Centrality task
  if (usedData == "ESD") {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *centralityTask = AddTaskCentrality(kTRUE);
  }

  // Compatibility task, only needed for skimmed ESD
  if (usedData == "sESD") {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalCompat.C");
    AliEmcalCompatTask *comptask = AddTaskEmcalCompat();
  }

  // Setup task
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
  AliEmcalSetupTask *setupTask = AddTaskEmcalSetup();
  setupTask->SetGeoPath("$ALICE_PHYSICS/OADB/EMCAL");
  
  // Tender Supplies
  if (useTender) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPreparation.C");
    //adjust pass when running locally. On grid give empty string, will be picked up automatically from path to ESD/AOD file
    AliAnalysisTaskSE *clusm = AddTaskEmcalPreparation(runPeriod);
  }

  // Names of the different objects passed around; these are the default names; added here mostly for documentation purposes
  // rhoName is only set if the background calculation is switched on (doBkg)
  TString tracksName = "PicoTracks";
  TString clustersName = "EmcCaloClusters";
  TString clustersCorrName = "CaloClustersCorr";
  TString rhoName = "";

  // ################# Now: Call jet preparation macro (picotracks, hadronic corrected caloclusters, ...) 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetPreparation.C");
  TString particlesMCName = "";
  if(isMC) particlesMCName = "MCParticlesSelected";
  AddTaskJetPreparation(runPeriod, tracksName, particlesMCName.Data(), clustersName, clustersCorrName);

  // ################# Now: Add jet finders+analyzers
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  AliEmcalJetTask* jetFinderTask = AddTaskEmcalJet(tracksName, clustersCorrName, kANTIKT, 0.2, kCHARGEDJETS, 0.150, 0.300);

  if (doBkg) {
    rhoName = "Rho";
    AliEmcalJetTask* jetFinderTaskKT = AddTaskEmcalJet(tracksName, clustersCorrName, kKT, 0.2, kCHARGEDJETS, 0.150, 0.300);

    TString kTpcKtJetsName = jetFinderTaskKT->GetName();
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRho.C");
    rhotask = (AliAnalysisTaskRho*) AddTaskRho(kTpcKtJetsName, tracksName, clustersCorrName, rhoName, 0.2, "TPC", 0.01, 0, 0, 2, kTRUE);
    //rhotask__->SetScaleFunction(sfunc);
    //rhotask->SelectCollisionCandidates(kPhysSel);
    rhotask->SetHistoBins(100,0,250);
  }
  // Here you can put in your AddTaskMacro for your task
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetSample.C");
  AliAnalysisTaskEmcalJetSample* anaTask = 0;
  AddTaskEmcalJetSample(tracksName, clustersCorrName, jetFinderTask->GetName(), rhoName, 4);

  // Set the physics selection for all given tasks
  TObjArray *toptasks = mgr->GetTasks();
  for (Int_t i=0; i<toptasks->GetEntries(); ++i) 
  {
    AliAnalysisTaskSE *task = dynamic_cast<AliAnalysisTaskSE*>(toptasks->At(i));
    if (!task)
      continue;
    if (task->InheritsFrom("AliPhysicsSelectionTask"))
      continue;
    ::Info("setPSel", "Set physics selection for %s (%s)", task->GetName(), task->ClassName());
    task->SelectCollisionCandidates(pSel);
  }

  mgr->SetUseProgressBar(1, 25);
        
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();

  if (useGrid) 
  {  // GRID CALCULATION

    AliAnalysisGrid *plugin = CreateAlienHandler(uniqueName, gridDir, gridMode, runNumbers, pattern, additionalCXXs, additionalHs, maxFilesPerWorker, workerTTL, isMC);
    mgr->SetGridHandler(plugin);

    // start analysis
    cout << "Starting GRID Analysis...";
    mgr->SetDebugLevel(0);
    mgr->StartAnalysis("grid");
  }
  else
  {  // LOCAL CALCULATION

    TChain* chain = 0;
    if (usedData == "AOD") 
    {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
      chain = CreateAODChain(localFiles.Data(), numLocalFiles);
    }
    else
    {  // ESD or skimmed ESD
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      chain = CreateESDChain(localFiles.Data(), numLocalFiles);
    }
    
    // start analysis
    cout << "Starting LOCAL Analysis...";
    mgr->SetDebugLevel(2);
    mgr->StartAnalysis("local", chain);
  }
}

//______________________________________________________________________________
void LoadLibs()
{
  // Load common libraries (better too many than too few)

  //load CGAL, Fastjet and SISCone
  if(gSystem->Load("/usr/lib/libCGAL") != 0 ) gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  gSystem->Load("libfastjettools");
  gSystem->Load("libfastjetcontribfragile");
  //
  gSystem->Load("libJETAN");
  gSystem->Load("libPWGJEEMCALJetTasks");

}

AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers, 
                                     const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, 
                                     Int_t workerTTL, Bool_t isMC)
{
  TDatime currentTime;
  TString tmpName(uniqueName);

  // Only add current date and time when not in terminate mode! In this case the exact name has to be supplied by the user
  if(strcmp(gridMode, "terminate"))
  {
    tmpName += "_";
    tmpName += currentTime.GetDate();
    tmpName += "_";
    tmpName += currentTime.GetTime();
  }

  TString tmpAdditionalLibs("");
  tmpAdditionalLibs = Form("libCGAL.so libJETAN.so libfastjet.so libsiscone.so libsiscone_spherical.so libfastjetplugins.so libfastjettools.so libfastjetcontribfragile.so libPWGJE.so libPWGmuon.so libPWGJEEMCALJetTasks.so %s %s",additionalCode.Data(),additionalHeaders.Data());


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
  plugin->SetAliPhysicsVersion("vAN-20151101");


  plugin->SetGridDataDir(gridDir); // e.g. "/alice/sim/LHC10a6"
  plugin->SetDataPattern(pattern); //dir structure in run directory
  if (!isMC)
   plugin->SetRunPrefix("000");

  plugin->AddRunList(runNumbers);

  plugin->SetGridWorkingDir(Form("work/%s",tmpName.Data()));
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output

  plugin->SetAnalysisSource(additionalCode.Data());
  plugin->SetAdditionalLibs(tmpAdditionalLibs.Data());
//  plugin->AddExternalPackage("boost::v1_43_0");
//  plugin->AddExternalPackage("cgal::v3.6");
//  plugin->AddExternalPackage("fastjet::v2.4.2");

  plugin->SetDefaultOutputs(kTRUE);
  //plugin->SetMergeExcludes("");
  plugin->SetAnalysisMacro(macroName.Data());
  plugin->SetSplitMaxInputFileNumber(maxFilesPerWorker);
  plugin->SetExecutable(execName.Data());
  plugin->SetTTL(workerTTL);
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName(jdlName.Data());
  plugin->SetPrice(1);      
  plugin->SetSplitMode("se");

  return plugin;
}
