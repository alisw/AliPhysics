// **** to run the macro: ***********
// alien-token-init username
// source /tmp/gclient_env_IDuser
// root -l runGrid.C

class AliAnalysisGrid;

void runGridCascadePbPb ( Bool_t   useMC             = kTRUE,  // kTRUE if analysing a MC sample 
                          Bool_t   runperformancetask = kTRUE,  
                          Float_t  centrlowlim       = 0.,
                          Float_t  centruplim        = 80.,
                          TString  centrest         = "V0M",
                          Short_t  lCollidingSystems = 1,       //0 = pp, 1 = AA
                          Float_t  vtxlim             = 15.,
                          Bool_t   usecfcontainers     = kTRUE,
                          Bool_t   kextrasel           = kFALSE,
                          Bool_t   acccut = kTRUE,
                          const char *plugin_mode="test") {

  // Load common libraries
  gSystem->Load("libCore.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so"); 
  gSystem->Load("libGui.so");
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libCDB.so");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libProof.so");
  gSystem->Load("libRAWDatabase.so");
  gSystem->Load("libSTEER.so");

   //__________________________________________________________________________
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  
  //__________________________________________________________________________
  // Create and configure the alien handler plugin
//  gROOT->LoadMacro("CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler(plugin_mode, runperformancetask, useMC);
  if (!alienHandler) return;
 
  //__________________________________________________________________________
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("CascadePbPbanalysis");
  
  //__________________________________________________________________________
  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);

  //__________________________________________________________________________
  // Add task 

  if (runperformancetask) {
    gROOT->LoadMacro("AliAnalysisTaskCheckPerformanceCascadePbPb.cxx++g");
    AliAnalysisTaskCheckPerformanceCascadePbPb *task = new AliAnalysisTaskCheckPerformanceCascadePbPb("TaskPerformanceCascade");
    task->SetDebugLevelCascade           (0);                 // not used in the task
    task->SetApplyAccCut                 (acccut);
    task->SetUseCFCont                   (usecfcontainers);
    task->SetAlephParamFor1PadTPCCluster (kTRUE);             // to set aleph param - ask which ones have to be used
  } else {
    gROOT->LoadMacro("AliAnalysisTaskCheckCascadePbPb.cxx++g");
    AliAnalysisTaskCheckCascadePbPb *task = new AliAnalysisTaskCheckCascadePbPb("TaskCascade");
    task->SetUseCFContCascadeCuts       (usecfcontainers);
  }

 
  task->SetCollidingSystems           (lCollidingSystems); // only for multiplicity binning
  task->SetAnalysisType               ("ESD");
  task->SetRelaunchV0CascVertexers    (0);                 // used but code is commented out
  task->SetQualityCutZprimVtxPos      (kTRUE);             // selects vertices in +-10cm
  task->SetQualityCutNoTPConlyPrimVtx (kTRUE);             // retains only events with tracking + SPD vertex
  task->SetQualityCutTPCrefit         (kTRUE);             // requires TPC refit flag to be true to select a track
  task->SetQualityCut80TPCcls         (kTRUE);             // rejects tracks that have less than 80 clusters in the TPC
  task->SetExtraSelections            (kextrasel);         // used to add other selection cuts
  task->SetCentralityLowLim           (centrlowlim);       // setting centrality selection vriables
  task->SetCentralityUpLim            (centruplim);
  task->SetCentralityEst              (centrest);
  task->SetVertexRange                (vtxlim);
 
  mgr->AddTask(task);

  // ESD case
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);
  if (runperformancetask) {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
  }

  // Physics selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSel = AddTaskPhysicsSelection(useMC);
  task->SelectCollisionCandidates();
  // Centrality selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentr = AddTaskCentrality();
  if (useMC) taskCentr->SetMCInput();
 
  //Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  if (runperformancetask) AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist", TList::Class(), AliAnalysisManager::kOutputContainer, "CascadePerformance.root");
  else AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist", TList::Class(), AliAnalysisManager::kOutputContainer, "Cascades.root");
  //__________________________________________________________________________
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  
  //__________________________________________________________________________
  // Disbale debug printouts
  mgr->SetDebugLevel(3);
  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  AliLog::SetGlobalDebugLevel(0);
  
  //__________________________________________________________________________
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");
};

//__________________________________________________________________________

AliAnalysisGrid* CreateAlienHandler(const char *plugin_mode, Bool_t runperformancetask, Bool_t useMC) {
  //__________________________________________________________________________
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
  AliAnalysisAlien *plugin= new AliAnalysisAlien();

  //__________________________________________________________________________
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(plugin_mode);

  plugin->SetNtestFiles(1);
  //__________________________________________________________________________
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");

  //__________________________________________________________________________
  // On GRID - current
  plugin->SetROOTVersion("v5-27-06d");
  plugin->SetAliROOTVersion("v4-21-16-AN");

  //__________________________________________________________________________
  // Declare input data to be processed.
  // Method 1: Create automatically XML collections using alien 'find' command.
  if (useMC) {
    plugin->SetGridDataDir("/alice/sim/LHC11a7");   // Define production directory 
    // Set data search pattern
    plugin->SetDataPattern("*ESDs.root"); 
    plugin->SetRunRange(137161,138225); 
  } else {
    plugin->SetGridDataDir("/alice/data/2010/LHC10h");   // Define production directory LFN
    plugin->SetDataPattern("ESDs/pass2_rev15/*/*ESDs.root");  // Set data search pattern
    // plugin->SetRunRange(80000,80000); // ...then add run numbers to be considered
    plugin->SetRunPrefix("000");
//    plugin->AddRunNumber(137366);
    plugin->AddRunNumber(138200);
//    plugin->AddRunNumber(139172);
  }
  // Method 2: Use your input collection
   
  //plugin->AddDataFile("/alice/cern.ch/user/d/dcolella/wn.xml");

  //__________________________________________________________________________
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  if (runperformancetask) plugin->SetGridWorkingDir("workperf");
  else plugin->SetGridWorkingDir("work");
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output

  if (runperformancetask) plugin->SetAnalysisSource("AliAnalysisTaskCheckPerformanceCascadePbPb.cxx");
  else plugin->SetAnalysisSource("AliAnalysisTaskCheckCascadePbPb.cxx");

  //__________________________________________________________________________
  //Enable same others packages
  //plugin->EnablePackage("PWG3dielectron.par");

  //__________________________________________________________________________
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("testmacro.C");
  //Add all extra files (task.cxx/.h)
  if (runperformancetask) {
    plugin->SetAdditionalLibs("AliAnalysisTaskCheckPerformanceCascadePbPb.h AliAnalysisTaskCheckPerformanceCascadePbPb.cxx");
    // Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable("CascadePerformancePbPb.sh");

  } else {
    plugin->SetAdditionalLibs("AliAnalysisTaskCheckCascadePbPb.h AliAnalysisTaskCheckCascadePbPb.cxx");
    // Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable("CascadePbPb.sh");

  }
  //__________________________________________________________________________
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  plugin->SetDefaultOutputs(kFALSE);
  if (runperformancetask) plugin->SetOutputFiles("CascadePerformance.root");
  else plugin->SetOutputFiles("Cascades.root");

  // Optionally define the files to be archived.
  plugin->SetOutputArchive("root_archive.zip:*.root log_archive.zip:stdout,stderr");

  //__________________________________________________________________________
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(100);
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  if (runperformancetask) plugin->SetJDLName("TaskCheckPerformanceCascadePbPb.jdl");
  else plugin->SetJDLName("TaskCheckCascadePbPb.jdl");

  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(100);

  // Optionally modify job price (default 1)
  plugin->SetPrice(1);

  // Merge via JDL
  // comment out the next line when using the "terminate" option, unless
  // you want separate merged files for each run
//  plugin->SetMergeViaJDL(kTRUE);  // run first in full mode, then in terminate
//  plugin->SetOneStageMerging(kFALSE);
//  plugin->SetMaxMergeStages(2);// to define the number of stages
  // Optionally set number of runs per master
  plugin->SetNrunsPerMaster(1);
  //
  plugin->SetOutputToRunNo();
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  plugin->SetUser("mnicassi");
  return plugin;
}

