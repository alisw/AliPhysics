/// **** to run the macro: ***********
// alien-token-init username
// root -l runGrid.C

class AliAnalysisGrid;

void runGridCascadePbPb( Bool_t   useMC               = kTRUE,  // kTRUE if analysing a MC sample 
                         Bool_t   runperformancetask  = kTRUE,  
                         Float_t  centrlowlim         = 0.,
                         Float_t  centruplim          = 90.,
                         TString  centrest            = "V0M",
                         Float_t  vtxlim              = 10.,
                         Int_t minnTPCcls             = 70, 
                         Bool_t   kextrasel           = kFALSE,
                         Bool_t   kacccut              = kFALSE,
                         Bool_t   krelaunchvertexers  = kFALSE,
                         TString  anatype             = "AOD",//"ESD",
                         TString  gridoutputdir       = "LHC10h_AOD086",
                         //the following are used for the Cascade task only
                         TString  datadir             = "/alice/data/2010/LHC10h",///alice/data/2011/LHC11h_2",
                         TString  datapattern         = "ESDs/pass2/AOD086/*/AliAOD.root", // "ESDs/pass2/*/*ESDs.root" // for data
                         const char *plugin_mode      ="full") {

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
  AliAnalysisGrid *alienHandler = CreateAlienHandler(plugin_mode, runperformancetask, useMC, anatype, gridoutputdir, datadir, datapattern);
  if (!alienHandler) return;
 
  //__________________________________________________________________________
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("CascadePbPbanalysis");
  
  //__________________________________________________________________________
  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);

  // Input handlers
  AliESDInputHandler* esdH = new AliESDInputHandler();
  AliAODInputHandler* aodH = new AliAODInputHandler();
  if (anatype=="ESD") mgr->SetInputEventHandler(esdH);
  else mgr->SetInputEventHandler(aodH);
  if (runperformancetask&&(anatype=="ESD")) {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
  }

  //__________________________________________________________________________
  // Add tasks

  if (anatype=="ESD") {
    // Physics selection
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *physSel = AddTaskPhysicsSelection(useMC);
    // Centrality selection
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentr = AddTaskCentrality();
    if (useMC){
      taskCentr->SetMCInput();
      taskCentr->DontUseCleaning();
    }
  }

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(useMC);

  if (runperformancetask) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/AliAnalysisTaskCheckPerformanceCascadePbPb.cxx++g");
    gROOT->LoadMacro("$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/macros/AddTaskCheckPerformanceCascadePbPb.C");
    AliAnalysisTaskCheckPerformanceCascadePbPb *task = AddTaskCheckPerformanceCascadePbPb(minnTPCcls, centrlowlim, centruplim, centrest,vtxlim,kextrasel ,kacccut ,krelaunchvertexers);

  } else {
    gROOT->LoadMacro("$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/AliAnalysisTaskCheckCascadePbPb.cxx++g");
    gROOT->LoadMacro("$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/macros/AddTaskCheckCascadePbPb.C");
    AliAnalysisTaskCheckCascadePbPb *task = AddTaskCheckCascadePbPb(minnTPCcls, centrlowlim, centruplim, centrest,vtxlim,kextrasel ,krelaunchvertexers);

  }
 


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

AliAnalysisGrid* CreateAlienHandler(const char *plugin_mode, Bool_t runperformancetask, Bool_t useMC, TString anatype, TString gridoutputdir, TString datadir, TString datapattern) {
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
  plugin->SetROOTVersion("v5-30-06-1");
  plugin->SetAliROOTVersion("v5-03-03-AN"); 

  //__________________________________________________________________________
  // Declare input data to be processed.
  // Method 1: Create automatically XML collections using alien 'find' command.
  if (useMC) {
    //plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis");   // Define production directory
    //plugin->SetGridDataDir("/alice/sim/LHC11a10b_plus");
    plugin->SetGridDataDir("/alice/sim/2011/LHC11f5");
    // Set data search pattern
    if (anatype == "ESD") plugin->SetDataPattern("*ESDs.root");
    else plugin->SetDataPattern("AOD090/*AOD.root");  
    plugin->AddRunNumber(137124);


  } else {
    plugin->SetGridDataDir(datadir.Data());   // Define production directory LFN
    plugin->SetDataPattern(datapattern.Data());  // Set data search pattern
    plugin->SetRunPrefix("000");
    //plugin->SetRunRange(80000,80000); // ...then add run numbers to be considered
    plugin->AddRunNumber(138534);
   
  }
  // Method 2: Use your input collection 
  //plugin->AddDataFile("/alice/cern.ch/user/m/mnicassi/139105.xml");
  //__________________________________________________________________________
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir(gridoutputdir.Data());
  plugin->SetGridOutputDir("output");

  //__________________________________________________________________________
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
    plugin->SetAdditionalLibs("AliAnalysisTaskCheckCascadePbPb.h AliAnalysisTaskCheckCascadePbPb.cxx ");
    // Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable("CascadePbPb.sh");

  }
  
  //__________________________________________________________________________
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  //plugin->SetDefaultOutputs(kFALSE);
  //if (runperformancetask) plugin->SetOutputFiles("CascadePerformance.root");
  //else plugin->SetOutputFiles("Cascades.root");

  // Optionally define the files to be archived.
  //plugin->SetOutputArchive("root_archive.zip:*.root log_archive.zip:stdout,stderr");

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

  // Optionally modify job price (default 1)
  plugin->SetPrice(1);

  // Merge via JDL
  // comment out the next line when using the "terminate" option, unless
  // you want separate merged files for each run
/*  plugin->SetMergeViaJDL(kTRUE);  // run first in full mode, then in terminate
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeFiles(50);
  plugin->SetMaxMergeStages(2);// to define the number of stages
*/
  // Optionally set number of runs per master
//  plugin->SetNrunsPerMaster(1);
  //
  plugin->SetOutputToRunNo();
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  plugin->SetUser("mnicassi");
  return plugin;
}

