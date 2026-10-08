class AliAnalysisGrid;
class AliAnalysisAlien;

void RunAnalysisAODVertexingHF()
{
  //
  // Test macro for AliAnalysisTaskSE's for heavy-flavour candidates
  // It has the structure of a Analysis Train:
  // - in this macro, change things related to running mode
  //   and input preparation 
  // - add your task using a AddTaskXXX macro 
  //
  // A.Dainese, andrea.dainese@lnl.infn.it
  // "grid" mode added by R.Bala, bala@to.infn.it
  //


  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/PWG -g"); 
  //
  TString trainName = "D2H";
  TString analysisMode = "grid"; // "local", "grid"
  TString inputMode    = "list"; // "list", "xml", or "dataset"
  Long64_t nentries=123567890,firstentry=0;
  Bool_t useParFiles=kFALSE;
  Bool_t useAlienPlugin=kTRUE;
  TString pluginmode="full";
  TString testfileslistWithPlugin="";
  TString loadMacroPath="$ALICE_PHYSICS/PWGHF/vertexingHF/macros/";
  //TString loadMacroPath="./"; // this is normally needed for CAF
  //

  if(analysisMode=="grid") {
    // Connect to AliEn
    TGrid::Connect("alien://");
  }


  // AliRoot libraries
  if(analysisMode=="local" || analysisMode=="grid") {
    TString loadLibraries="LoadLibraries.C"; loadLibraries.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(loadLibraries.Data());
    LoadLibraries(useParFiles);
    gSystem->Load("libGui");
    gSystem->Load("libRAWDatabase");
    gSystem->Load("libCDB");
    gSystem->Load("libSTEER");
    gSystem->Load("libTRDbase");
    gSystem->Load("libPWGTRD");
  }


  // Create Alien plugin, if requested
  if(useAlienPlugin) {  
    //    if(analysisMode!="grid") {printf("Analysis mode must be grid, to use alien plugin\n"); return;}
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,useParFiles,testfileslistWithPlugin);  
    if(!alienHandler) return;
  }


  //-------------------------------------------------------------------
  // Prepare input
  TChain *chainAOD = 0;
  TString dataset; // for proof

  if(!useAlienPlugin) {
    TString makeAODInputChain="../MakeAODInputChain.C"; makeAODInputChain.Prepend(loadMacroPath.Data());
    if(inputMode=="list") {
      // Local files
      gROOT->LoadMacro(makeAODInputChain.Data());
      chainAOD = MakeAODInputChain();// with this it reads ./AliAOD.root and ./AliAOD.VertexingHF.root
      //chainAOD = MakeAODInputChain("alien:///alice/cern.ch/user/r/rbala/newtrain/out_lhc08x/180100/",1,1);
      printf("ENTRIES %d\n",chainAOD->GetEntries());
    } else if(inputMode=="xml") {
      // xml
      gROOT->LoadMacro(makeAODInputChain.Data());
      chainAOD = MakeAODInputChain("collection_aod.xml","collection_aodHF.xml");
    } else if(inputMode=="dataset") {
      // CAF dataset
      dataset="/ITS/dainesea/AODVertexingHF_LHC08x_180100";
    }
  }

  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  // Connect plug-in to the analysis manager
  if(useAlienPlugin) mgr->SetGridHandler(alienHandler);

  // Input
  AliAODInputHandler *inputHandler = new AliAODInputHandler("handler","handler for D2H");

  mgr->SetInputEventHandler(inputHandler);
  //-------------------------------------------------------------------

  
  //-------------------------------------------------------------------
  // Analysis tasks (wagons of the train)   
  //
  // First add the task for the PID response setting
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(kFALSE,kTRUE);

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  AliAnalysisTaskPIDqa *pidQA = AddTaskPIDqa();


  TString taskName;
  
  ////// ADD THE FULL D2H TRAIN
  /*taskName="../AddD2HTrain.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  Bool_t readMC=kFALSE;
  AddD2HTrain(readMC);//,1,0,0,0,0,0,0,0,0,0,0);*/
  
  ////// OR ADD INDIVIDUAL TASKS
  
  taskName="AddTaskHFQA.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAnalysisTaskSEHFQA *qaTask = AddTaskHFQA(0,"",kFALSE, kTRUE, 0 , "QApp",kTRUE,kTRUE,kFALSE, kTRUE, kFALSE,kTRUE);

  /*  taskName="AddTaskCompareHF.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSECompareHF *cmpTask = AddTaskCompareHF();
    */  
    taskName="AddTaskD0Mass.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSED0Mass *d0massTask = AddTaskD0Mass();
    AliAnalysisTaskSED0Mass *d0massLikeSignTask = AddTaskD0Mass(1); 
    
  //-------------------------------------------------------------------
  
  //
  // Run the analysis
  //    
  if(chainAOD) printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainAOD->GetEntries());
  
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
  {

    mgr->StartAnalysis(analysisMode.Data(),chainAOD,nentries,firstentry);

}
  
  return;
}
//_____________________________________________________________________________
//
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="test",Bool_t useParFiles=kFALSE, TString testfileslistWithPlugin="")
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(pluginmode.Data());
   plugin->SetUser();
   // Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion();
   plugin->SetAliROOTVersion();
   plugin->SetNtestFiles(1);
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/AddGoodRuns.C");

   // Declare input data to be processed.
   //************************************************
   // Set data file list to test on local mode
   //************************************************  
   plugin->SetFileForTestMode(testfileslistWithPlugin.Data());

   //************************************************
   // Set data search pattern for DATA
   //************************************************  
   //Method 1: To create automatically xml through plugin  
   plugin->SetGridDataDir("/alice/data/2010/LHC10d"); // specify LHC period
   plugin->SetDataPattern("pass2/AOD018/*AliAOD.root"); // specify reco pass and AOD set
   plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
   // OR plugin->SetFriendChainName("deltas/AliAOD.VertexingHF.root");
   // Adds only the good runs from the Monalisa Run Condition Table
   // More than one period can be added but the period name has to be removed from GridDataDir (to be tested)
   Int_t totruns=0;
   //totruns += AddGoodRuns(plugin,"LHC10b"); // specify LHC period
   //totruns += AddGoodRuns(plugin,"LHC10c"); // specify LHC period
   totruns += AddGoodRuns(plugin,"LHC10d"); // specify LHC period
   plugin->SetNrunsPerMaster(totruns);

   // Method 2: Declare existing data files (e.g xml collections)

   //plugin->AddDataFile("/alice/cern.ch/user/r/rbala/000168068_000170593.xml");
   //  plugin->SetDataPattern("*AliAOD.root");
   //  plugin->SetFriendChainName("./AliAOD.VertexingHF.root"); 

   //************************************************
   // Set data search pattern for MONTECARLO
   //************************************************
   /* 
   plugin->SetGridDataDir("/alice/sim/LHC10d3"); // specify MC sample
   plugin->SetDataPattern("AOD005/*AliAOD.root"); // specify AOD set
   plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
   // OR plugin->SetFriendChainName("deltas/AliAOD.VertexingHF.root");
   // Adds only the good runs from the Monalisa Run Condition Table 
   // More than one period can be added!
   Int_t totruns=0;
   totruns += AddGoodRuns(plugin,"LHC10b","LHC10d3"); // specify LHC period for anchor runs; and the name of the MC production
   //totruns += AddGoodRuns(plugin,"LHC10c","LHC10f7"); // specify LHC period for anchor runs;  and the name of the MC production
   //totruns += AddGoodRuns(plugin,"LHC10d","LHC10f7"); // specify LHC period for anchor runs;  and the name of the MC production
   plugin->SetNrunsPerMaster(totruns);
   */
   //
   // Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("myHFanalysis");
   // Name of executable
   plugin->SetExecutable("myHFanalysis.sh");
   // Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
   // Declare the analysis source files names separated by blancs. To be compiled runtime
   // using ACLiC on the worker nodes.
   //plugin->SetAnalysisSource("AliDStarJets.cxx");
   // Declare all libraries (other than the default ones for the framework. These will be
   // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->SetAdditionalLibs("libPWGflowBase.so libPWGflowTasks.so libPWGHFbase.so libPWGHFvertexingHF.so libGui.so libRAWDatabase.so libCDB.so libSTEER.so libTRDbase.so libPWGTRD.so");
   // use par files
   if(useParFiles) {
     plugin->EnablePackage("STEERBase.par");
     plugin->EnablePackage("ESD.par");
     plugin->EnablePackage("AOD.par");
     plugin->EnablePackage("ANALYSIS.par");
     plugin->EnablePackage("OADB.par");
     plugin->EnablePackage("ANALYSISalice.par");
     plugin->EnablePackage("CORRFW.par");
     plugin->EnablePackage("PWGHFbase.par");
     plugin->EnablePackage("PWGHFvertexingHF.par");
   }
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/PWG -g"); 

   plugin->SetDefaultOutputs(kTRUE);
   // merging via jdl
   plugin->SetMergeViaJDL(kTRUE);
   plugin->SetOneStageMerging(kFALSE);
   plugin->SetMaxMergeStages(2);

   // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("AnalysisHF.C");
   // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   // Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("TaskHF.jdl");

   return plugin;
}
