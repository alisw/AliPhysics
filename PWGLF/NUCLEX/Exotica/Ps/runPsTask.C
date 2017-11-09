const TString kTaskName = "PsEfficiency";
const int     kJobSplitting = 20;
const TString kGridWorkDir = "pPb";
const TString kGridOutputDir = "output";
const Int_t source = 0;

void runPsTask(const char* mode = "test", bool isMC=true){
  // Include AliPhysics and AliRoot libraries
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  //Check the run mode
  if(std::strcmp(mode,"local")!=0 && std::strcmp(mode,"grid")!=0 && std::strcmp(mode,"test")!=0){
    Printf("\n\tIncorrect run option, check first argument of run macro");
    Printf("\tmode = local, test or grid\n");
    return;
  }
  Printf("%s analysis chosen",mode);
  //Define the plugin
  AliAnalysisGrid *plugin = 0x0;
  TChain *chain = 0x0;
  if(std::strcmp(mode,"local")!=0) {
    plugin = CreateAlienHandler(mode,isMC);
    if(!plugin) return;
  }
  else {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
    chain = CreateAODChain("localdata.txt");
  }
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager(kTaskName);
  if(plugin) mgr->SetGridHandler(plugin);
  // Define input handlers
  if(isMC){
    AliMCEventHandler*  mcHandler = new AliMCEventHandler();
    if (plugin) mgr->SetMCtruthEventHandler(mcHandler);
  }
  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  // Compile the class
  gROOT->LoadMacro("AliAnalysisTaskPsEfficiency.cxx++g");
  // PID response
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(isMC); // useMC
  // Ps efficiency task
  gROOT->LoadMacro("AddTaskPsEfficiency.C");
  AliAnalysisTaskPsEfficiency *task = AddTaskPsEfficiency(isMC);
  // Disbale debug printouts
  mgr->SetDebugLevel(3);
  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  AliLog::SetGlobalDebugLevel(0);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  // Start analysis in grid.
  if (std::strcmp(mode,"local")==0)
    mgr->StartAnalysis(mode,chain);
  else
    mgr->StartAnalysis("grid");
}

AliAnalysisGrid* CreateAlienHandler(const char* mode, bool isMC){
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(gridmode.Data());
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20170425-1");
  plugin->SetExecutableCommand("aliroot -b -q");
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  //plugin->SetAdditionalLibs("libPWGLFnuclex.so");
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  // Declare all libraries (other than the default ones for the framework. These will benloaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.

  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  // To only save certain files, use SetDefaultOutputs(kFALSE), and then
  // SetOutputFiles("list.root other.filename") to choose which files to save
  plugin->SetDefaultOutputs();
  //plugin->SetOutputFiles("list.root");

  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro(Form("%s.C",kTaskName.Data()));

  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(kJobSplitting);

  // Optionally modify the executable name (default analysis.sh)
  plugin->SetExecutable(Form("%s.sh",kTaskName.Data()));

  // set number of test files to use in "test" mode
  plugin->SetNtestFiles(2);

  // file containing a list of chuncks to be used for testing
  plugin->SetFileForTestMode("testdata");

  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);

  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);

  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");

  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName(Form("%s.jdl",kTaskName.Data()));

  // Optionally modify job price (default 1)
  plugin->SetPrice(1);

  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");

  return plugin;
}
