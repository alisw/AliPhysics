const TString kTaskName = "Filter";
const int     kJobSplitting = 20;
const TString kGridWorkDir = "PbPb";
const TString kGridOutputDir = "output";
const Int_t source = 1;

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const TString gridmode = "test", bool isMC = kFALSE){
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(gridmode.Data());
  // set number of test files to use in "test" mode
  plugin->SetNtestFiles(2);
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20170221-1");
  plugin->SetExecutableCommand("aliroot -b -q");
 // Declare input data to be processed.
  TArrayI runlist;
  Int_t start = 0, nrun = 0;
  //Setting the run ranges
  start = 0;
  nrun = 12;
  //Setting the run list
  if (isMC){
  }
  else {
    plugin->SetGridDataDir("/alice/data/2015/LHC15o/");
    plugin->SetDataPattern("*/pass3_lowIR_pidfix/*AliESDs.root");
    plugin->SetRunPrefix("000");
    Int_t rl[12] = {
      244918, 244975, 244980, 244982, 244983, 245061, 245064, 245066, 245068, 246390, 246391, 246392
  	};
    runlist.Set(12,rl);
  }
  Info("CreateAlienHandler", "Setting %s run list from %s %s", isMC ? "MC" : "Data", plugin->GetGridDataDir(), plugin->GetDataPattern());
  for(Int_t i = start; i < start + nrun; i++) {
    if(i >= runlist.GetSize()) Fatal("CreateAlienHandler", "Adding runs which are out of bound");
    plugin->AddRunNumber(runlist[i]);
  }
  plugin->SetNrunsPerMaster(1);
  plugin->SetOutputToRunNo();
  plugin->SetGridWorkingDir(Form("%s_%s",kGridWorkDir.Data(),isMC ? "MC" : ""));
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir(kGridOutputDir.Data());
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  // plugin->SetAdditionalLibs("libPWGLFnuclex.so");
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  //plugin->SetAdditionalLibs("AliAnalysisCODEXtask.cxx AliAnalysisCODEXtask.h AliAnalysisCODEX.cxx AliAnalysisCODEX.h");
  //plugin->SetAnalysisSource("AliAnalysisCODEX.cxx AliAnalysisCODEXtask.cxx");
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
  // Modify merging options
  plugin->SetMergeViaJDL(kFALSE);
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  return plugin;
}

void runFilter(TString gridmode = "test", bool isMC = false) {
  // check run type
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  if(!gridmode.Contains("local") && !gridmode.Contains("test") && !gridmode.Contains("grid")) {
    Printf("\n\tIncorrect run option, check first argument of run macro");
    Printf("\tint gridmode = local, test or grid\n");
    return;
  }
  Printf("%s analysis chosen",gridmode.Data());
  AliAnalysisGrid *plugin = 0x0;
  TChain *chain = 0x0;
  if (gridmode != "local") {
    plugin = CreateAlienHandler(gridmode,isMC);
    if(!plugin) return;
  } else {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
    chain = CreateAODChain("localdata.txt");
  }
  //---- Create the analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(kTaskName);
  if(plugin) mgr->SetGridHandler(plugin);
  //  Input
  if (isMC) {
    AliMCEventHandler*  mcHandler = new AliMCEventHandler();
    if (plugin) mgr->SetMCtruthEventHandler(mcHandler);
  }
  AliESDInputHandler* iH = new AliESDInputHandler("handler","handler for my analisys");
  mgr->SetInputEventHandler(iH);
  // Other tasks
  // Physics selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSel = AddTaskPhysicsSelection(isMC); // useMC
  // Centrality selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask * multSelTask = AddTaskMultSelection(kFALSE); // user mode:
  // PID response
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(isMC); // useMC
  // PID QA
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  AliAnalysisTaskPIDqa *pidQATask = AddTaskPIDqa();
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/NUCLEX/Utils/CODEX/AddCODEXtask.C");
  AliAnalysisCODEXtask * filter = AddCODEXtask();
  //Setting custom option for each data period
  filter->mCentralityMode = 1;
  filter->Discard("");
  // Disable debug printouts
  mgr->SetDebugLevel(3);
  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  AliLog::SetGlobalDebugLevel(0);
  if (!mgr->InitAnalysis()) return;

  mgr->PrintStatus();
  // Start analysis in grid.
  if (gridmode.Contains("local"))
    mgr->StartAnalysis(gridmode,chain);
  else
    mgr->StartAnalysis("grid");
}
