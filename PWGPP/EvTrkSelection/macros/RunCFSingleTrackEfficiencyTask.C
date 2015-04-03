class AliAnalysisGrid;
class AliAnalysisAlien;

//_______________________________| Loading Libraries |________________________________
void Load() {
    
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS/TPC -I$ALICE_PHYSICS/CONTAINERS -I$ALICE_PHYSICS/STEER/STEER -I$ALICE_PHYSICS/STEER/STEERBase -I$ALICE_PHYSICS/STEER/ESD -I$ALICE_PHYSICS/STEER/AOD -I$ALICE_PHYSICS/TRD -I$ALICE_PHYSICS/macros -I$ALICE_PHYSICS/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGPP -g"); 

  //load the required aliroot libraries
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libPhysics");
  gSystem->Load("libVMC");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGPP");
}

//_______________________________| Running on Grid |________________________________
void RunCFSingleTrackEfficiencyTask()
{

  const Bool_t readAOD = 1;

  TString        analysisMode    =  "local"; // "local", "grid", or "proof"
  TString           inputMode    =  "list"; // "list", "xml", or "dataset"
  Long64_t           nentries    =   123567890,firstentry=0;
  Bool_t       useAlienPlugin    =   kTRUE;
  TString          pluginmode    =   "test";
  TString testfileslistWithPlugin="filesAOD.txt"; // list of local files to test a la "files.txt" to use by the plugin

  TBenchmark benchmark;
  benchmark.Start("AliCFSingleTrackEfficiencyTask");

  Load();

  if(analysisMode=="grid") {
    TGrid::Connect("alien://") ;    //  Create an AliRunTagCuts and an AliEventTagCuts Object
  }

  if(useAlienPlugin) {
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,testfileslistWithPlugin,readAOD);
    if(!alienHandler) return;
  }

  printf("CREATE ANALYSIS MANAGER\n");
  AliAnalysisManager *mgr = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  if(useAlienPlugin) mgr->SetGridHandler(alienHandler);

  AliMCEventHandler*  mcHandler = new AliMCEventHandler();
  if (!readAOD) mgr->SetMCtruthEventHandler(mcHandler);
	 
  AliInputEventHandler* dataHandler;
  if   (readAOD) dataHandler = new AliAODInputHandler();
  else           dataHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(dataHandler);
	 
  //CREATE THE TASK
  printf("Prepare to create the task\n");

  // Run physics selection if not reading AODs
  if (!readAOD) {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kTRUE);
  }

  // Add new tasks
  //
  // First add the task for the PID response setting
  gROOT->LoadMacro("$ALICE_PHYSICS/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(kTRUE,kTRUE);
	 
  gROOT->LoadMacro("$ALICE_PHYSICS/ANALYSIS/macros/AddTaskPIDqa.C");
  AliAnalysisTaskPIDqa *pidQA = AddTaskPIDqa();

  // gROOT->LoadMacro("AliSingleTrackEffCuts.cxx++g");
  // gROOT->LoadMacro("AliCFSingleTrackEfficiencyTask.cxx++g");
  gROOT->LoadMacro("AddSingleTrackEfficiencyTask.C");
  AliCFSingleTrackEfficiencyTask *task = AddSingleTrackEfficiencyTask(readAOD,"Nch");
  AliCFSingleTrackEfficiencyTask *taskPi = AddSingleTrackEfficiencyTask(readAOD,"Pion",AliPID::kPion,211);
  AliCFSingleTrackEfficiencyTask *taskKa = AddSingleTrackEfficiencyTask(readAOD,"Kaon",AliPID::kKaon,321);
 
  // Run the analysis
  TChain * analysisChain=0;
	 
  if(analysisChain) printf("CHAIN HAS %d ENTRIES\n",(Int_t)analysisChain->GetEntries());
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
  if(analysisMode!="proof") {
    mgr->StartAnalysis(analysisMode.Data(),analysisChain,nentries,firstentry);
  }
	 
  benchmark.Stop("AliCFSingleTrackEfficiencyTask");
  benchmark.Show("AliCFSingleTrackEfficiencyTask");
	 
  return;  
}


//_______________________________| CreateAlienHandler |________________________________
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="test", TString testfileslistWithPlugin="", const Bool_t readAOD=kTRUE)
{

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetRunMode(pluginmode.Data());
  plugin->SetUser("zconesa");
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-08");
  plugin->SetAliROOTVersion("v5-05-03-AN");
  plugin->SetNtestFiles(1);

  // Set data file list to test on local mode  
  plugin->SetFileForTestMode(testfileslistWithPlugin.Data());

  // Set data search pattern for DATA grid Mode
  plugin->SetGridDataDir("/alice/sim/2013/LHC13d3"); // specify MC sample
  if(readAOD) plugin->SetDataPattern("AOD/*AliAOD.root");// specify AOD set
  else plugin->SetDataPattern("*/AliESDs.root");
  
  Int_t totruns=0;
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/AddGoodRuns.C");
  totruns += AddGoodRuns(plugin,"LHC13b","LHC13b");
  totruns += AddGoodRuns(plugin,"LHC13c","LHC13c");
  plugin->SetNrunsPerMaster(totruns);
  
  // plugin->AddDataFile("/alice/cern.ch/user/z/zconesa/sim/LHC13d3/195483_195529.xml");

  plugin->SetGridWorkingDir("sim/LHC13d3/ST290713");
  plugin->SetGridOutputDir("out"); 
  
  plugin->SetExecutable("ST290713.sh");
  plugin->SetAnalysisSource("AliSingleTrackEffCuts.cxx AliCFSingleTrackEfficiencyTask.cxx");
  plugin->SetAdditionalLibs("AliSingleTrackEffCuts.h AliSingleTrackEffCuts.cxx AliCFSingleTrackEfficiencyTask.cxx AliCFSingleTrackEfficiencyTask.h");  
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS/TPC -I$ALICE_PHYSICS/CONTAINERS -I$ALICE_PHYSICS/STEER/STEER -I$ALICE_PHYSICS/STEER/STEERBase -I$ALICE_PHYSICS/STEER/ESD -I$ALICE_PHYSICS/STEER/AOD -I$ALICE_PHYSICS/TRD -I$ALICE_PHYSICS/macros -I$ALICE_PHYSICS/ANALYSIS -I$ALICE_PHYSICS/OADB -g"); 
  
  plugin->SetDefaultOutputs(kTRUE);
  // merging via jdl
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(2);
  
  plugin->SetSplitMaxInputFileNumber(5);
  
  plugin->SetAnalysisMacro("ST290713.C");
  plugin->SetJDLName("TaskHFST290713.jdl");
  
  return plugin;
}
