const TString kTaskName = "ampt_validation";
const int     kJobSplitting = 100;
const TString kGridWorkDir = "AMPT_string_melting_off";
const TString kGridOutputDir = "output_ycut1";
const bool    kUseXml = true;

// note on this macro: to use the task you need at least AliPhysics::vAN-20160906-1

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const TString gridmode = "test", bool isMC = kFALSE)
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetRunMode(gridmode.Data());

  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20160906-1");
  plugin->SetExecutableCommand("aliroot -b -q");

  if (!kUseXml) plugin->SetGridDataDir("/alice/sim/2016/LHC16f4b/");
  for (int i = 1; i <= 12; ++i)
    if (!kUseXml) plugin->AddRunNumber(((i - 1) / 6) * 10 + (i - 1) % 6 + 1);
    else plugin->AddDataFile(Form("/alice/cern.ch/user/m/mpuccio/%s/%i.xml",kGridWorkDir.Data(),((i - 1) / 6) * 10 + (i - 1) % 6 + 1));
  plugin->SetUseMCchain();

  plugin->SetNrunsPerMaster(1);
  plugin->SetOutputToRunNo();

  plugin->SetGridWorkingDir(kGridWorkDir.Data());

  plugin->SetGridOutputDir(kGridOutputDir.Data());

  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  plugin->SetAdditionalLibs("libpythia6.so libAliPhytia6.so");
  //plugin->SetAnalysisSource("AliAnalysisTaskNucleiKine.cxx");

  plugin->SetDefaultOutputs();

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

void runESD (TString gridmode = "test", bool isMC = true) {
  // check run type
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  if(!gridmode.Contains("test") && !gridmode.Contains("full")) {
    Printf("\n\tIncorrect run option, check first argument of run macro");
    Printf("\tint gridmode = full or test\n");
    return;
  }
  Printf("%s analysis chosen",gridmode.Data());

  AliAnalysisGrid *plugin = 0x0;
  TChain *chain = 0x0;
  plugin = CreateAlienHandler(gridmode,isMC);

  //---- Create the analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(kTaskName);
  if(plugin) mgr->SetGridHandler(plugin);
  AliESDEvent *esdE = new AliESDEvent();
  esdE->CreateStdContent();
  AliESDVertex *vtx = new AliESDVertex(0.,0.,100);
  vtx->SetName("VertexTracks");
  vtx->SetTitle("VertexTracks");
  esdE->SetPrimaryVertexTracks(vtx);
  AliDummyHandler *dumH = new AliDummyHandler;//  static_cast<AliDummyHandler*>mgr->GetInputEventHandler();
  dumH->SetEvent(esdE);
  mgr->SetInputEventHandler(dumH);

  AliMCEventHandler*  mcHandler = new AliMCEventHandler();
  mcHandler->SetReadTR(false);
  if (plugin) mgr->SetMCtruthEventHandler(mcHandler);

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/NUCLEX/Nuclei/NucleiKine/AddTaskNucleiKine.C");
  AliAnalysisTaskNucleiKine * filter = AddTaskNucleiKine();//kTRUE);

  // Disable debug printouts
  mgr->SetDebugLevel(3);
  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  AliLog::SetGlobalDebugLevel(0);
  if (!mgr->InitAnalysis()) return;

  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");
}
