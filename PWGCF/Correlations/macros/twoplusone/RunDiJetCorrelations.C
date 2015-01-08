
/* $Id: AliAnalysisTaskDiJetCorrelations.cxx */
class AliAnalysisGrid;
class AliAnalysisAlien;

void RunDiJetCorrelations()
{
  
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGCF -I$ALICE_ROOT/PWG/FLOW/Tasks -g");
  
  //TString mySE=" ALICE::CERN::SE";
  //gSystem->Setenv("alien_CLOSE_SE",mySE.Data());
  //Use AliRoot includes to compile our task
  //gROOT->ProcessLine(".include $ALICE_ROOT/include");
  //gROOT->ProcessLine(".include $PWD/.");
    
  // Seeting for analysis run
  TString      analysisMode = "grid"; // "local", "grid", or "proof"
  TString        pluginmode = "terminate"; // full terminate or test mode
  Bool_t     useAlienPlugin =  kTRUE;
  TString         inputMode = "list"; // "list", "xml", or "dataset"
  Long64_t         nentries =  123567890,  firstentry=0;
  TString testfileslistWithPlugin = "/Users/varmaraghava/Greeshma_Analysis/MECheck/filesAOD.txt";
    
  // Setting the AddTaskOptions but NOT Used for the moment
  Bool_t             readMC =  kFALSE; // DATA or MC data set ON<->OFF
  Bool_t              genMC =  kFALSE; // MC Generation ON<->OFF
  Bool_t             tracks =  kTRUE; // DATA or MC Reco ON<->OFF
  Bool_t             mixing =  kFALSE; // Event Mixing ON<->OFF
  Bool_t           trackEff =  kTRUE; //  track eff ON<->OFF
  Bool_t               Cent =  kFALSE; //system kTRUE means pA/PbPb
  
  
  if(analysisMode=="grid") {
    TGrid::Connect("alien://");
  } 
  
  // AliRoot libraries
  if(analysisMode=="local" || analysisMode=="grid") {
    TString loadLibraries="./LoadLibraries.C"; //loadLibraries.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(loadLibraries.Data());
    LoadLibraries();
  }
  
   
  // Create Alien plugin, if requested 
  // Prepare input
  TChain *chainAOD = 0;  
  if(useAlienPlugin) {  
    //if(analysisMode!="grid") {printf("Analysis mode must be grid, to use alien plugin\n"); return;}
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,testfileslistWithPlugin.Data());
    if(!alienHandler) return;
  }  else{
    TString makeAODInputChain="../MakeAODInputChain.C"; makeAODInputChain.Prepend(loadMacroPath.Data());
    if(inputMode=="list") {
      // Local files
      gROOT->LoadMacro(makeAODInputChain.Data());
      chainAOD = MakeAODInputChain();// with this it reads ./AliAOD.root and ./AliAOD.VertexingHF.root
      printf("ENTRIES %d\n",chainAOD->GetEntries());
    } else if(inputMode=="xml") {
      // xml
      gROOT->LoadMacro(makeAODInputChain.Data());
      chainAOD = MakeAODInputChain("collection_aod.xml","collection_aodHF.xml");
    } 
  }
  
    
  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  if(useAlienPlugin) mgr->SetGridHandler(alienHandler);
  
  // Input
  AliAODInputHandler *inputHandler = new AliAODInputHandler("handler","handler for Di-Jets Correlatiosn");
  mgr->SetInputEventHandler(inputHandler);
  

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  gROOT->LoadMacro("AliAnalysisTaskDiJetCorrelations.cxx++g");   
  gROOT->LoadMacro("AddTaskDiJetCorrelations.C");
  AliAnalysisTaskDiJetCorrelations *taskSE1  = AddTaskDiJetCorrelations("SE", kFALSE, kTRUE, 12.0, 16.0, 5.0, 8.0);
  AliAnalysisTaskDiJetCorrelations *taskME1  = AddTaskDiJetCorrelations("ME", kTRUE,  kTRUE, 12.0, 16.0, 5.0, 8.0);

  //AliAnalysisTaskDiJetCorrelations *taskSE3  = AddTaskDiJetCorrelations("SE_3", kFALSE, kTRUE, 12.0, 30.0, 5.0, 10.0);
  //AliAnalysisTaskDiJetCorrelations *taskME4  = AddTaskDiJetCorrelations("ME_4", kTRUE,  kTRUE, 12.0, 30.0, 5.0, 10.0);

  // Run the analysis   
  if(chainAOD) printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainAOD->GetEntries());
  
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
    
  
  if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
 
  if(analysisMode!="proof") {
    mgr->StartAnalysis(analysisMode.Data(),chainAOD,nentries,firstentry);
  } else {
    mgr->StartAnalysis(analysisMode.Data(),dataset.Data(),nentries,firstentry);
  }
  
  return;
  
}

//___________________________________| Alien Handler function and plugin|______________________________________
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="test", TString testfileslistWithPlugin="")
{
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(pluginmode.Data());
  plugin->SetUser("rvarma");
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-08-6");
  plugin->SetAliROOTVersion("vAN-20140702");
  plugin->SetFileForTestMode(testfileslistWithPlugin.Data());
  plugin->SetNtestFiles(1); // if "test mode"
  
  
  //________________________Set data search pattern for DATA and MC
  //Method 1: To create automatically xml through plugin
  gROOT->LoadMacro("./AddGoodRuns.C");
  plugin->SetGridDataDir("/alice/data/2010/LHC10h/"); // specify LHC period
  plugin->SetDataPattern("*/AOD086/*/AliAOD.root"); // specify reco pass and AOD set
  Int_t totruns=0;
  totruns += AddGoodRuns(plugin,"LHC10h_Gir",""); // specify LHC period
  plugin->SetNrunsPerMaster(totruns);    
  
  plugin->SetGridWorkingDir("TEST/Dcorr/Dec8/");
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  plugin->SetExecutable("AnalysisDiJetCorr.sh");
  plugin->SetAnalysisSource("AliAnalysisTaskDiJetCorrelations.cxx");

  // Declare all addtional libraries  here.  
  plugin->SetAdditionalLibs("AliAnalysisTaskDiJetCorrelations.h AliAnalysisTaskDiJetCorrelations.cxx libCORRFW.so libEMCALUtils.so libJETAN.so THnEfficiency.root"); 
  plugin->AddIncludePath("-I. -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/CORRFW -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/JETAN -I$ALICE_ROOT/PWG4");

  plugin->SetDefaultOutputs(kTRUE);
  plugin->SetSplitMaxInputFileNumber(10);
  
  // merging via jdl
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(2);
  plugin->SetMaxMergeFiles(10);
  plugin->SetTTL(70000);
  plugin->SetKeepLogs();
  plugin->SetAnalysisMacro("AnalysisDiJetCorr.C");
  plugin->SetJDLName("DiJetCorr.jdl");
  
  return plugin;

}
