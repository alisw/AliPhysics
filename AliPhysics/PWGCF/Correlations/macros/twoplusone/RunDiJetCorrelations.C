

class AliAnalysisGrid;
class AliAnalysisAlien;


//______________| Loading Proper Libraries
void Load() {
    
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGCF -I$ALICE_ROOT/PWGCF/Correlations -I$ALICE_ROOT/PWGCF/Correlations/DPhi -I$ALICE_ROOT/PWGCF/Correlations/macros -g");
    
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    
    //load the required aliroot libraries
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libMinuit.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libCORRFW.so");
    gSystem->Load("libPWGTools.so");
    gSystem->Load("libPWGCFCorrelationsBase.so");
    gSystem->Load("libPWGCFCorrelationsDPhi.so");


}

//______________| Run Di-Jet Function
void RunDiJetCorrelations()
{
    
    
  // Seeting for analysis run
  Load();
  TString      analysisMode = "grid"; // "local", "grid", or "proof"
  TString        pluginmode = "full"; // full terminate or test mode
  Bool_t     useAlienPlugin =  kTRUE;
  TString         inputMode = "list"; // "list", "xml", or "dataset"
  Long64_t         nentries =  123567890,  firstentry=0;
  TString testfileslistWithPlugin = "./filesAOD.txt";
  TString loadMacroPath="$ALICE_ROOT/PWGCF/Correlations/macros/";
  
  if(analysisMode=="grid"){
    TGrid::Connect("alien://");
  } 
   
  // Create Alien plugin, if requested 
  TChain *chainAOD = 0;
  if(useAlienPlugin){
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
  AliAODInputHandler *inputHandler = new AliAODInputHandler("handler","handler for Di-Jets Correlations");
  mgr->SetInputEventHandler(inputHandler);
  
  gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
    
  gROOT->LoadMacro("AliAnalysisTaskDiJetCorrelations.cxx++g");
  gROOT->LoadMacro("AddTaskDiJetCorrelations.C");
    AliAnalysisTaskDiJetCorrelations *taskSE1  = AddTaskDiJetCorrelations("SE_1", kFALSE, kTRUE, 7.0, 12.0, 3.0, 8.0);
    AliAnalysisTaskDiJetCorrelations *taskME1  = AddTaskDiJetCorrelations("ME_1", kTRUE, kTRUE,  7.0, 12.0, 3.0, 8.0);
  
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
  plugin->SetUser("gkoyitha");
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-08-6");
  plugin->SetAliROOTVersion("vAN-20141230");
    //plugin->SetAliROOTVersion("vAN-20141021");
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
  
  plugin->SetGridWorkingDir("/Jan13/CrossCheck/");
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  plugin->SetExecutable("AnalysisDiJetCorr.sh");
  plugin->SetAnalysisSource("AliAnalysisTaskDiJetCorrelations.cxx");
                            

  // Declare all addtional libraries  here.  
  plugin->SetAdditionalLibs("libCORRFW.so libEMCALUtils.so libJETAN.so AliAnalysisTaskDiJetCorrelations.h AliAnalysisTaskDiJetCorrelations.cxx");
   
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGCF -I$ALICE_ROOT/PWGCF/Correlations -I$ALICE_ROOT/PWGCF/Correlations/DPhi -I$ALICE_ROOT/PWGCF/Correlations/macros -g");

  plugin->SetDefaultOutputs(kTRUE);
  plugin->SetSplitMaxInputFileNumber(10);
  
  // merging via jdl
  //plugin->SetMergeViaJDL(kFALSE);
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(2);
  plugin->SetMaxMergeFiles(10);
  plugin->SetTTL(30000);
  plugin->SetKeepLogs();
  plugin->SetAnalysisMacro("AnalysisDiJetCorr.C");
  plugin->SetJDLName("DiJetCorr.jdl");
  
  return plugin;

}




// Loading Libs



