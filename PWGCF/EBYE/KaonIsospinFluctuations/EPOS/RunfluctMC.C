class AliAnalysisGrid;
class AliAnalysisAlien;

void RunfluctMC()
{
 

gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I\"/usr/local/CERN/root/include\" -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OCDB -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGCF -I$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/AliFemto -I$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/AliFemtoUser -g");
  
  //
  TString trainName = "D2H";
  TString analysisMode = "grid"; // "local", "grid", or "proof"
  TString inputMode    = "list"; // "list", "xml", or "dataset"
  Long64_t nentries=12356789,firstentry=0;
  //Long64_t nentries=1000,firstentry=0;
  Bool_t useParFiles=kFALSE;
  Bool_t useAlienPlugin=kTRUE;
  TString pluginmode="test"; //test, full , terminate
  Bool_t saveProofToAlien=kFALSE;
  TString proofOutdir = "";
  TString loadMacroPath="/afs/cern.ch/work/s/sbaidyan/nudyn_EPOS/EPOSGen/";
  

  if(analysisMode=="grid") {
    // Connect to AliEn
    TGrid::Connect("alien://");
  } 


  // AliRoot libraries
  if(analysisMode=="local" || analysisMode=="grid") {
    TString loadLibraries="LoadLibraries.C"; loadLibraries.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(loadLibraries.Data());
    //  LoadLibraries(useParFiles);
  } 

  gROOT->ProcessLine(".include $ALICE_ROOT/include");


  // Create Alien plugin, if requested
  if(useAlienPlugin) {  
    if(analysisMode!="grid") {printf("Analysis mode must be grid, to use alien plugin\n"); return;}
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,useParFiles);  
    if(!alienHandler) return;
  }


  //-------------------------------------------------------------------
  // Prepare input
  //TChain *chainAOD = 0;
  TString dataset; // for proof

 

  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  // Connect plug-in to the analysis manager
  if(useAlienPlugin) mgr->SetGridHandler(alienHandler);

  // Input
   AliESDInputHandler *inputHandler = new AliESDInputHandler("handler","handler for D2H");
  //For AMPT
  AliMCEventHandler* mchandler = new AliMCEventHandler();
  mchandler->SetReadTR(kFALSE); 
mgr->SetMCtruthEventHandler(mchandler); // set MC handler for AliAnalysisManager
mgr->SetInputEventHandler(inputHandler);
  //-------------------------------------------------------------------

  
  //-------------------------------------------------------------------
  // Analysis tasks (wagons of the train)   
  //
  // First add the task for the PID response setting
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
//    AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(kFALSE,kTRUE);
    AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(kTRUE,kTRUE);
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask* multSelectionTask = AddTaskMultSelection(kFALSE);
   	bool isMonteCarlo = kTRUE;
	bool isMonteCarlo = 1;
       

  TString taskName;
  // gROOT->LoadMacro("AliHelperPID.cxx++g");


  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TOF");
  gSystem->AddIncludePath("-I$ALICE_ROOT/ANALYSIS");

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");

  gROOT->LoadMacro("AliAnalysisTaskFluctMCTEPOS.cxx++g");

  
  gROOT->LoadMacro("AddTaskMC.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  // gROOT->ProcessLine(".L AddTaskFluct.C");
  //Int_t nbin = 0; // will contain the number of centrality bins
AliCentralitySelectionTask *centralityTask = AddTaskCentrality();
 
 
      AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection(kTRUE);
      centralityTask->SetMCInput();

  
  AliAnalysisTaskFluctMCTEPOS *task =  AddTaskMC("lambdak0TEST",isMonteCarlo); 	
 
  //-------------------------------------------------------------------
  
  //
  // Run the analysis
  //    
  // if(chainESD) printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainESD->GetEntries());
  
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
   if(analysisMode!="proof") {
     mgr->StartAnalysis(analysisMode.Data(),nentries,firstentry);
  } else {
     // proof
    mgr->StartAnalysis(analysisMode.Data(),dataset.Data(),nentries,firstentry);
    }
  
  return;
}
//_____________________________________________________________________________
//
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="test",Bool_t useParFiles=kFALSE)
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(pluginmode.Data());
   plugin->SetUser("sbaidyan");
   // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
   /* plugin->SetROOTVersion("v5-34-30-1");
      plugin->SetAliROOTVersion("v5-06-38");*/
     plugin->SetAliPhysicsVersion("vAN-20200731-1");

   plugin->SetNtestFiles(3);
   // gROOT->LoadMacro("AddGoodRuns.C");

    
   plugin->SetGridDataDir("/alice/sim/2016/LHC16d2/");   
   plugin->SetDataPattern("/*/AliESDs.root");
   //   plugin->SetTreeName("TE");

   //plugin->AddRunNumber( 11 ); // take the most central bImp bin separately
   //for (Int_t i=2; i<4; i++) // take other bins
       //plugin->AddRunNumber( i );
   plugin->AddRunNumber(245064);
  
  
  
   plugin->SetGridWorkingDir("EPOSGen");
   // Name of executable
   plugin->SetExecutable("myCVEanalysis.sh");
   // Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  

   plugin->SetAnalysisSource("AliAnalysisTaskFluctMCTEPOS.cxx");

   plugin->SetAdditionalLibs("libPWGTools.so libPWGLFSTRANGENESS.so AliAnalysisTaskFluctMCTEPOS.h AliAnalysisTaskFluctMCTEPOS.cxx");
 

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
     //  plugin->EnablePackage("PWGHFhfe.par");
     //plugin->EnablePackage("PWGHFcorrelationHF.par");
     

   }
    plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_PHYSICS -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/PWGCF -I$ALICE_PHYSICS/PWGCF/FEMTOSCOPY -I$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/AliFemto -I$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/AliFemtoUser -g");

  
   plugin->SetDefaultOutputs(kTRUE);
   // merging via jdl
    plugin->SetMergeViaJDL(kFALSE);
   plugin->SetOneStageMerging(kFALSE);
     plugin->SetMaxMergeStages(2);

   plugin-> SetSplitMaxInputFileNumber(50);
   
   plugin->SetAnalysisMacro("AnalysisHF.C");
 
   //  plugin->SetJDLName("TaskHF.jdl");

   return plugin;
}
