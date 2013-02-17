class AliAnalysisAlien;

void runGridLambdaOverK0sJets(TString  runMode    = "test", 
			      TString  alirootVer = "v5-04-31LF-AN"/*"v5-04-14-AN"*/,
			      TString  rootVer    = "v5-34-02-1",
			      TString  dataPath   = "ESDs/pass2/AOD086/*/AliAOD.root",
			      TString  dataDir    = "/alice/data/2010/LHC10h",
			      TString  workDir    = "test",
			      TString  name       = "LambdaOverK0sRatio", 
			      Int_t    data       = 2010,
			      Float_t  minCen     = 0.,
			      Float_t  maxCen     = 90.,
			      Float_t  ptMinTrig  = 8.,
			      Float_t  ptMaxTrig  = 20.,
			      Float_t  etaMaxTrig = 0.75,
			      Float_t  rapMaxV0   = 0.75,
			      Bool_t   sepInjec   = kTRUE,
			      Bool_t   isMC       = kFALSE,
			      Bool_t   usePID     = kTRUE,
			      Bool_t   doQA       = kFALSE,
			      Int_t    run        = 137530/*138624*/){
  
  Printf("   \nThe parameters of the programm are : \n ");
  Printf(" \t Analysis mode:\t %s\n \t Centrality:\t %.1lf - %.1lf\n \t Use MC Data?:\t %s\n \t Use PID?:\t %s\n",
	 "Grid",minCen,maxCen,
	 (isMC) ? "Yes" : "No",
	 (usePID) ? "Yes" : "No");
  
  // _____________________________________________________ //
  
  InitAndLoadLibs();
  
  AliAnalysisManager *mgr = new AliAnalysisManager("Manager");
  
  AliAnalysisGrid *alienHandler = CreateAlienHandler(runMode,alirootVer,rootVer,dataPath,dataDir,workDir,isMC,run);  
  if (!alienHandler) return;
  mgr->SetGridHandler(alienHandler);
    
  AliAODInputHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  
  //PID
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTask *pidTask = AddTaskPIDResponse(isMC);
  //AliAnalysisTask *pidTask = AddTaskPIDResponse(isMC,kTRUE);
  if(!pidTask) { printf("no PIDtask\n"); return; }
 
  Float_t checkIDTrig= kTRUE;

  // My task
  gROOT->LoadMacro("AliAnalysisTaskLambdaOverK0sJets.cxx++g"); 
  //gSystem->Load("libPWGLFSTRANGENESS");
  gROOT->LoadMacro("AddTaskLambdaOverK0sJets.C");
  AliAnalysisTaskLambdaOverK0sJets *task = AddTaskLambdaOverK0sJets(name,data,minCen,maxCen,ptMinTrig,ptMaxTrig,etaMaxTrig,checkIDTrig,rapMaxV0,sepInjec,isMC,usePID,doQA);
   // _____________________________________________________ //
 
   if (!mgr->InitAnalysis()) return;
   mgr->PrintStatus();
   
   mgr->StartAnalysis("grid");
}

// ______________________________________________________________

void InitAndLoadLibs() {
  
    gSystem->Load("libCore.so"); 
    gSystem->Load("libTree.so");                 
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics");
    gSystem->Load("libMinuit.so");  
    gSystem->Load("libProof.so");
    gSystem->Load("libGui.so");
    gSystem->Load("libXMLParser.so");
    gSystem->Load("libProofPlayer.so");
    gSystem->Load("libXMLIO.so");

    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libCDB.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libCORRFW.so");
    gSystem->Load("libJETAN.so");
    gSystem->Load("libRAWDatabase.so");
    gSystem->Load("libSTEER.so");
    gSystem->Load("libCORRFW.so");

   
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
}

// ___________________________________________________________________ //

// **** It is change 'AliAnalysisGrid' by 'AliAnalysisAlien'
AliAnalysisAlien* CreateAlienHandler(TString runMode,TString alirootVer,
				     TString rootVer,TString dataPath,
				     TString dataDir,TString workDir,
				     Bool_t isMC,Int_t kRun) {
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetCheckCopy(kFALSE);

  plugin->SetRunMode(runMode);
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion(rootVer);
  plugin->SetAliROOTVersion(alirootVer);
  
  /////////////////////////////////////////////////////////////////

  // Declare input data to be processed.
  // Method 1: Create automatically XML collections using alien 'find' command.
  plugin->SetGridDataDir(dataDir);
  
  if (!isMC) 
    plugin->SetRunPrefix("000");
  plugin->SetDataPattern(dataPath);
    
  plugin->AddRunNumber(kRun);

  // Method 2: Declare existing data files (raw collections, xml collections, root file)
  const char working_dir[250];

  sprintf(working_dir, "%s/%d",workDir.Data(),kRun);
  TString path = TString(working_dir);
  plugin->SetGridWorkingDir(path);

  //plugin->SetGridWorkingDir(workDir);
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  plugin->SetAnalysisSource("AliAnalysisTaskLambdaOverK0sJets.cxx");
  //plugin->SetAdditionalLibs("AliAnalysisTaskMultiplicity.h AliAnalysisTaskMultiplicity.cxx");
  plugin->SetAdditionalLibs("AliAnalysisTaskLambdaOverK0sJets.cxx AliAnalysisTaskLambdaOverK0sJets.h");
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("AnalysisProduction.C");
  //plugin->SetAnalysisMacro("mytask.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(50);
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  plugin->SetMaxInitFailed(12);
  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("TaskProduction.jdl");
  //plugin->SetJDLName("mytask.jdl");
  plugin->SetMergeViaJDL(kTRUE);
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");
  return plugin;

}
