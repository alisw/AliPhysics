//====================================================================================================================================================

Bool_t RunAnalysisTaskMFTExample(const Char_t *runType="local",       // "grid" and "local" types have been tested
				 const Char_t *runMode="full") {
  
  //  enum {kGenerated, kReconstructed};

  LoadLibs();

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr) delete mgr;
  mgr = new AliAnalysisManager("AM","Analysis Manager");
  
  AliAODInputHandler* inputHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(inputHandler);
  
  if (!strcmp(runType,"grid")) mgr->SetGridHandler(CreateAlienHandler(runMode));
  
  gROOT->LoadMacro("AliAnalysisTaskMFTExample.cxx++g");   
  AliAnalysisTaskMFTExample *task = new AliAnalysisTaskMFTExample("AliAnalysisTaskMFTExample");

  // in cm, taken from Fig. 7.4 of the ITS Upgrade TDR, in the hypothesis of ~80 tracks participating to the vtx
  task -> SetVtxResolutionITS(5.e-4, 5.e-4, 4.e-4);
  task -> SetVertexMode(AliAnalysisTaskMFTExample::kGenerated);

  mgr->AddTask(task);
  
  // create output container(s)
  AliAnalysisDataContainer *histogramList = mgr->CreateContainer("list", TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisMFT_Output.root");
  
  // connect input and output
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, histogramList);

  // RUN ANALYSIS

  TStopwatch timer;
  timer.Start();
  mgr->InitAnalysis();
  mgr->PrintStatus();
  
  if (!strcmp(runType,"grid")) {
    printf("Starting MFT analysis on the grid");
    mgr->StartAnalysis("grid");
  }

  else if (!strcmp(runType,"local")) {
    printf("Starting MFT analysis locally");
    mgr->StartAnalysis("local", GetInputLocalData());
  }

  else AliError(Form("Specified run type %s is not supported", runType));
  
  timer.Stop();
  timer.Print();
  return kTRUE;

}

//====================================================================================================================================================

AliAnalysisGrid* CreateAlienHandler(const Char_t *runMode) {

  // Set up the analysis plugin in case of a grid analysis

  TString rootVersion = "v5-34-08-6";
  TString alirootVersion = "vAN-20140727";

  AliAnalysisAlien *analysisPlugin = new AliAnalysisAlien();
  if (!analysisPlugin) { Printf("Error : analysisPlugin is null !!!"); return kFALSE; }
  analysisPlugin->SetAPIVersion("V1.1x");
  analysisPlugin->SetROOTVersion(rootVersion.Data());
  analysisPlugin->SetAliROOTVersion(alirootVersion.Data());
  analysisPlugin->SetExecutableCommand("aliroot -b -q");

  // Overwrite all generated files, datasets and output results from a previous session
  analysisPlugin->SetOverwriteMode();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  analysisPlugin->SetRunMode(runMode);  // VERY IMPORTANT 
  
  analysisPlugin->SetAdditionalRootLibs("CORRFW");
  analysisPlugin->SetAdditionalRootLibs("PWGmuon");
  
  // Location of Data and Working dir
  analysisPlugin->SetGridDataDir("/alice/cern.ch/user/a/auras/MFT/simulations_2014/PbPb/jpsi_prompt/pix20um20um_plane400um/");
  analysisPlugin->SetDataPattern("*/AliAOD.Muons.root");
  analysisPlugin->SetRunPrefix("");
  analysisPlugin->SetGridWorkingDir("MFT/analysis_2014/PbPb/jpsi_prompt/pix20um20um_plane400um/");
  
  // Declare alien output directory. Relative to working directory.
  analysisPlugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime using ACLiC on the worker nodes.
  analysisPlugin->SetAnalysisSource("AliAnalysisTaskMFTExample.cxx");
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  analysisPlugin->SetAdditionalLibs("libCORRFW.so libPWGmuon.so AliAnalysisTaskMFTExample.h AliAnalysisTaskMFTExample.cxx");
  
  analysisPlugin->AddIncludePath("-I.");

  analysisPlugin->SetOutputToRunNo();
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  analysisPlugin->SetAnalysisMacro("MFTAnalysis.C");

  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  analysisPlugin->SetSplitMaxInputFileNumber(100);
  // Number of runs per master job
  analysisPlugin->SetNrunsPerMaster(1);
  
  // Optionally modify the executable name (default analysis.sh)
  analysisPlugin->SetExecutable("MFTAnalysis.sh");

  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  analysisPlugin->SetMaxInitFailed(5);
  // Optionally resubmit threshold.
  analysisPlugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  analysisPlugin->SetTTL(30000);
  // Optionally set input format (default xml-single)
  analysisPlugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  analysisPlugin->SetJDLName("MFTAnalysis.jdl");
  // Optionally modify job price (default 1)
  analysisPlugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  analysisPlugin->SetSplitMode("se");
  
  analysisPlugin->SetNtestFiles(5);
  //   analysisPlugin->SetMergeViaJDL(1);
  analysisPlugin->SetOverwriteMode(kTRUE);
  
  return analysisPlugin;

}

//====================================================================================================================================================

TChain* GetInputLocalData() {

  // Set up the chain of input events in case of a local analysis

  TChain *chain = new TChain("aodTree");
  chain->Add("./AliAOD.Muons.root");

  return chain;

}

//====================================================================================================================================================

void LoadLibs() {

  gSystem->AddIncludePath("-I$ALICE_ROOT/include ");
  
  gSystem->Load("libTree")          ;
  gSystem->Load("libGeom")          ;
  gSystem->Load("libVMC")           ;
  gSystem->Load("libMinuit")        ;
  gSystem->Load("libPhysics")       ;
  gSystem->Load("libSTEERBase")     ;
  gSystem->Load("libESD")           ;
  gSystem->Load("libAOD")           ;
  gSystem->Load("libANALYSIS")      ;
  gSystem->Load("libOADB")          ;
  gSystem->Load("libANALYSISalice") ;
  gSystem->Load("libCORRFW")        ;
  gSystem->Load("libPWGmuon")       ;

}

//====================================================================================================================================================
