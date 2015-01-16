Bool_t RunAnalysisTaskMuonHadronCorrelations(Int_t runNumber = 188362, const char * runMode="full", TString centMethod="V0M") {
  
  TString rootVersion = "v5-34-02-1";
  TString alirootVersion = "v5-04-22-AN";
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
  
  // -------------------------------
  // ANALYSIS MANAGER
  // -------------------------------

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr) delete mgr;
  mgr = new AliAnalysisManager("AM","Analysis Manager");
  
  // -------------------------------
  // INPUT EVENT HANDLER
  // -------------------------------

  AliAODInputHandler* inputHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(inputHandler);
  
  // -------------------------------
  // ANALYSIS PLUGIN CONFIGURATION
  // -------------------------------

  AliAnalysisAlien *analysisPlugin = new AliAnalysisAlien();
  if (!analysisPlugin) { Printf("Error : analysisPlugin is null !!!"); return kFALSE; }
  analysisPlugin->SetAPIVersion("V1.1x");
  //analysisPlugin->SetROOTVersion(rootVersion.Data());
  analysisPlugin->SetAliROOTVersion(alirootVersion.Data());
  analysisPlugin->SetExecutableCommand("aliroot -b -q");

  // Overwrite all generated files, datasets and output results from a previous session
  analysisPlugin->SetOverwriteMode();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  analysisPlugin->SetRunMode(runMode);  // VERY IMPORTANT 
  
  analysisPlugin->SetAdditionalRootLibs("CORRFW");
  analysisPlugin->SetAdditionalRootLibs("PWGmuon");
  
  // -------------------------------
  // GRID CONFIG
  // -------------------------------

  // Location of Data
  analysisPlugin->SetGridDataDir("/alice/data/2012/LHC12g");
  analysisPlugin->SetDataPattern("*ESDs/pass2/AOD/*/AliAOD.root");
  analysisPlugin->SetRunPrefix("000");   // real data

  analysisPlugin->AddRunNumber(runNumber);

  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  //   	   analysisPlugin->SetGridWorkingDir("LHC11h_MUON");
  analysisPlugin->SetGridWorkingDir("MuonHadronCorrelations_LHC12g");
  
  // Declare alien output directory. Relative to working directory.
  analysisPlugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output

  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  analysisPlugin->SetAnalysisSource("AliAnalysisTaskMuonHadronCorrelations.cxx");

  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  analysisPlugin->SetAdditionalLibs("libCORRFW.so libPWGmuon.so AliAnalysisTaskMuonHadronCorrelations.h AliAnalysisTaskMuonHadronCorrelations.cxx");
  
  analysisPlugin->AddIncludePath("-I.");

  //les repertoires individuels d'output des jobs sur la grille ont comme nom le numero du run
  analysisPlugin->SetOutputToRunNo();
  
  // Optionally define the files to be archived.
  //   analysisPlugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
  //   analysisPlugin->SetOutputArchive("log_archive.zip:stdout,stderr");

  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  analysisPlugin->SetAnalysisMacro("MuonHadronCorr.C");

  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  //	   analysisPlugin->SetSplitMaxInputFileNumber(20);
  analysisPlugin->SetSplitMaxInputFileNumber(40);
  // Number of runs per master job
  analysisPlugin->SetNrunsPerMaster(1);
  
  // Optionally modify the executable name (default analysis.sh)
  analysisPlugin->SetExecutable("MuonHadronCorr.sh");

  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  //   analysisPlugin->SetMaxInitFailed(5);

  // Optionally resubmit threshold.
  //   analysisPlugin->SetMasterResubmitThreshold(90);

  // Optionally set time to live (default 30000 sec)
  analysisPlugin->SetTTL(30000);
  // Optionally set input format (default xml-single)
  analysisPlugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  analysisPlugin->SetJDLName("MuonHadronCorr.jdl");
  // Optionally modify job price (default 1)
  analysisPlugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  analysisPlugin->SetSplitMode("se");
  
  analysisPlugin->SetNtestFiles(5);
  //   analysisPlugin->SetMergeViaJDL(1);
  analysisPlugin->SetOverwriteMode(kTRUE);
  
  mgr->SetGridHandler(analysisPlugin);
  
  // -------------------------------
  // PHYSICS SELECTION
  // -------------------------------
  //  gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
  //  AddTaskPhysicsSelection(0);
  mgr->AddStatisticsTask(AliVEvent::kAny);
  
  gROOT->LoadMacro("AddAnalysisTaskMuonHadronCorrelations.C");
  AddAnalysisTaskMuonHadronCorrelations();

  TStopwatch timer;
  timer.Start();
  // -------------------------------
  // RUN ANALYSIS
  // -------------------------------
  mgr->InitAnalysis();
  mgr->PrintStatus();
  
  mgr->StartAnalysis("grid");
  
  timer.Stop();
  timer.Print();
  return kTRUE;

};
