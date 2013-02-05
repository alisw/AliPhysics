Bool_t RunAnalysisTaskMuonHadronCorrelations(Int_t runNumber = 188362, const char * runMode="full", TString centMethod="V0M") {
  
  TString rootVersion = "v5-34-02-1";
  TString alirootVersion = "v5-04-22-AN";
  gSystem->AddIncludePath("-I$ALICE_ROOT/include ");
  
  gSystem->Load("libTree.so")          ;
  gSystem->Load("libGeom.so")          ;
  gSystem->Load("libVMC.so")           ;
  gSystem->Load("libMinuit.so")        ;
  gSystem->Load("libPhysics.so")       ;
  gSystem->Load("libSTEERBase.so")     ;
  gSystem->Load("libESD.so")           ;
  gSystem->Load("libAOD.so")           ;
  gSystem->Load("libANALYSIS.so")      ;
  gSystem->Load("libOADB.so")          ;
  gSystem->Load("libANALYSISalice.so") ;
  gSystem->Load("libCORRFW.so")        ;
  gSystem->Load("libPWGmuon.so")       ;
  
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
  //  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  //  AddTaskPhysicsSelection(0);
  mgr->AddStatisticsTask(AliVEvent::kAny);
  
  // -------------------------------
  // MUON TRACK CUTS CONFIGURATION
  // -------------------------------

  gROOT->LoadMacro("AliAnalysisTaskMuonHadronCorrelations.cxx++g");
  AliAnalysisTaskMuonHadronCorrelations *task = new AliAnalysisTaskMuonHadronCorrelations("AliAnalysisTaskMuonHadronCorrelations");

  // Set analysis cuts   
  task->SetFilterBitCentralBarrel(7);  // -> 128
  task->SetMaxEtaCentralBarrel(1.0);
  task->SetTriggerMatchLevelMuon(1);

  const Int_t nBinCent = 4;
  Double_t centLimits[nBinCent+1] = {0., 20., 40, 60., 100.};
  task->SetCentBinning(nBinCent, centLimits);

  task->SetCentMethod(centMethod);

  const Int_t nBinPt = 3;
  Double_t ptLimits[nBinPt+1] = {0., 1., 2., 4.};
  task->SetPtBinning(nBinPt, ptLimits);

  mgr->AddTask(task);

  // create output container
  AliAnalysisDataContainer *output1 = mgr->CreateContainer("list", TList::Class(), AliAnalysisManager::kOutputContainer, "MuonHadronCorrelations.root");
  
  // finaly connect input and output
  mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output1);
    
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
