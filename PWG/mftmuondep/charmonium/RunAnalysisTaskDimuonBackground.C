//====================================================================================================================================================

Bool_t RunAnalysisTaskDimuonBackground(const Char_t *runType="local",       // "grid" and "local" types have been tested
				       const Char_t *runMode="full") {
  
  //  enum {kGenerated, kReconstructed};

  LoadLibs();

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr) delete mgr;
  mgr = new AliAnalysisManager("AM","Analysis Manager");
  
  AliMultiInputEventHandler *handler = new AliMultiInputEventHandler();
  handler->AddInputEventHandler(new AliAODInputHandler());

  //-------------------- Setup of the mixed event handler ---------------
  
  const Int_t mixNum = 20;      // Number of times the main event is used for a mixing procedure
  const Int_t bufferSize = 1;   // Number of events involved in each mixing procedure with the main one (typically one)

  AliMixInputEventHandler *mixHandler = new AliMixInputEventHandler(bufferSize, mixNum);
  mixHandler->SetInputHandlerForMixing(dynamic_cast<AliMultiInputEventHandler*> handler);
  AliMixEventPool *evPool = new AliMixEventPool();

  AliMixEventCutObj *zVertex = new AliMixEventCutObj(AliMixEventCutObj::kZVertex, -10, 10, 1);

  evPool->AddCut(zVertex);
  mixHandler->SetEventPool(evPool);

  mixHandler->SelectCollisionCandidates(AliVEvent::kAny);
  handler->AddInputEventHandler(mixHandler);

  //---------------------------------------------------------------------

  mgr->SetInputEventHandler(handler);
  
  if (!strcmp(runType,"grid")) mgr->SetGridHandler(CreateAlienHandler(runMode));
  
  gROOT->LoadMacro("AliAnalysisTaskDimuonBackground.cxx++g");   
  AliAnalysisTaskDimuonBackground *task = new AliAnalysisTaskDimuonBackground("AliAnalysisTaskDimuonBackground");

  // in cm, taken from Fig. 7.4 of the ITS Upgrade TDR, in the hypothesis of ~80 tracks participating to the vtx
  // task -> SetVtxResolutionITS(5.e-4, 5.e-4, 4.e-4);

  task -> SetVertexMode(AliAnalysisTaskDimuonBackground::kReconstructed);

  task -> SetMinTriggerMatch(1);
  task -> SetSingleMuonMinPt(1.0);
  task -> SetSingleMuonMinEta(-3.6);
  task -> SetSingleMuonMaxEta(-2.5);
  task -> SetSingleMuonMaxChi2(9999.);
  
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
  TString alirootVersion = "vAN-20141128";

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
  
  // Method1: Location of Data and Working dir
  // analysisPlugin->SetGridDataDir("/alice/sim/2014/LHC14j5_new_plus/");
  // analysisPlugin->SetDataPattern("AliAOD.root");
  // analysisPlugin->AddRunNumber(137161);

  // Method2: Use existing xml collection

  analysisPlugin->AddDataFile("../../collections/137844_1.xml"); 
  analysisPlugin->AddDataFile("../../collections/137848_1.xml"); 
  analysisPlugin->AddDataFile("../../collections/138190_1.xml"); 
  analysisPlugin->AddDataFile("../../collections/138192_1.xml"); 
  analysisPlugin->AddDataFile("../../collections/138197_1.xml"); 
  analysisPlugin->AddDataFile("../../collections/138201_1.xml"); 
  analysisPlugin->AddDataFile("../../collections/138225_1.xml"); 
  analysisPlugin->AddDataFile("../../collections/138275_1.xml"); 
  analysisPlugin->AddDataFile("../../collections/138364_1.xml"); 
  analysisPlugin->AddDataFile("../../collections/138396_1.xml");

  // analysisPlugin->AddDataFile("137230_2.xml");
  // analysisPlugin->AddDataFile("137162_2.xml");
  // analysisPlugin->AddDataFile("137161_2.xml");
  // analysisPlugin->AddDataFile("137231_2.xml");
  // analysisPlugin->AddDataFile("137430_2.xml");
  // analysisPlugin->AddDataFile("137366_2.xml");
  // analysisPlugin->AddDataFile("137243_2.xml");
  // analysisPlugin->AddDataFile("137236_2.xml");
  // analysisPlugin->AddDataFile("137232_2.xml");
  // analysisPlugin->AddDataFile("137440_2.xml");
  // analysisPlugin->AddDataFile("137434_2.xml");
  // analysisPlugin->AddDataFile("137432_2.xml");
  // analysisPlugin->AddDataFile("137431_2.xml");
  // analysisPlugin->AddDataFile("137443_2.xml");
  // analysisPlugin->AddDataFile("137441_2.xml");
  // analysisPlugin->AddDataFile("137544_2.xml");
  // analysisPlugin->AddDataFile("137541_2.xml");
  // analysisPlugin->AddDataFile("137539_2.xml");
  // analysisPlugin->AddDataFile("137549_2.xml");
  // analysisPlugin->AddDataFile("137608_2.xml");
  // analysisPlugin->AddDataFile("137595_2.xml");
  // analysisPlugin->AddDataFile("137686_2.xml");
  // analysisPlugin->AddDataFile("137639_2.xml");
  // analysisPlugin->AddDataFile("137638_2.xml");
  // analysisPlugin->AddDataFile("137692_2.xml");
  // analysisPlugin->AddDataFile("137691_2.xml");
  // analysisPlugin->AddDataFile("137718_2.xml");
  // analysisPlugin->AddDataFile("137704_2.xml");
  // analysisPlugin->AddDataFile("137751_2.xml");
  // analysisPlugin->AddDataFile("137724_2.xml");
  // analysisPlugin->AddDataFile("137722_2.xml");
  // analysisPlugin->AddDataFile("137752_2.xml");
  // analysisPlugin->AddDataFile("137844_2.xml");
  // analysisPlugin->AddDataFile("138190_2.xml");
  // analysisPlugin->AddDataFile("137848_2.xml");
  // analysisPlugin->AddDataFile("138197_2.xml");
  // analysisPlugin->AddDataFile("138192_2.xml");
  // analysisPlugin->AddDataFile("138225_2.xml");
  // analysisPlugin->AddDataFile("138201_2.xml");
  // analysisPlugin->AddDataFile("138275_2.xml");
  // analysisPlugin->AddDataFile("138364_2.xml");
  // analysisPlugin->AddDataFile("138396_2.xml");
  
  // analysisPlugin->AddDataFile("137161_3.xml");
  // analysisPlugin->AddDataFile("137162_3.xml");
  // analysisPlugin->AddDataFile("137230_3.xml");
  // analysisPlugin->AddDataFile("137231_3.xml");
  // analysisPlugin->AddDataFile("137232_3.xml");
  // analysisPlugin->AddDataFile("137236_3.xml");
  // analysisPlugin->AddDataFile("137243_3.xml");
  // analysisPlugin->AddDataFile("137366_3.xml");
  // analysisPlugin->AddDataFile("137430_3.xml");
  // analysisPlugin->AddDataFile("137431_3.xml");
  // analysisPlugin->AddDataFile("137432_3.xml");
  // analysisPlugin->AddDataFile("137434_3.xml");
  // analysisPlugin->AddDataFile("137440_3.xml");
  // analysisPlugin->AddDataFile("137441_3.xml");
  // analysisPlugin->AddDataFile("137443_3.xml");
  // analysisPlugin->AddDataFile("137539_3.xml");
  // analysisPlugin->AddDataFile("137541_3.xml");
  // analysisPlugin->AddDataFile("137544_3.xml");
  // analysisPlugin->AddDataFile("137549_3.xml");
  // analysisPlugin->AddDataFile("137595_3.xml");
  // analysisPlugin->AddDataFile("137608_3.xml");
  // analysisPlugin->AddDataFile("137638_3.xml");
  // analysisPlugin->AddDataFile("137639_3.xml");
  // analysisPlugin->AddDataFile("137686_3.xml");
  // analysisPlugin->AddDataFile("137691_3.xml");
  // analysisPlugin->AddDataFile("137692_3.xml");
  // analysisPlugin->AddDataFile("137704_3.xml");
  // analysisPlugin->AddDataFile("137718_3.xml");
  // analysisPlugin->AddDataFile("137722_3.xml");
  // analysisPlugin->AddDataFile("137724_3.xml");
  // analysisPlugin->AddDataFile("137751_3.xml");
  // analysisPlugin->AddDataFile("137752_3.xml");
  // analysisPlugin->AddDataFile("137844_3.xml");
  // analysisPlugin->AddDataFile("137848_3.xml");
  // analysisPlugin->AddDataFile("138190_3.xml");
  // analysisPlugin->AddDataFile("138192_3.xml");
  // analysisPlugin->AddDataFile("138197_3.xml");
  // analysisPlugin->AddDataFile("138201_3.xml");
  // analysisPlugin->AddDataFile("138225_3.xml");
  // analysisPlugin->AddDataFile("138275_3.xml");
  // analysisPlugin->AddDataFile("138364_3.xml");
  // analysisPlugin->AddDataFile("138396_3.xml");
  
  analysisPlugin->SetGridWorkingDir("MFT/analysis/LHC14cj5_PbPb_0_5/charmonia/background");
  
  // Declare alien output directory. Relative to working directory.
  analysisPlugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime using ACLiC on the worker nodes.
  analysisPlugin->SetAnalysisSource("AliAnalysisTaskDimuonBackground.cxx");
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  analysisPlugin->SetAdditionalLibs("libCORRFW.so libPWGmuon.so libEventMixing.so AliAnalysisTaskDimuonBackground.h AliAnalysisTaskDimuonBackground.cxx");
  
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
  chain->Add("./test_03Nov2014_zVtx_Rndm/AliAOD.Muons.root");

  return chain;

}

//====================================================================================================================================================

void LoadLibs() {

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

}

//====================================================================================================================================================
