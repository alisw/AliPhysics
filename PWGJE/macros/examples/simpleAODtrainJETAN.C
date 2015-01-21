AliAnalysisGrid* CreateAlienHandler(const Char_t *mode);

void simpleAODtrainJETAN(const Char_t *mode="full"){
  //
  //  Run analysis, mode can be "full", "test", "submit" or "terminate",
  //           see CreateAlienHandler for details.
  //

  // Load common libraries
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   

  gSystem->Load("libCORRFW");
  gSystem->Load("libJETAN");
  gSystem->Load("libfastjet");
  gSystem->Load("libSISConePlugin");
  gSystem->Load("libFASTJETAN");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGJE");

  // Use AliRoot includes to compile our task
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/JETAN -I$ALICE_PHYSICS/PWG4/JetTasks");

  // Create and configure the alien handler plugin
  AliAnalysisGrid *alienHandler = CreateAlienHandler(mode);  
  if (!alienHandler) return;

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);

  // Standard train runs on ESD; AOD may not work?
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  // Need output handler, but tree is not filled...
  gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddAODOutputHandler.C"));
  AliVEventHandler* handler = AddAODOutputHandler();
  //handler->SetFillAODForRun(kFALSE);

  mgr->SetCommonFileName("AnalysisResult.root");

  // Only for ESD?
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  //AddTaskPhysicsSelection(kFALSE, kTRUE);


  // gROOT->Macro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C"); // Only for ESD

  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/ConfigLegoTrainPWGJE.C");
  ConfigLegoTrainPWGJE(1108); // For DeltaAODName; LHC11h settings

  const Int_t kHighPtFilterMask = 768;  // LHC11h
  const Int_t kHighPtFilterMaskBest = 256; // LHC11h
  const Float_t kR = 0.3;

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/AddTaskJetCluster.C");
  jetclu = AddTaskJetCluster("AOD","",kHighPtFilterMask,AliVEvent::kAny,"ANTIKT",kR,0,kTRUE, AliAnalysisManager::GetGlobalStr("kJetDeltaAODName",gDebug),0.15,0.9,10,0);

  jetclu->SetEventSelection(kTRUE); 
  jetclu->SetJetTypes(1<<0|1<<2|1<<3);
  jetclu->SetNRandomCones(1);
  jetclu->SetCentralityCut(0.,80.);
  jetclu->SetJetOutputMinPt(0);
  jetclu->SetFilterMaskBestPt(kHighPtFilterMaskBest);
  //jetclu->SetDebugLevel(1);

  // output dir: PWGJE_cluster_AOD__ANTIKT03_B0_Filter00272_Cut00150_Skip00

  jetclukt = AddTaskJetCluster("AOD","",kHighPtFilterMask,AliVEvent::kAny,"KT",kR,0,kTRUE,AliAnalysisManager::GetGlobalStr("kJetDeltaAODName",gDebug),0.15,0.9,10,0);
  jetclukt->SetBackgroundCalc(kTRUE);
  jetclukt->SetEventSelection(kTRUE); 
  jetclukt->SetCentralityCut(0.,80.);
  jetclukt->SetJetOutputMinPt(0);
  jetclukt->SetFilterMaskBestPt(kHighPtFilterMaskBest);
  //jetclukt->SetDebugLevel(1);

  AliAnalysisManager::SetGlobalInt("kPhysicsSelectionFlag",AliVEvent::kMB|AliVEvent::kCentral|AliVEvent::kSemiCentral|AliVEvent::kEMCEJE|AliVEvent::kEMCEGA|AliVEvent::kEMC1);

  Int_t kRpar = Int_t(10*kR + 0.001);
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/AddTaskJetBackgroundSubtract.C");
  jetbkg = AddTaskJetBackgroundSubtract(TString(Form("clustersAOD_ANTIKT%02d_B0_Filter%05d_Cut00150_Skip00",kRpar,kHighPtFilterMask)),2);
  jetbkg->SetBackgroundBranch(Form("jeteventbackground_clustersAOD_KT%02d_B0_Filter%05d_Cut00150_Skip00",kRpar,kHighPtFilterMask));
  jetbkg->SelectCollisionCandidates(AliAnalysisManager::GetGlobalInt("kPhysicsSelectionFlag", gDebug));
  jetbkg->SetKeepJets(kTRUE);
  jetbkg->SetNonStdOutputFile(AliAnalysisManager::GetGlobalStr("kJetDeltaAODName",gDebug));
  //jetbkg->SetDebugLevel(1);

  // Looking for input dir: clustersAOD_ANTIKT03_B2_Filter00768_Cut00150_Skip00
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/AddTaskFragmentationFunction.C");
  Int_t kEvtClass = 1; // Central events

  // Need to add 'type = AOD' as third argument -- only for MC...
  fftask = AddTaskFragmentationFunction("clustersAOD_ANTIKT", "", "", "", "", kHighPtFilterMask, kR,2,150.,kEvtClass, "_Skip00");
  //fftask->SetPhiCorrHistoBins(); 
  fftask->SetEventSelectionMask(AliVEvent::kMB|AliVEvent::kCentral|AliVEvent::kSemiCentral);
  fftask->UseAODInputJets(kFALSE);  // to pick up jets from the output AOD
  //fftask->SetDebugLevel(1);


  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  // Start analysis via AlienPlugin; test runs etc are handled via flag to CreateAlienHandler
  mgr->StartAnalysis("grid"); 
};


AliAnalysisGrid* CreateAlienHandler(const Char_t *mode)
{

   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   plugin->SetOverwriteMode();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(mode); //terminate");
// Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-34-08");
   plugin->SetAliROOTVersion("v5-05-11-AN");
   plugin->AddExternalPackage("cgal::v3.6 fastjet::v2.4.2");

   // Declare input data to be processed.
   // Method 1: Create automatically XML collections using alien 'find' command.
   // Define production directory LFN
   // One file:
   //plugin->SetGridDataDir("/alice/data/2011/LHC11h_2/000170081/ESDs/pass2/AOD115/0829/");
   // One run:
   plugin->SetGridDataDir("/alice/data/2011/LHC11h_2/000170081/ESDs/pass2/AOD115/");
   // Set data search pattern
   plugin->SetDataPattern("*AOD.root");

   // Can be extended to multiple runs:
   //
   //   plugin->SetRunPrefix("000");   // real data
   // ...then add run numbers to be considered
   //plugin->AddRunNumber(191027);
   //   plugin->AddRunNumber(104065);  // real data
   //   plugin->SetOutputSingleFolder("output");
   //   plugin->SetOutputToRunNo();

   
   // Method 2: Declare existing data files (raw collections, xml collections, root file)
   // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
   // XML collections added via this method can be combined with the first method if
   // the content is compatible (using or not tags)
   //   plugin->AddDataFile("tag.xml");
   //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");


   // Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("AOD_jets_PbPb");
   // Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
   // Declare the analysis source files names separated by blancs. To be compiled runtime
   // using ACLiC on the worker nodes.
   // plugin->SetAnalysisSource("AliAnalysisAODJetHists.cxx");
   // Declare all libraries (other than the default ones for the framework. These will be
   // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->SetAdditionalLibs("libCORRFW.so libJETAN.so libCGAL.so libfastjet.so libSISConePlugin.so libFASTJETAN.so libPWGTools.so libPWGJE.so");
   plugin->AddIncludePath("-I$ALICE_ROOT/JETAN -I$ALICE_PHYSICS/PWGJE");
   // Declare the output file names separated by blancs.
   // (can be like: file.root or file.root@ALICE::Niham::File)
   //   plugin->SetOutputFiles("Pt.ESD.1.root");
   plugin->SetDefaultOutputs(kFALSE);
   plugin->SetOutputFiles("AnalysisResult.root EventStat_temp.root");
   //plugin->SetDefaultOutputs();
   // Optionally define the files to be archived.
   //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
   //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
   // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("AODJetsPbPb.C");
   // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(100);
   // Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable("AODJetsPbPb.sh");
   // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
   //   plugin->SetMaxInitFailed(5);
   // Optionally resubmit threshold.
   //   plugin->SetMasterResubmitThreshold(90);
   // Optionally set time to live (default 30000 sec)
   plugin->SetTTL(30000);
   // Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
   // Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("AODJetsPbPb.jdl");
   // Optionally modify job price (default 1)
   plugin->SetPrice(1);      
   // Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   return plugin;
}
