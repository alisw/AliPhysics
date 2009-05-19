class AliAnalysisGrid;

void RunAnalysisAODVertexingHF()
{
  //
  // Test macro for AliAnalysisTaskSE's for heavy-flavour candidates
  // It has the structure of a Analysis Train:
  // - in this macro, change things related to running mode
  //   and input preparation 
  // - add your task using a AddTaskXXX macro 
  //
  // A.Dainese, andrea.dainese@lnl.infn.it
  // "grid" mode added by R.Bala, bala@to.infn.it
  //

  //
  TString analysisMode = "grid"; // "local", "grid", or "proof" (not yet)
  TString inputMode    = "list"; // "list", "xml", or "dataset" (not yet)
  Long64_t nentries=1234567890,firstentry=0;
  Bool_t useParFiles=kTRUE;
  Bool_t useAlienPlugin=kTRUE;
  TString pluginmode="full";
  TString loadMacroPath="$ALICE_ROOT/PWG3/vertexingHF/";
  //

  if(analysisMode=="grid") {
    // Connect to AliEn
    TGrid::Connect("alien://");
  } else if(analysisMode=="proof") {
    // Connect to the PROOF cluster
    printf("PROOF mode not yet functional..\n");
    return;
    TProof::Open("alicecaf");
    //TProof::Reset("alicecaf");
  }


  // AliRoot libraries
  if(analysisMode=="local" || analysisMode=="grid") {
    TString loadLibraries="LoadLibraries.C"; loadLibraries.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(loadLibraries.Data());
    LoadLibraries(useParFiles);
  } else if (analysisMode=="proof") {
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libVMC.so");    
    // Enable the needed packages
    //gProof->ClearPackages();
    if(!useParFiles) {
      gProof->UploadPackage("AF-v4-16");
      gProof->EnablePackage("AF-v4-16");
    } else {
      TString parDir="/afs/cern.ch/user/d/dainesea/code/";
      TString parFile;
      // --- Enable the STEERBase Package
      parFile="STEERBase.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("STEERBase");
      // --- Enable the ESD Package
      parFile="ESD.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("ESD");
      // --- Enable the AOD Package
      parFile="AOD.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("AOD");
      // --- Enable the ANALYSIS Package
      parFile="ANALYSIS.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("ANALYSIS");
      // --- Enable the ANALYSISalice Package
      parFile="ANALYSISalice.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("ANALYSISalice");
      // --- Enable the CORRFW Package
      parFile="CORRFW.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("CORRFW");
      // --- Enable the PWG3vertexingHF Package
      parFile="PWG3vertexingHF.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("PWG3vertexingHF");
    }
    gProof->ShowEnabledPackages(); // show a list of enabled packages
  }


  // Create Alien plugin, if requested
  if(useAlienPlugin) {  
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,useParFiles);  
    if(!alienHandler) return;
  }


  //-------------------------------------------------------------------
  // Prepare input chain
  TChain *chainAOD = 0;
  TString dataset; // for proof

  if(!useAlienPlugin) {
    TString makeAODInputChain="MakeAODInputchain.C"; loadLibraries.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(makeAODInputChain.Data());
    if(inputMode=="list") {
      // Local files
      chainAOD = MakeAODInputChain();// with this it reads ./AliAOD.root and ./AliAOD.VertexingHF.root
      chainAOD = MakeAODInputChain("alien:///alice/cern.ch/user/r/rbala/newtrain/out_lhc08x/180100/",1,1);
    } else if(inputMode=="xml") {
      // xml
      chainAOD = MakeAODInputChain("collection_aod.xml","collection_aodHF.xml");
    } else if(inputMode=="dataset") {
      // CAF dataset
      //gProof->ShowDataSet();
      dataset="/ITS/dainesea/AODVertexingHF_LHC08x_180100_small";
      chainAOD = MakeAODInputChainCAF(dataset.Data());
    }
  }

  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  // Connect plug-in to the analysis manager
  if(useAlienPlugin) mgr->SetGridHandler(alienHandler);


  // Input
  AliAODInputHandler *inputHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(inputHandler);
  //-------------------------------------------------------------------

  
  //-------------------------------------------------------------------
  // Analysis tasks (wagons of the train)   
  //
  TString taskName;

  taskName="AddTaskCompareHF.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAnalysisTaskSECompareHF *cmpTask = AddTaskCompareHF();

  //taskName="AddTaskSelectHF.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliAnalysisTaskSESelectHF *seleTask = AddTaskSelectHF();

  //taskName="AddTaskBkgLikeSign.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliAnalysisTaskSEBkgLikeSignJPSI *lsTask = AddTaskBkgLikeSign();

  //taskName="AddTaskBtoJPSItoEle.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliAnalysisTaskSEBtoJPSItoEle *jpsiTask = AddTaskBtoJPSItoEle();

  //taskName="AddTaskCF.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliCFHeavyFlavourTask *cfTask = AddTaskCF();

  //taskName="AddTaskCFMultiVar.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliCFHeavyFlavourTaskMultiVar *cfmvTask = AddTaskCFMultiVar();

  //taskName="AddTaskMultiVarMultiStep.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliCFHeavyFlavourTaskMultiVarMultiStep *cfmvmsTask = AddTaskCFMultiVarMultiStep();


  //-------------------------------------------------------------------

  //
  // Run the analysis
  //    
  if(chainAOD) printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainAOD->GetEntries());

  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
  mgr->StartAnalysis(analysisMode.Data(),chainAOD,nentries,firstentry);

  return;
}
//_____________________________________________________________________________
//
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="test",Bool_t useParFiles=kFALSE)
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(pluginmode.Data());
   plugin->SetUser("dainesea");
   plugin->SetNtestFiles(1);
   // Set versions of used packages
   plugin->SetAPIVersion("V2.4");
   plugin->SetROOTVersion("v5-23-02");
   plugin->SetAliROOTVersion("v4-17-00");
   // Declare input data to be processed.
   // Method 1: Create automatically XML collections using alien 'find' command.
   // Define production directory LFN
   plugin->SetGridDataDir("/alice/cern.ch/user/r/rbala/newtrain/out_lhc08x/");
   // Set data search pattern
   plugin->SetDataPattern("AliAOD.root");
   plugin->SetFriendChainName("AliAOD.VertexingHF.root");
   // ...then add run numbers to be considered
   plugin->AddRunNumber(180100);
   // Method 2: Declare existing data files (raw collections, xml collections, root file)
   // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
   // XML collections added via this method can be combined with the first method if
   // the content is compatible (using or not tags)
   //plugin->AddDataFile("/alice/cern.ch/user/r/rbala/newtrain/collection/collection_aod_lhc08w.xml");
   //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
   // Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("work");
   // Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
   // Declare the analysis source files names separated by blancs. To be compiled runtime
   // using ACLiC on the worker nodes.
   //plugin->SetAnalysisSource("$ALICE_ROOT/PWG3/vertexingHF/AliAnalysisTaskSECompareHF.cxx");
   // Declare all libraries (other than the default ones for the framework. These will be
   // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->SetAdditionalLibs("libPWG3vertexingHF.so");
   // use par files
   if(useParFiles) {
     plugin->EnablePackage("STEERBase.par");
     plugin->EnablePackage("ESD.par");
     plugin->EnablePackage("AOD.par");
     plugin->EnablePackage("ANALYSIS.par");
     plugin->EnablePackage("ANALYSISalice.par");
     plugin->EnablePackage("CORRFW.par");
     plugin->EnablePackage("PWG3vertexingHF.par");
   }
   // Declare the output file names separated by blancs.
   // (can be like: file.root or file.root@ALICE::Niham::File)
   plugin->SetOutputFiles("CmpHF.root");
   // Optionally define the files to be archived.
   //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
   // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("AnalysisHF.C");
   // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(5);
   // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
   //plugin->SetMaxInitFailed(5);
   // Optionally resubmit threshold.
   //plugin->SetMasterResubmitThreshold(90);
   // Optionally set time to live (default 30000 sec)
   //plugin->SetTTL(20000);
   // Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
   // Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("TaskHF.jdl");
   // Optionally modify job price (default 1)
   //plugin->SetPrice(1);      
   // Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   // Optionally set the preferred SE    
   plugin->SetPreferedSE("ALICE::Legnaro::SE");
   
   return plugin;
}
