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


  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWG3 -I$ALICE_ROOT/PWG3/vertexingHF -g"); 
  //
  TString trainName = "D2H";
  TString analysisMode = "grid"; // "local", "grid", or "proof"
  TString inputMode    = "list"; // "list", "xml", or "dataset"
  Long64_t nentries=123567890,firstentry=0;
  Bool_t useParFiles=kFALSE;
  Bool_t useAlienPlugin=kTRUE;
  TString pluginmode="full";
  Bool_t saveProofToAlien=kFALSE;
  TString proofOutdir = "";
  TString loadMacroPath="$ALICE_ROOT/PWG3/vertexingHF/macros/";
  //TString loadMacroPath="./"; // this is normally needed for CAF
  //

  if(analysisMode=="grid") {
    // Connect to AliEn
    TGrid::Connect("alien://");
  } else if(analysisMode=="proof") {
    // Connect to the PROOF cluster
    if(inputMode!="dataset") {printf("Input mode must be dataset, for proof analysis\n"); return;}
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    TProof::Open("alicecaf");
    //TProof::Reset("alicecaf");
    if(saveProofToAlien) {
      TGrid::Connect("alien://");
      if(gGrid) {
	TString homedir = gGrid->GetHomeDirectory();
	TString workdir = homedir + trainName;
	if(!gGrid->Cd(workdir)) {
	  gGrid->Cd(homedir);
	  if(gGrid->Mkdir(workdir)) {
	    gGrid->Cd(trainName);
	    ::Info("VertexingTrain::Connect()", "Directory %s created", gGrid->Pwd());
	  }
	}	   
	gGrid->Mkdir("proof_output");
	gGrid->Cd("proof_output");
	proofOutdir = Form("alien://%s", gGrid->Pwd());
      } 
    }
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
    gSystem->Load("libMinuit.so");    
    // Enable the needed packages
    //gProof->ClearPackages();
    TString parDir="/afs/cern.ch/user/d/dainesea/code/";
    TString parFile;
    if(!useParFiles) {
      gProof->UploadPackage("AF-v4-17");
      gProof->EnablePackage("AF-v4-17");
      // --- Enable the PWG3vertexingHF Package
      parFile="PWG3vertexingHF.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("PWG3vertexingHF");
    } else {
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
      // --- Enable the PWG3base Package
      parFile="PWG3base.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("PWG3base");
      // --- Enable the PWG3vertexingHF Package
      parFile="PWG3vertexingHF.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("PWG3vertexingHF");
      // --- Enable the PWG3muon Package
      parFile="PWG3muon.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("PWG3muon");
    }
    gProof->ShowEnabledPackages(); // show a list of enabled packages
  }


  // Create Alien plugin, if requested
  if(useAlienPlugin) {  
    if(analysisMode!="grid") {printf("Analysis mode must be grid, to use alien plugin\n"); return;}
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,useParFiles);  
    if(!alienHandler) return;
  }


  //-------------------------------------------------------------------
  // Prepare input
  TChain *chainAOD = 0;
  TString dataset; // for proof

  if(!useAlienPlugin) {
    TString makeAODInputChain="MakeAODInputChain.C"; makeAODInputChain.Prepend(loadMacroPath.Data());
    if(inputMode=="list") {
      // Local files
      gROOT->LoadMacro(makeAODInputChain.Data());
      chainAOD = MakeAODInputChain();// with this it reads ./AliAOD.root and ./AliAOD.VertexingHF.root
      //chainAOD = MakeAODInputChain("alien:///alice/cern.ch/user/r/rbala/newtrain/out_lhc08x/180100/",1,1);
      printf("ENTRIES %d\n",chainAOD->GetEntries());
    } else if(inputMode=="xml") {
      // xml
      gROOT->LoadMacro(makeAODInputChain.Data());
      chainAOD = MakeAODInputChain("collection_aod.xml","collection_aodHF.xml");
    } else if(inputMode=="dataset") {
      // CAF dataset
      //gProof->ShowDataSets();
      dataset="/ITS/dainesea/AODVertexingHF_LHC08x_180100";
    }
  }

  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  // Connect plug-in to the analysis manager
  if(useAlienPlugin) mgr->SetGridHandler(alienHandler);

  // Input
  AliAODInputHandler *inputHandler = new AliAODInputHandler();
  inputHandler->AddFriend("./AliAOD.VertexingHF.root");
  //inputHandler->AddFriend("deltas/AliAOD.VertexingHF.root");
  if(analysisMode=="proof" ) {
    if(saveProofToAlien) mgr->SetSpecialOutputLocation(proofOutdir);
  }
  mgr->SetInputEventHandler(inputHandler);
  //-------------------------------------------------------------------

  
  //-------------------------------------------------------------------
  // Analysis tasks (wagons of the train)   
  //
  TString taskName;
  
  ////// ADD THE FULL D2H TRAIN
  taskName="AddD2HTrain.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  Bool_t readMC=kFALSE;
  AddD2HTrain(readMC);//,1,0,0,0,0,0,0,0,0,0,0);
  
  ////// OR ADD INDIVIDUAL TASKS
  
  /*
    taskName="AddTaskCompareHF.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSECompareHF *cmpTask = AddTaskCompareHF();
    
    taskName="AddTaskD0Mass.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSED0Mass *d0massTask = AddTaskD0Mass();
    AliAnalysisTaskSED0Mass *d0massLikeSignTask = AddTaskD0Mass(1); 
  
    taskName="AddTaskDplus.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEDplus *dplusTask = AddTaskDplus();
  
    taskName="AddTaskDs.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEDs *dsTask = AddTaskDs();
    
    //taskName="AddTaskSelectHF.C"; taskName.Prepend(loadMacroPath.Data());
    //gROOT->LoadMacro(taskName.Data());
    //AliAnalysisTaskSESelectHF *seleTask = AddTaskSelectHF();
    
    taskName="AddTaskBkgLikeSignD0.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEBkgLikeSignD0 *lsD0Task = AddTaskBkgLikeSignD0();
    
    taskName="AddTaskBkgLikeSignJPSI.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEBkgLikeSignJPSI *lsJPSITask = AddTaskBkgLikeSignJPSI();
    
    //taskName="AddTaskBtoJPSItoEle.C"; taskName.Prepend(loadMacroPath.Data());
    //gROOT->LoadMacro(taskName.Data());
    //AliAnalysisTaskSEBtoJPSItoEle *jpsiTask = AddTaskBtoJPSItoEle();
    
    taskName="AddTaskCFMultiVarMultiStep.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliCFHeavyFlavourTaskMultiVarMultiStep *cfmvmsTask = AddTaskCFMultiVarMultiStep();
    
    taskName="AddTaskCharmFraction.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    Int_t switchMC[5]={1,1,1,1,1};
    AliAnalysisTaskSECharmFraction *cFractTask = AddTaskCharmFraction("d0D0",switchMC);
    
    // attach a private task (not committed)
    // (the files MyTask.h MyTask.cxx AddMyTask.C have to be declared in plugin
    // configuration, see below)
    
    if(analysisMode.Data()=="proof") {
    gProof->LoadMacro("MyTask.cxx++g");
    } else {
    gROOT->LoadMacro("MyTask.cxx++g");
    }
    gROOT->LoadMacro("AddMyTask.C");
    MyTask *myTask = AddMyTask();
    
    
    if(analysisMode.Data()=="proof") {
    gProof->LoadMacro("AliDStarJets.cxx++g");
    } else {
    gROOT->LoadMacro("AliDStarJets.cxx++g");
    }
    gROOT->LoadMacro("AddTaskDStarJets.C");
    AliDStarJets *myTask = AddTaskDStarJets();
  */
  //-------------------------------------------------------------------
  
  //
  // Run the analysis
  //    
  if(chainAOD) printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainAOD->GetEntries());
  
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
  if(analysisMode!="proof") {
    mgr->StartAnalysis(analysisMode.Data(),chainAOD,nentries,firstentry);
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
   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(pluginmode.Data());
   plugin->SetUser("rbala");
   plugin->SetNtestFiles(1);
   // Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-26-00b-6");
   plugin->SetAliROOTVersion("v4-19-18-AN");
   // Declare input data to be processed.
   // Method 1: Create automatically XML collections using alien 'find' command.
   // Define production directory LFN
   //  plugin->SetGridDataDir("/alice/cern.ch/user/r/rbala/data_pass4_good_runCINT1B_8thfeb/");
   //plugin->SetGridDataDir("/alice/sim/PDC_09/LHC09a4/AOD3/");
   // Set data search pattern
   plugin->SetGridDataDir("/alice/data/2010/LHC10c");
   plugin->SetDataPattern("pass2/*AliAOD.root");
   // Adds only the good runs from the Monalisa Run Condition Table
   AddGoodRuns(plugin,"LHC10c");
   // ...then add run numbers to be considered
   plugin->SetMaxMergeFiles(100);
   plugin->SetNrunsPerMaster(100);
   plugin->SetNumberOfReplicas(2);
   //  or
   //plugin->SetRunRange(529000,529007);
   // Method 2: Declare existing data files (raw collections, xml collections, root file)
   // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
   // XML collections added via this method can be combined with the first method if
   // the content is compatible (using or not tags)
   //plugin->AddDataFile("/alice/cern.ch/user/r/rbala/newtrain/collection/collection_aod_lhc08w.xml");
   //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
   // Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("myHFanalysis");
   // Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
   // Declare the analysis source files names separated by blancs. To be compiled runtime
   // using ACLiC on the worker nodes.
   //plugin->SetAnalysisSource("AliDStarJets.cxx");
   // Declare all libraries (other than the default ones for the framework. These will be
   // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->SetAdditionalLibs("libPWG3vertexingHF.so libPWG3base.so libPWG3muon.so");
   // use par files
   if(useParFiles) {
     plugin->EnablePackage("STEERBase.par");
     plugin->EnablePackage("ESD.par");
     plugin->EnablePackage("AOD.par");
     plugin->EnablePackage("ANALYSIS.par");
     plugin->EnablePackage("ANALYSISalice.par");
     plugin->EnablePackage("CORRFW.par");
     plugin->EnablePackage("PWG3base.par");
     plugin->EnablePackage("PWG3vertexingHF.par");
     plugin->EnablePackage("PWG3muon.par");
   }
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWG3 -I$ALICE_ROOT/PWG3/vertexingHF -g");
   // Declare the output file names separated by blancs.
   // (can be like: file.root or file.root@ALICE::Niham::File)
   plugin->SetDefaultOutputs(kTRUE);
   //plugin->SetOutputFiles("output.root CmpHF.root CmpHFnt.root D0InvMass.root InvMassDplus.root InvMassDplus_nt1.root InvMassDplus_nt2.root");
   // Optionally define the files to be archived.
   //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
   //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
   // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("AnalysisHF.C");
   // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(10);
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


   return plugin;
}
//----------------------------------------------------------------------------
void AddGoodRuns(AliAnalysisAlien* plugin,TString lhcPeriod) {
  //
  // Adds good runs from the Monalisa Run Condition Table
  //
  plugin->SetRunPrefix("000");

  if(lhcPeriod=="LHC10b") {
    Int_t nruns=52;
    Int_t runlist[52]={117222, 117220, 117120, 117118, 117116, 117112, 117109, 117099, 117092, 117086, 117082, 117077, 117063, 117060, 117059, 117054, 117053, 117052, 117050, 117048, 116684, 116643, 116642, 116611, 116610, 116609, 116574, 116571, 116562, 116561, 116559, 116403, 116402, 116401, 116288, 116287, 116102, 115514, 115414, 115413, 115401, 115393, 115345, 115335, 115328, 115325, 115322, 115318, 115310, 115193, 115186, 114931};
   
    for(Int_t k=0;k<nruns;k++){
      plugin->AddRunNumber(runlist[k]);
    }
    plugin->SetNRuns(nruns);
  }

  if(lhcPeriod=="LHC10c") { 
    Int_t nruns=57;
    Int_t runlist[57]={120829, 120825, 120824, 120823, 120822, 120821, 120820, 120758, 120750, 120741, 120671, 120617, 120616, 120505, 120504, 120503, 120244, 120079, 120076, 120073, 120072, 120069, 120067, 120066, 120064, 119971, 119969, 119965, 119961, 119952, 119948, 119941, 119935, 119934, 119926, 119924, 119923, 119917, 119913, 119909, 119907, 119904, 119862, 119859, 119856, 119853, 119849, 119846, 119845, 119844, 119842, 119841, 119838, 119837, 119163, 119161, 119159};
   
    for(Int_t k=0;k<nruns;k++){
      plugin->AddRunNumber(runlist[k]);
    }
    plugin->SetNRuns(nruns);
  }

  return;
}
