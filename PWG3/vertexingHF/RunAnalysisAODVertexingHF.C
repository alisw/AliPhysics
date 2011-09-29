class AliAnalysisGrid;
class AliAnalysisAlien;

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


  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWG3 -I$ALICE_ROOT/PWG3/base -I$ALICE_ROOT/PWG3/vertexingHF -I$ALICE_ROOT/PWG2/FLOW/AliFlowCommon -I$ALICE_ROOT/PWG2/FLOW/AliFlowTasks -g"); 
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
    TString makeAODInputChain="../MakeAODInputChain.C"; makeAODInputChain.Prepend(loadMacroPath.Data());
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
  if(analysisMode=="proof" ) {
    inputHandler->AddFriend("./AliAOD.VertexingHF.root");
    //inputHandler->AddFriend("deltas/AliAOD.VertexingHF.root");
    if(saveProofToAlien) mgr->SetSpecialOutputLocation(proofOutdir);
  }
  mgr->SetInputEventHandler(inputHandler);
  //-------------------------------------------------------------------

  
  //-------------------------------------------------------------------
  // Analysis tasks (wagons of the train)   
  //
  // First add the task for the PID response setting
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(kFALSE,kTRUE);

  TString taskName;
  
  ////// ADD THE FULL D2H TRAIN
  /*taskName="../AddD2HTrain.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  Bool_t readMC=kFALSE;
  AddD2HTrain(readMC);//,1,0,0,0,0,0,0,0,0,0,0);*/
  
  ////// OR ADD INDIVIDUAL TASKS
  
  
  /*  taskName="AddTaskCompareHF.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSECompareHF *cmpTask = AddTaskCompareHF();
  */ 
    taskName="AddTaskD0Mass.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSED0Mass *d0massTask = AddTaskD0Mass();
    AliAnalysisTaskSED0Mass *d0massLikeSignTask = AddTaskD0Mass(1); 
    /*
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
    
    taskName="AddTaskCFMultiVarMultiStep.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliCFHeavyFlavourTaskMultiVarMultiStep *cfmvmsTask = AddTaskCFMultiVarMultiStep();
    

    taskName="AddTaskSECharmFraction.C"; 
    taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    Int_t switchMC[5]={0,0,0,0,0};
    Int_t ppPbPb=1;// 0 for pp, 1 for PbPb, used to siwtch on/off the removal of daughters from the primary vertex
    AliAnalysisTaskSECharmFraction *cFractTask = AddTaskSECharmFraction("standard",switchMC,readMC,kTRUE,kFALSE,"D0toKpiCharmFractCuts.root","c",ppPbPb);
    // arguments: filename,switchMC,readmc,usepid,likesign,cutfilename,containerprefix

    
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
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(pluginmode.Data());
   plugin->SetUser();
   // Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion();
   plugin->SetAliROOTVersion();

   gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/AddGoodRuns.C");

   // Declare input data to be processed.
   //************************************************
   // Set data search pattern for DATA
   //************************************************   
   plugin->SetGridDataDir("/alice/data/2010/LHC10d"); // specify LHC period
   plugin->SetDataPattern("pass2/AOD018/*AliAOD.root"); // specify reco pass and AOD set
   plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
   // OR plugin->SetFriendChainName("deltas/AliAOD.VertexingHF.root");
   // Adds only the good runs from the Monalisa Run Condition Table
   // More than one period can be added but the period name has to be removed from GridDataDir (to be tested)
   Int_t totruns=0;
   //totruns += AddGoodRuns(plugin,"LHC10b"); // specify LHC period
   //totruns += AddGoodRuns(plugin,"LHC10c"); // specify LHC period
   totruns += AddGoodRuns(plugin,"LHC10d"); // specify LHC period
   plugin->SetNrunsPerMaster(totruns);
   
   //************************************************
   // Set data search pattern for MONTECARLO
   //************************************************
   /* 
   plugin->SetGridDataDir("/alice/sim/LHC10d3"); // specify MC sample
   plugin->SetDataPattern("AOD005/*AliAOD.root"); // specify AOD set
   plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
   // OR plugin->SetFriendChainName("deltas/AliAOD.VertexingHF.root");
   // Adds only the good runs from the Monalisa Run Condition Table 
   // More than one period can be added!
   Int_t totruns=0;
   totruns += AddGoodRuns(plugin,"LHC10b","LHC10d3"); // specify LHC period for anchor runs; and the name of the MC production
   //totruns += AddGoodRuns(plugin,"LHC10c","LHC10f7"); // specify LHC period for anchor runs;  and the name of the MC production
   //totruns += AddGoodRuns(plugin,"LHC10d","LHC10f7"); // specify LHC period for anchor runs;  and the name of the MC production
   plugin->SetNrunsPerMaster(totruns);
   */
   //
   // Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("myHFanalysis");
   // Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
   // Declare the analysis source files names separated by blancs. To be compiled runtime
   // using ACLiC on the worker nodes.
   //plugin->SetAnalysisSource("AliDStarJets.cxx");
   // Declare all libraries (other than the default ones for the framework. These will be
   // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->SetAdditionalLibs("libPWG3base.so libPWG3vertexingHF.so");
   // use par files
   if(useParFiles) {
     plugin->EnablePackage("STEERBase.par");
     plugin->EnablePackage("ESD.par");
     plugin->EnablePackage("AOD.par");
     plugin->EnablePackage("ANALYSIS.par");
     plugin->EnablePackage("OADB.par");
     plugin->EnablePackage("ANALYSISalice.par");
     plugin->EnablePackage("CORRFW.par");
     plugin->EnablePackage("PWG3base.par");
     plugin->EnablePackage("PWG3vertexingHF.par");
   }
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWG3 -I$ALICE_ROOT/PWG3/base -I$ALICE_ROOT/PWG3/vertexingHF -g");

   plugin->SetDefaultOutputs(kTRUE);
   // merging via jdl
   plugin->SetMergeViaJDL(kTRUE);
   plugin->SetOneStageMerging(kFALSE);
   plugin->SetMaxMergeStages(2);

   // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("AnalysisHF.C");
   // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   // Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("TaskHF.jdl");

   return plugin;
}
