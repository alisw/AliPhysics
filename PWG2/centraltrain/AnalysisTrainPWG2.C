//===================== ANALYSIS TRAIN ==============================
// To use: copy this macro to your work directory, modify the global
// part to match your needs, then run root.
// 
//    root[0] .L AnalysisTrain.C
// Grid full mode as below (other modes: test, offline, submit, terminate)
//    root[1] AnalysisTrainNew("grid", "full")
// CAF mode (requires root v5-23-02 + aliroot v4-16-Rev08)
//    root[2] AnalysisTrainNew("proof")
// Local mode requires AliESds.root or AliAOD.root in ./data directory
//    root[3] AnalysisTrainNew("local")
// 
// In proof and grid modes, a token is needed and sourcing the
// produced environment file.
//
// If 'saveTrain' flag is set, the train will generate a directory
// name and run in this directory. A configuration file
// 'ConfigTrain.C' will be generated.  One can replay at any time the
// train via:
// 
//    root[1] AnalysisTrainNew(ana_mode, plugin_mode, "train_default_<date>/ConfigTrain.C")

//==================   TRAIN NAME   ==================================
TString     train_name         = "PWG2";    // local folder name
TString     train_tag          = "_Pb-Pb_"; // Train special tag
					    // appended to visible
					    // name. ("data", "sim",
					    // "pp", "highmult", ...) 
// Name in train page (DON'T CHANGE)
TString     visible_name       = Form("PWG2%s$2_$3", train_tag.Data()); //# FIXED #
// Add train composition and other comments
TString     job_comment        = "PWG2 tasks on ESDs";
TString     job_tag            = Form("%s: %s", visible_name.Data(), 
				      job_comment.Data());
//====================================================================

// ### Settings that make sense in PROOF only
//====================================================================
TString     proof_cluster      = "alice-caf.cern.ch";
Bool_t      useAFPAR           = kFALSE;  // use AF special par file
TString     AFversion          = "AF-v4-17";
// Change CAF dataset here
TString     proof_dataset      = "/COMMON/COMMON/LHC09a4_run8100X#/esdTree";
TString     proof_outdir       = "";

// ### Settings that make sense when using the Alien plugin
//====================================================================
Int_t       runOnData          = 1;       // Set to 1 if processing real data
Int_t       iCollision         = 1;       // 0=pp, 1=Pb-Pb
Bool_t      usePLUGIN          = kTRUE;   // do not change
Bool_t      useProductionMode  = kTRUE;   // use the plugin in production mode

// Usage of par files ONLY in grid mode and ONLY if the code is not
// available in the deployed AliRoot versions. Par file search path:
// local dir, if not there $ALICE_ROOT.  To refresh par files, remove
// the ones in the workdir, then do "make <target.par>" in AliRoot.
Bool_t      usePAR             = kFALSE;  // use par files for extra libs
Bool_t      useCPAR            = kFALSE;  // use par files for common libs
TString     root_version       = "v5-27-06c";  // *CHANGE ME IF MORE
					       // *RECENT IN GRID* 
TString     aliroot_version    = "v4-21-15-AN";  // *CHANGE ME IF MORE
						 // *RECENT IN GRID* 
// Change production base directory here (test mode)
TString     alien_datadir      = "/alice/data/2010/LHC10h";
// Work directory in GRID (DON'T CHANGE)
TString     grid_workdir       = "/alice/cern.ch/user/a/alidaq/PWG2/PWG2_$2";
// Data pattern - change as needed for test mode
// TString     data_pattern       = "*ESDs/pass1/AOD033/*AOD.root";
TString     data_pattern       = "*ESDs/pass1/*ESDs.root";
// Set the run range
Int_t run_numbers[10] = {137844}; // **********************!!!!!!!
// AliEn output directory. If blank will become output_<train_name>
// Output directory (DON'T CHANGE)
TString     alien_outdir       = "$1/PWG2OUT$2";
// Input collection (production mode)
TString     data_collection    = "$1/input.xml";
// Output folder to write delta AOD's. Considered if not null.
TString     outputSingleFolder = "";
//TString     outputSingleFolder = "deltas";
// Number of files merged in a chunk
Int_t       maxMergeFiles      = 20;
// Files that should not be merged
TString     mergeExclude       = ("AliAOD.root "
				  "AliAOD.VertexingHF.root "
				  "AliAOD.Jets.root "
				  "deltaAODPartCorr.root "
				  "AliAOD.Muons.root "
				  "AliAOD.Dimuons.root "
				  "AliAOD.Dielectron.root "
				  "AliAODCentrality.root");
TString     mergeDirName       = "PWG2OUT$2";
// Make replicas on the storages below
TString     outputStorages      = "disk=4";
// Number of runs per master job
Int_t       nRunsPerMaster     = 10;
// Maximum number of files per job (gives size of AOD)
Int_t       nFilesPerJob       = 5;
// Int_t       nFilesPerJob       = 1; (AOD->delta AOD production case)
// ### Settings that make sense only for local analysis
//==============================================================================
// Change local xml dataset for local interactive analysis
TString     local_xmldataset   = "";

// ### Other flags to steer the analysis
//==============================================================================
Bool_t      usePhysicsSelection = kTRUE;  // use physics selection
Bool_t      useBKrejection      = kFALSE; // use BK rejection
Bool_t      useCentrality       = kTRUE;  // centrality delta AOD
Bool_t      useTender           = kFALSE; // use tender wagon
Bool_t      useV0tender         = kFALSE; // use V0 correction in tender
Bool_t      useMergeViaJDL      = kTRUE;  // merge via JDL
Bool_t      useFastReadOption   = kFALSE; // use xrootd tweaks
Bool_t      useOverwriteMode    = kTRUE;  // overwrite existing collections
Bool_t      useDATE             = kFALSE; // use date in train name
Bool_t      useDBG              = kTRUE;  // activate debugging
Bool_t      useMC               = kFALSE; // use MC info
Bool_t      useTAGS             = kFALSE; // use ESD tags for selection
Bool_t      useKFILTER          = kFALSE; // use Kinematics filter
Bool_t      useTR               = kFALSE; // use track references
Bool_t      useCORRFW           = kTRUE;  // do not change
Bool_t      useAODTAGS          = kFALSE; // use AOD tags
Bool_t      saveTrain           = kTRUE;  // save train configuration as:
Bool_t      saveCanvases        = kFALSE; // save canvases created in Terminate
Bool_t      saveProofToAlien    = kFALSE; // save proof outputs in AliEn

// ### Analysis modules to be included. Some may not be yet fully implemented.
//==============================================================================
Int_t       iAODanalysis       = 0;      // Analysis on input AOD's
Int_t       iAODhandler        = 1;      // Analysis produces an AOD or dAOD's

Int_t       iPWG2fmd           = 1;      // FMD forward analysis (PWG2)
Int_t       iPWG2femto         = 1;      // Femtoscopy two-track analysis (PWG2)
Int_t       iPWG2spectra       = 1;      // Charge hadron spectra (PWG2)

// Temporaries.
TString anaPars = "";
TString anaLibs = "";
// Function signatures
class AliAnalysisAlien;

//______________________________________________________________________________
void AnalysisTrainPWG2(const char *analysis_mode="local",
		       const char *plugin_mode="full",
		       const char *config_file="")
{
  // Main analysis train macro. If a configuration file is provided,
  // all parameters are taken from there but may be altered by
  // CheckModuleFlags.
  if (strlen(config_file) && !LoadConfig(config_file)) return;
  TString smode(analysis_mode);
  smode.ToUpper();
  TString spmode(plugin_mode);
  spmode.ToLower();
  if (spmode == "test") useProductionMode = kFALSE;
  // Check compatibility of selected modules
  CheckModuleFlags(smode);
  if (saveTrain)              WriteConfig();

  printf("=================================================================\n");
  printf("===========    RUNNING ANALYSIS TRAIN %s IN %s MODE   ==========\n", 
	 train_name.Data(),smode.Data());
  printf("=================================================================\n");
  printf("=  Configuring analysis train for:                              =\n");
  if (iAODanalysis)        printf("=  %-60s  =\n", "AOD analysis");
  else                     printf("=  %-60s  =\n", "ESD analysis");
  if (usePhysicsSelection) printf("=  %-60s  =\n", "Physics selection");
  if (useTender)           printf("=  %-60s  =\n", "TENDER");
  if (iPWG2fmd)            printf("=  %-60s  =\n", "PWG2 FMD");
  if (iPWG2femto)          printf("=  %-60s  =\n", "PWG2 Femtoscopy");
  if (iPWG2spectra)        printf("=  %-60s  =\n", "PWG2 Charged Spectra");
  printf("=================================================================\n");
  printf(":: use physics selection: %d\n", (UInt_t)usePhysicsSelection);
  printf(":: use xrootd tweaks:     %d\n", (UInt_t)useFastReadOption);
  printf(":: use overwrite xml    : %d\n", (UInt_t)useOverwriteMode);
  printf(":: use merge via JDL:     %d\n", (UInt_t)useMergeViaJDL);
  printf(":: use MC truth:          %d\n", (UInt_t)useMC);
  printf(":: use KINE filter:       %d\n", (UInt_t)useKFILTER);
  printf(":: use track references:  %d\n", (UInt_t)useTR);
  printf(":: use tags:              %d\n", (UInt_t)useTAGS);
  printf(":: use AOD tags:          %d\n", (UInt_t)useAODTAGS);
  printf(":: use debugging:         %d\n", (UInt_t)useDBG);
  printf(":: use PAR files:         %d\n", (UInt_t)usePAR);
  printf(":: use AliEn plugin:      %d\n", (UInt_t)usePLUGIN);

  //==================================================================
  // Connect to back-end system
  // if (!Connect(smode)) {
  //   ::Error("AnalysisTrain", "Could not connect to %s back-end", 
  //           analysis_mode);
  //   return;
  // }

  // Load common libraries and set include path
  if (!LoadCommonLibraries(smode)) {
    ::Error("AnalysisTrain", "Could not load common libraries");
    return;
  }

  // Make the analysis manager and connect event handlers
  AliAnalysisManager *mgr  = 
    new AliAnalysisManager("Analysis Train", "Production train");
  if (saveProofToAlien) mgr->SetSpecialOutputLocation(proof_outdir);
  if (!strcmp(plugin_mode, "test")) mgr->SetNSysInfo(1);
  // Load analysis specific libraries
  if (!LoadAnalysisLibraries(smode)) {
    ::Error("AnalysisTrain", "Could not load analysis libraries");
    return;
  }

  // Create input handler (input container created automatically)
  if (iAODanalysis) {
    // AOD input handler
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
  } else {
    // ESD input handler
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    if (useTAGS) esdHandler->SetReadTags();
    mgr->SetInputEventHandler(esdHandler);
  }
  // Monte Carlo handler
  if (useMC && !iAODanalysis) {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
    mcHandler->SetReadTR(useTR);
  }
  // AOD output container, created automatically when setting an AOD handler
  if (iAODhandler) {
    // AOD output handler
    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName("AliAOD.root");
    mgr->SetOutputEventHandler(aodHandler);
    if (iAODanalysis) {
      aodHandler->SetFillAOD(kFALSE);
      aodHandler->SetCreateNonStandardAOD();
      // if (iPWG3vertexing) 
      //   aodHandler->SetOutputFileName("AliAOD.VertexingHF.root");
    }
  }
  // Debugging if needed
  if (useDBG) mgr->SetDebugLevel(3);
  if (saveCanvases) mgr->SetSaveCanvases(kTRUE);

  //==========================================================================
  // Create the chain. In this example it is created only from ALIEN
  // files but can be done to work in batch or grid mode as well.
  TChain *chain = CreateChain(smode, plugin_mode);

  //==========================================================================
  // Load the tasks configuration macros for all wagons. These files
  // are supposed now to be in the current workdir, but in AliEn they
  // will be in the file catalog, mapped from AliRoot and pecified in
  // the jdl input list.
  // 
  // For now connection to top input container and common AOD output
  // container is done in this macro, but in future these containers
  // will be connected from each task configuration macro.
  AddAnalysisTasks();

  // Run the analysis
  //
  if (usePLUGIN) {
    AliAnalysisGrid *alienHandler = CreateAlienHandler(plugin_mode);
    AliAnalysisManager::GetAnalysisManager()->SetGridHandler(alienHandler);
  }

  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    if (saveTrain || strlen(config_file)) gSystem->ChangeDirectory(train_name);
    StartAnalysis(smode, chain);
  }
}

//____________________________________________________________________
void AddToMacroPath(const char* path)
{
  TString macroPath(gROOT->GetMacroPath());
  TString tgt(gSystem->ExpandPathName(path));
  if (macroPath.Contains(tgt)) return;
  macroPath.Append(":");
  macroPath.Append(tgt);
  gROOT->SetMacroPath(macroPath.Data());
}

//____________________________________________________________________
void AddAnalysisTasks()
{
  // Add all analysis task wagons to the train
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  //
  // Tender and supplies. Needs to be called for every event.
  //
  if (useTender) {
    AddToMacroPath("$ALICE_ROOT/ANALYSIS/TenderSupplies/");
    gROOT->LoadMacro("AddTaskTender.C");
    // IF V0 tender needed, put kTRUE below
    AliAnalysisTaskSE *tender = AddTaskTender(useV0tender);
    // tender->SetDebugLevel(2);
  }
  
  if (usePhysicsSelection) {
    // Physics selection task
    AddToMacroPath("$ALICE_ROOT/ANALYSIS/macros/");
    gROOT->LoadMacro("AddTaskPhysicsSelection.C");
    mgr->RegisterExtraFile("event_stat.root");
    AliPhysicsSelectionTask *physSelTask =
      AddTaskPhysicsSelection(useMC,useBKrejection);
    mgr->AddStatisticsTask(AliVEvent::kMB);
  }
  
  if (useCentrality) {
    // Common Centrality task
    AddToMacroPath("$ALICE_ROOT/ANALYSIS/macros/");
    gROOT->LoadMacro("AddTaskCentrality.C");
    AliCentralitySelectionTask* ctask = AddTaskCentrality();
    if (ctask) {
      Int_t pass = 0;
      if      (data_pattern.Contains("pass1")) pass = 1;
      else if (data_pattern.Contains("pass2")) pass = 2;
      else if (data_pattern.Contains("pass3")) pass = 3;
      ctask->SetPass(pass);
      if (useMC) ctask->SetMCInput();
    }
    else
      ::Warning("AnalysisTrainNew", "AliCentralitySelectionTask cannot "
		"run for this train conditions - EXCLUDED");
  }
  
  // AOD tags
  if (useAODTAGS) {
    AliAnalysisTaskTagCreator* tagTask =
      new AliAnalysisTaskTagCreator("AOD Tag Creator");
    mgr->AddTask(tagTask);
    AliAnalysisDataContainer *coutTags =
      mgr->CreateContainer("cTag",  TTree::Class(),
			   AliAnalysisManager::kOutputContainer,
			   "AOD.tag.root");
    mgr->ConnectInput (tagTask, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(tagTask, 1, coutTags);
  }


  // ********** PWG2 wagons *****************************************
  AliAnalysisManager::SetCommonFileName("PWG2histograms.root");
  
  // PWG2 FMD
  if (iPWG2fmd) {
    AddToMacroPath("$ALICE_ROOT/PWG2/FORWARD/analysis2/");
    
    gROOT->LoadMacro("AddTaskForwardMult.C");
    gROOT->LoadMacro("AddTaskCentralMult.C");
    /* Input parameters for the analysis task.
     *
     * - sys, collision system 0: get from data, 1: pp, 2: AA
     * - sNN, Center of mass energy (per nucleon pair) 0: get from
     *        data, otherwise sqrt(sNN) in GeV
     * - fld, L3 magnetic field. if sNN=0: get from data, otherwise
     *        magnetic field in kG
     */
    UShort_t sys = iCollision+1;
    UShort_t sNN = 0;
    Short_t  fld = 0.5;
    AliAnalysisTask *taskForward = AddTaskForwardMult(useMC, sys, sNN, fld);
    if (!taskForward)
      ::Warning("AnalysisTrainNew", "AliForwardMultiplicityTask cannot "
		"run for this train conditions - EXCLUDED");
    AliAnalysisTask *taskCentral = AddTaskCentralMult(useMC, sys, sNN, fld);
    if (!taskCentral)
      ::Warning("AnalysisTrainNew", "AliCentralMultiplicityTask cannot "
		"run for this train conditions - EXCLUDED");
  }
  
  // PWG2 Femtoscopy
  if (iPWG2femto) {
    AddToMacroPath("$ALICE_ROOT/PWG2/FEMTOSCOPY/macros");
    // AddToMacroPath("$ALICE_ROOT/PWG2/FEMTOSCOPY/macros/Train/TwoTrackQA");
    gROOT->LoadMacro("AddTaskFemto.C");
    AliAnalysisTask *taskFemto = 
      AddTaskFemto("$ALICE_ROOT/PWG2/FEMTOSCOPY/macros/Train/TwoTrackQA/ConfigFemtoAnalysisCentral.C",
		   Form("%i, %i, %i", iAODanalysis?0:1, iCollision, 0));
    if (!taskFemto) 
      ::Warning("AnalysisTrainNew", "AliAnalysisTaskFemto cannot run for "
		"this train conditions - EXCLUDED");
  }

  // PWG2 Spectra
  if (iPWG2spectra) {
    AddToMacroPath("$ALICE_ROOT/PWG2/SPECTRA/macros");
    gROOT->LoadMacro("AddTaskITSsaSpectra.C");
    AliAnalysisTask *taskSpectraITSsa = 
      AddTaskITSsaSpectra(0,useMC,-1,-1,iCollision);
    if (!taskSpectraITSsa) 
      ::Warning("AnalysisTrainNew", "AliAnalysisTaskSpectraITSsa cannot "
		"run for this train conditions - EXCLUDED");

    gROOT->LoadMacro("AddTaskChargedHadronSpectraITSTruncatedMean.C");
    AliAnalysisTask *taskSpectraITSTPC = 
      AddTaskChargedHadronSpectraITSTruncatedMean(-1,-1,useMC,iCollision,
						  "configChargedHadronSpectraITSTruncatedMeanTask.C");
    if (!taskSpectraITSTPC) 
      ::Warning("AnalysisTrainNew", "AliAnalysisTaskSpectraITSTPC "
		"cannot run for this train conditions - EXCLUDED");
  }
}

//______________________________________________________________________________
void StartAnalysis(const char *mode, TChain *chain) {
  // Start analysis.
  Int_t imode = -1;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!strcmp(mode, "LOCAL")) imode = 0;
  if (!strcmp(mode, "PROOF")) imode = 1;
  if (!strcmp(mode, "GRID"))  imode = 2;
  switch (imode) {
  case 0:
    if (!chain) {
      ::Error("AnalysisTrainNew.C::StartAnalysis", "Cannot create the chain");
      return;
    }
    mgr->StartAnalysis(mode, chain);
    return;
  case 1:
    if (!proof_dataset.Length()) {
      ::Error("AnalysisTrainNew.C::StartAnalysis", "proof_dataset is empty");
      return;
    }
    mgr->StartAnalysis(mode, proof_dataset, 1000);
    return;
  case 2:
    if (usePLUGIN) {
      if (!mgr->GetGridHandler()) {
	::Error("AnalysisTrainNew.C::StartAnalysis", 
		"Grid plugin not initialized");
	return;
      }
      mgr->StartAnalysis("grid");
    } else {
      if (!chain) {
	::Error("AnalysisTrainNew.C::StartAnalysis", "Cannot create the chain");
	return;
      }
      mgr->StartAnalysis(mode, chain);
    }
    return;
  }
}

//______________________________________________________________________________
void CheckModuleFlags(const char *mode) 
{
  // Checks selected modules and insure compatibility
  Int_t imode = -1;
  if (!strcmp(mode, "LOCAL")) imode = 0;
  if (!strcmp(mode, "PROOF")) imode = 1;
  if (!strcmp(mode, "GRID"))  imode = 2;
  if (imode==1) {
    if (!usePAR) {
      ::Info("AnalysisTrainNew.C::CheckModuleFlags", 
	     "PAR files enabled due to PROOF analysis");
      usePAR = kTRUE;
    }
  }
  if (imode != 2) {
    ::Info("AnalysisTrainNew.C::CheckModuleFlags", 
	   "AliEn plugin disabled since not in GRID mode");
    usePLUGIN = kFALSE;
  }
  if (iAODanalysis) {
    // AOD analysis
    if (useMC)
      ::Info("AnalysisTrainNew.C::CheckModuleFlags", 
	     "MC usage disabled in analysis on AOD's");
    if (useAODTAGS)
      ::Info("AnalysisTrainNew.C::CheckModuleFlags", 
	     "AOD tags usage disabled in analysis on AOD's");
    useMC = kFALSE;
    useTR = kFALSE;
    useAODTAGS = kFALSE;
  } else {
    // ESD analysis
    if (!useMC) useTR = kFALSE;
    if (!useTR) {
      ::Info("AnalysisTrainNew.C::CheckModuleFlags", 
	     "iPWG2evchar disabled if not reading track references");
      iPWG2evchar = 0;
    }
  }
  if (useKFILTER && !useMC) useKFILTER = kFALSE;
  if (useAODTAGS && !iAODhandler) useAODTAGS = kFALSE;
}

//______________________________________________________________________________
Bool_t Connect(const char *mode) 
{
  // Connect <username> to the back-end system.
  Int_t imode = -1;
  if (!strcmp(mode, "LOCAL")) imode = 0;
  if (!strcmp(mode, "PROOF")) imode = 1;
  if (!strcmp(mode, "GRID"))  imode = 2;
  TString username = gSystem->Getenv("alien_API_USER");
  switch (imode) {
  case 0:
    break;
  case 1:
    if  (!username.Length()) {
      ::Error(Form("AnalysisTrainNew.C::Connect <%s>", mode), 
	      "Make sure you:\n "
	      "\t1. Have called: alien-token-init <username>\n"
	      "\t2. Have called: >source /tmp/gclient_env_$UID");
      return kFALSE;
    }
    ::Info("AnalysisTrainNew.C::Connect", 
	   "Connecting user <%s> to PROOF cluster <%s>",
	   username.Data(), proof_cluster.Data());
    gEnv->SetValue("XSec.GSI.DelegProxy", "2");
    //  TProof::Open(Form("%s@%s:31093", username.Data(), 
    //               proof_cluster.Data()));
    TProof::Open(Form("%s@%s", username.Data(), proof_cluster.Data()));
    if (!gProof) {
      if (strcmp(gSystem->Getenv("XrdSecGSISRVNAMES"), "lxfsrd0506.cern.ch"))
	::Error(Form("AnalysisTrainNew.C::Connect <%s>", mode), 
		"Environment XrdSecGSISRVNAMES different from "
		"lxfsrd0506.cern.ch");
      return kFALSE;
    }
    TGrid::Connect("alien://");
    if (gGrid) {
      TString homedir = gGrid->GetHomeDirectory();
      TString workdir = homedir + train_name;
      if (!gGrid->Cd(workdir)) {
	gGrid->Cd(homedir);
	if (gGrid->Mkdir(workdir)) {
	  gGrid->Cd(train_name);
	  ::Info("AnalysisTrainNew::Connect()", "Directory %s created", 
		 gGrid->Pwd());
	}
      }
      gGrid->Mkdir("proof_output");
      gGrid->Cd("proof_output");
      proof_outdir = Form("alien://%s", gGrid->Pwd());
    }
    break;
  case 2:
    if (usePLUGIN && !gSystem->Getenv("alien_CLOSE_SE")) {
      ::Error(Form("AnalysisTrainNew.C::Connect <%s>", mode),
	      "When using the AliEn plugin it is preferable to define the "
	      "variable alien_CLOSE_SE in your environment.");
      return kFALSE;
    }
    ::Info("AnalysisTrainNew.C::Connect", "Connecting user <%s> to AliEn ...",
	   username.Data());
    TGrid::Connect("alien://");
    if (!gGrid || !gGrid->IsConnected()) return kFALSE;
    break;
  default:
    ::Error("AnalysisTrainNew.C::Connect", "Unknown run mode: %s", mode);
    return kFALSE;
  }
  ::Info("AnalysisTrainNew.C::Connect","Connected in %s mode", mode);
  return kTRUE;
}

//______________________________________________________________________________
Bool_t LoadCommonLibraries(const char *mode)
{
  // Load common analysis libraries.
  Int_t imode = -1;
  if (!strcmp(mode, "LOCAL")) imode = 0;
  if (!strcmp(mode, "PROOF")) imode = 1;
  if (!strcmp(mode, "GRID"))  imode = 2;
  if (!gSystem->Getenv("ALICE_ROOT")) {
    ::Error("AnalysisTrainNew.C::LoadCommonLibraries", 
	    "Analysis train requires that analysis libraries are "
	    "compiled with a local AliRoot");
    return kFALSE;
  }
  Bool_t success = kTRUE;
  // ROOT libraries
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so");

  // Load framework classes. Par option ignored here.
  switch (imode) {
  case 0:
  case 2:
    if (useCPAR) {
      success &= LoadLibrary("STEERBase", mode, kTRUE);
      success &= LoadLibrary("ESD", mode, kTRUE);
      success &= LoadLibrary("AOD", mode, kTRUE);
      success &= LoadLibrary("ANALYSIS", mode, kTRUE);
      success &= LoadLibrary("ANALYSISalice", mode, kTRUE);
      if (useCORRFW) success &= LoadLibrary("CORRFW", mode, kTRUE);
    } else {
      success &= LoadLibrary("libSTEERBase.so", mode);
      success &= LoadLibrary("libESD.so", mode);
      success &= LoadLibrary("libAOD.so", mode);
      success &= LoadLibrary("libANALYSIS.so", mode);
      success &= LoadLibrary("libANALYSISalice.so", mode);
      if (useCORRFW) success &= LoadLibrary("libCORRFW.so", mode);
      gROOT->ProcessLine(".include $ALICE_ROOT/include");
    }
    break;
  case 1:
    Int_t ires = -1;
    if (useAFPAR && !gSystem->AccessPathName(AFversion)) 
      ires = gProof->UploadPackage(AFversion);
    if (ires < 0) {
      success &= LoadLibrary("STEERBase", mode);
      success &= LoadLibrary("ESD", mode);
      success &= LoadLibrary("AOD", mode);
      success &= LoadLibrary("ANALYSIS", mode);
      success &= LoadLibrary("ANALYSISalice", mode);
      if (useCORRFW) success &= LoadLibrary("CORRFW", mode);
    } else {
      ires = gProof->EnablePackage(AFversion);
      if (ires<0) success = kFALSE;
      if (useCORRFW) success &= LoadLibrary("CORRFW", mode);
    }
    break;
  default:
    ::Error("AnalysisTrainNew.C::LoadCommonLibraries", 
	    "Unknown run mode: %s", mode);
    return kFALSE;
  }
  if (success) {
    ::Info("AnalysisTrainNew.C::LoadCommodLibraries", 
	   "Load common libraries:    SUCCESS");
    ::Info("AnalysisTrainNew.C::LoadCommodLibraries", 
	   "Include path for Aclic compilation:\n%s",
	   gSystem->GetIncludePath());
  } else {
    ::Info("AnalysisTrainNew.C::LoadCommodLibraries", 
	   "Load common libraries:    FAILED");
  }

  return success;
}

//______________________________________________________________________________
Bool_t LoadAnalysisLibraries(const char *mode)
{
  // Load common analysis libraries.
  Bool_t success = kTRUE;
  if (useTender) {
    if (!LoadLibrary("TENDER", mode, kTRUE) ||
	!LoadLibrary("TENDERSupplies", mode, kTRUE)) return kFALSE;
  }
  // PWG2 fmd
  if (iPWG2fmd) {
    if (!LoadLibrary("PWG2forward2", mode, kTRUE)) return kFALSE;
  }
  // PWG2 femtoscopy
  if (iPWG2femto) {
    if (!LoadLibrary("PWG2AOD", mode, kTRUE) ||
	!LoadLibrary("PWG2femtoscopy", mode, kTRUE) ||
	!LoadLibrary("PWG2femtoscopyUser", mode, kTRUE)) return kFALSE;
  }
  // PWG2 spectra
  if (iPWG2spectra) {
    if (!LoadLibrary("PWG2spectra", mode, kTRUE)) return kFALSE;
  }
  ::Info("AnalysisTrainNew.C::LoadAnalysisLibraries", 
	 "Load other libraries:   SUCCESS");
  return kTRUE;
}

//______________________________________________________________________________
Bool_t LoadLibrary(const char *module, const char *mode, Bool_t rec=kFALSE)
{
  // Load a module library in a given mode. Reports success.
  Int_t imode = -1;
  Int_t result;
  TString smodule(module);
  if (!strcmp(mode, "LOCAL")) imode = 0;
  if (!strcmp(mode, "PROOF")) imode = 1;
  if (!strcmp(mode, "GRID"))  imode = 2;
  TString mod(module);
  if (!mod.Length()) {
    ::Error("AnalysisTrainNew.C::LoadLibrary", "Empty module name");
    return kFALSE;
  }
  // If a library is specified, just load it
  if (smodule.EndsWith(".so")) {
    mod.Remove(mod.Index(".so"));
    result = gSystem->Load(mod);
    if (result < 0) {
      ::Error("AnalysisTrainNew.C::LoadLibrary", 
	      "Could not load library %s", module);
      return kFALSE;
    }
    if (rec) anaLibs += Form("%s.so ",mod.Data());
    return kTRUE;
  }
  // Check if the library is already loaded
  if (strlen(gSystem->GetLibraries(Form("%s.so", module), "", kFALSE)) > 0)
    return kTRUE;
  switch (imode) {
  case 0:
  case 2:
    if (usePAR) {
      result = SetupPar(module);
      if (rec) anaPars += Form("%s.par ", module);
    } else {
      result = gSystem->Load(Form("lib%s.so", module));
      if (rec) anaLibs += Form("lib%s.so ", module);
    }
    break;
  case 1:
    result = gProof->UploadPackage(module);
    if (result<0) {
      result = 
	gProof->UploadPackage(gSystem->ExpandPathName(Form("$ALICE_ROOT/%s.par",
							   module)));
      if (result<0) {
	::Error("AnalysisTrainNew.C::LoadLibrary", 
		"Could not find module %s.par in current directory "
		"nor in $ALICE_ROOT", module);
	return kFALSE;
      }
    }
    result = gProof->EnablePackage(module);
    break;
  default:
    return kFALSE;
  }
  if (result < 0) {
    ::Error("AnalysisTrainNew.C::LoadLibrary", 
	    "Could not load module %s", module);
    return kFALSE;
  }
  return kTRUE;
}


//______________________________________________________________________________
TChain *CreateChain(const char *mode, const char *plugin_mode)
{
  // Create the input chain
  Int_t imode = -1;
  if (!strcmp(mode, "LOCAL")) imode = 0;
  if (!strcmp(mode, "PROOF")) imode = 1;
  if (!strcmp(mode, "GRID"))  imode = 2;
  TChain *chain = NULL;
  // Local chain
  switch (imode) {
  case 0:
    if (iAODanalysis) {
      if (!local_xmldataset.Length()) {
	// Local AOD
	chain = new TChain("aodTree");
	if (gSystem->AccessPathName("data/AliAOD.root"))
	  ::Error("AnalysisTrainNew.C::CreateChain", 
		  "File: AliAOD.root not in ./data dir");
	else {
	  if (!saveTrain) chain->Add("data/AliAOD.root");
	  else            chain->Add("../data/AliAOD.root");
	}
      } else {
	// Interactive AOD
	chain = CreateChainSingle(local_xmldataset, "aodTree");
      }
    } else {
      if (!local_xmldataset.Length()) {
	// Local ESD
	chain = new TChain("esdTree");
	if (gSystem->AccessPathName("data/AliESDs.root"))
	  ::Error("AnalysisTrainNew.C::CreateChain", 
		  "File: AliESDs.root not in ./data dir");
	else {
	  if (!saveTrain) chain->Add("data/AliESDs.root");
	  else            chain->Add("../data/AliESDs.root");
	}
      } else {
	// Interactive ESD
	chain = CreateChainSingle(local_xmldataset, "esdTree");
      }
    }
    break;
  case 1:
    break;
  case 2:
    if (usePLUGIN) {
      // AliAnalysisGrid *alienHandler = CreateAlienHandler(plugin_mode);
      // AliAnalysisManager::GetAnalysisManager()->SetGridHandler(alienHandler);
    } else {
      TString           treeName = "esdTree";
      if (iAODanalysis) treeName = "aodTree";
      chain = CreateChainSingle("wn.xml", treeName);
    }
    break;
  default:
  }
  if (chain && chain->GetNtrees()) return chain;
  return NULL;
}

//______________________________________________________________________________
TChain* CreateChainSingle(const char* xmlfile, const char *treeName)
{
  printf("*******************************\n");
  printf("*** Getting the ESD Chain   ***\n");
  printf("*******************************\n");
  TAlienCollection * myCollection  = TAlienCollection::Open(xmlfile);

  if (!myCollection) {
    ::Error("AnalysisTrainNew.C::CreateChainSingle", 
	    "Cannot create an AliEn collection from %s", xmlfile) ;
    return NULL ;
  }

  TChain* chain = new TChain(treeName);
  myCollection->Reset() ;
  while ( myCollection->Next() ) chain->Add(myCollection->GetTURL("")) ;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
Int_t SetupPar(char* pararchivename)
{
  if (!pararchivename || !strlen(pararchivename)) return -1;
  char processline[1024];
  if (gSystem->AccessPathName(Form("%s.par", pararchivename))) {
    if (!gSystem->AccessPathName(Form("%s/%s.par", 
				      gSystem->Getenv("ALICE_ROOT"),
				      pararchivename))) {
      ::Info("AnalysisTrainNew.C::SetupPar", 
	     "Getting %s.par from $ALICE_ROOT", pararchivename);
      TFile::Cp(gSystem->ExpandPathName(Form("$ALICE_ROOT/%s.par", 
					     pararchivename)),
		Form("%s.par",pararchivename));
    } else {
      ::Error("AnalysisTrainNew.C::SetupPar", "Cannot find %s.par", 
	      pararchivename);
      return -1;
    }
  }
  if (usePLUGIN && saveTrain) 
    gSystem->Exec(Form("ln -s ../%s.par %s",pararchivename, train_name.Data()));
  gSystem->Exec(Form("tar xvzf %s.par", pararchivename));

  TString ocwd = gSystem->WorkingDirectory();
  if (!gSystem->ChangeDirectory(pararchivename)) return -1;

  // check for BUILD.sh and execute
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    printf("*******************************\n");
    if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
      Error("runProcess","Cannot Build the PAR Archive! - Abort!");
      return -1;
    }
  }

  // check for SETUP.C and execute
  if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
    printf("*******************************\n");
    printf("*** Setup PAR archive       ***\n");
    printf("*******************************\n");
    gROOT->Macro("PROOF-INF/SETUP.C");
  }
  if (!gSystem->ChangeDirectory(ocwd.Data())) return -1;
  return 0;
}

//______________________________________________________________________________
AliAnalysisAlien* CreateAlienHandler(const char *plugin_mode)
{
  // Check if user has a valid token, otherwise make one. This has
  // limitations.  One can always follow the standard procedure of
  // calling alien-token-init then
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  // Set the run mode (can be "full", "test", "offline", "submit" or
  // "terminate")
  plugin->SetRunMode(plugin_mode);
  if (useProductionMode) {
    plugin->SetProductionMode();
    plugin->AddDataFile(data_collection);
  }

  if (!outputSingleFolder.IsNull()) {
    plugin->SetOutputSingleFolder(outputSingleFolder);
    plugin->SetOutputToRunNo();
  }
  plugin->SetJobTag(job_tag);
  plugin->SetNtestFiles(5);
  plugin->SetCheckCopy(kFALSE);
  plugin->SetMergeDirName(mergeDirName);
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion(root_version);
  plugin->SetAliROOTVersion(aliroot_version);
  // Declare input data to be processed.
  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  plugin->SetGridDataDir(alien_datadir);
  // Set data search pattern
  plugin->SetDataPattern(data_pattern);
  if (!useProductionMode) {
    if (runOnData) {
      plugin->SetRunPrefix("%09d");
    }
    //   if (!iAODanalysis) plugin->SetRunRange(run_range[0], run_range[1]);
    for (Int_t i=0; i<10; i++) {
      if (run_numbers[i]==0) break;
      plugin->AddRunNumber(run_numbers[i]);
    }
  }
  // Define alien work directory where all files will be
  // copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir(grid_workdir);
  // Declare alien output directory. Relative to working directory.
  if (alien_outdir.IsNull()) alien_outdir = Form("output_%s",train_name.Data());
  plugin->SetGridOutputDir(alien_outdir);

  TString ana_sources = "";
  TString ana_add = "";
  if (usePAR && anaPars.Length()) {
    printf("%s\n", anaPars.Data());
    TObjArray *arr;
    TObjString *objstr;
    arr = anaPars.Tokenize(" ");
    TIter next(arr);
    while ((objstr=(TObjString*)next())) plugin->EnablePackage(objstr->GetString());
    delete arr;
  }

  // Declare the analysis source files names separated by blancs. To
  // be compiled runtime using ACLiC on the worker nodes.
  ana_sources = ana_sources.Strip();
  // Declare all libraries (other than the default ones for the
  // framework. These will be loaded by the generated analysis
  // macro. Add all extra files (task .cxx/.h) here.
  anaLibs     = anaLibs.Strip();
  if (ana_sources.Length()) plugin->SetAnalysisSource(ana_sources);
  if (anaLibs.Length())     plugin->SetAdditionalLibs(anaLibs);

  // Declare the output file names separated by blancs.  (can be like:
  // file.root or file.root@ALICE::Niham::File)
  plugin->SetDefaultOutputs();
  plugin->SetMergeExcludes(mergeExclude);
  plugin->SetMaxMergeFiles(maxMergeFiles);
  plugin->SetNrunsPerMaster(nRunsPerMaster);
  // Optionally define the files to be archived.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:AliAOD.root,AOD.tag.root@ALICE::NIHAM::File");


  // Put default output files to archive
  TString listhists = "";
  TString listaods  = "";
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TIter next(mgr->GetOutputs());
  AliAnalysisDataContainer *output;
  while ((output=(AliAnalysisDataContainer*)next())) {
    const char *filename = output->GetFileName();
    if (!(strcmp(filename, "default"))) {
      if (!mgr->GetOutputEventHandler()) continue;
      filename = mgr->GetOutputEventHandler()->GetOutputFileName();
      if (listaods.Length()) listaods += ",";
      listaods += filename;
      listaods += ",";
      listaods += "pyxsec_hists.root";
    } else {
      if (!strcmp(filename, "pyxsec_hists.root")) continue;
      if (listhists.Contains(filename)) continue;
      if (listhists.Length()) listhists += ",";
      listhists += filename;
    }
  }
  if (mgr->GetExtraFiles().Length()) {
    if (listaods.Length()) listaods += ",";
    listaods += mgr->GetExtraFiles();
    listaods.ReplaceAll(" ", ",");
  }
  if (listhists.Length()) listhists = Form("hist_archive.zip:%s@%s", 
					   listhists.Data(), 
					   outputStorages.Data());
  if (listaods.Length())  listaods  = Form("aod_archive.zip:%s@%s", 
					   listaods.Data(), 
					   outputStorages.Data());
  if (!listhists.Length() && !listaods.Length()) {
    ::Fatal("AnalysisTrainNew", "No task output !");
  }
  TString outputArchive = Form("log_archive.zip:stderr@%s", 
			       outputStorages.Data());
  if (listaods.Length()) {
    outputArchive += " ";
    outputArchive += listaods;
  }
  if (listhists.Length()) {
    outputArchive += " ";
    outputArchive += listhists;
  }
  // Set friends
  //   if (iAODanalysis && iPWG3d2h)
  //      plugin->SetFriendChainName("AliAOD.VertexingHF.root");
  //   plugin->SetOutputArchive(outputArchive);
  
  // Optionally set a name for the generated analysis macro (default
  // MyAnalysis.C)
  plugin->SetAnalysisMacro(Form("%s.C", train_name.Data()));
  // Optionally set a name for the generated validation script
  plugin->SetValidationScript("PWG2validation.sh");
  // Optionally set maximum number of input files/subjob (default 100,
  // put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(nFilesPerJob);
  // Optionally set number of failed jobs that will trigger killing
  // waiting sub-jobs.
  //   plugin->SetMaxInitFailed(5);
  // Optionally modify the number of replicas
  plugin->SetNumberOfReplicas(4);
  // Optionally resubmit threshold.
  //   plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(70000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName(Form("%s.jdl", train_name.Data()));
  // Optionally modify the executable name (default analysis.sh)
  plugin->SetExecutable(Form("%s.sh", train_name.Data()));
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  // Merge via JDL
  plugin->SetMergeViaJDL(useMergeViaJDL);
  // Use fastread option
  plugin->SetFastReadOption(useFastReadOption);
  // UseOverwrite mode
  plugin->SetOverwriteMode(useOverwriteMode);
  plugin->SetExecutableCommand("aliroot -b -q");
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  return plugin;
}

//______________________________________________________________________________
void WriteConfig()
{
  // Write train configuration in a file. The file name has the
  // format: train_[trainName]_ddMonthyyyy_time.C
  if (useDATE) {
    gSystem->Exec("date +%d%b%Y_%Hh%M > date.tmp");
    ifstream fdate("date.tmp");
    if (!fdate.is_open()) {
      ::Error("AnalysisTrainNew.C::Export","Could not generate file name");
      return;
    }
    const char date[64];
    fdate.getline(date,64);
    fdate.close();
    gSystem->Exec("rm date.tmp");
    train_name = Form("%s_%s", train_name.Data(), date);
  }
  TString cdir = gSystem->WorkingDirectory();
  gSystem->MakeDirectory(train_name);
  gSystem->ChangeDirectory(train_name);
  ofstream out;
  TString outName(Form("%sConfig.C",train_name.Data()));
  out.open(outName.Data(), ios::out);
  if (out.bad()) {
    ::Error("AnalysisTrainNew.C::Export", 
	    "Cannot open %s for writing", outName.Data());
    return;
  }
  out << "{" << endl;
  out << "   train_name      = " << "\"" << train_name.Data() << "\";" << endl;
  out << "   proof_cluster   = " << "\"" << proof_cluster.Data() << "\";" << endl;
  out << "   useAFPAR        = " << useAFPAR << ";" << endl;
  if (useAFPAR)
    out << "   AFversion       = " << AFversion.Data() << ";" << endl;
  out << "   proof_dataset   = " << "\"" << proof_dataset.Data() << "\";" << endl;
  out << "   usePLUGIN       = " << usePLUGIN << ";" << endl;
  out << "   usePAR          = " << usePAR << ";" << endl;
  out << "   useCPAR         = " << useCPAR << ";" << endl;
  out << "   root_version    = " << "\"" << root_version.Data() << "\";" << endl;
  out << "   aliroot_version = " << "\"" << aliroot_version.Data() << "\";" << endl;
  out << "   alien_datadir   = " << "\"" << alien_datadir.Data() << "\";" << endl;
  if (!alien_outdir.Length()) alien_outdir = Form("output_%s",train_name.Data());
  out << "   alien_outdir    = " << "\"" << alien_outdir.Data() << "\";" << endl;
  out << "   maxMergeFiles   = " << maxMergeFiles << ";" << endl;
  out << "   mergeExclude    = " << "\"" << mergeExclude.Data() << "\";" << endl;
  out << "   nRunsPerMaster  = " << nRunsPerMaster << ";" << endl;
  out << "   nFilesPerJob    = " << nFilesPerJob << ";" << endl;
  //   for (Int_t i=0; i<10; i++) {
  //      if (run_numbers[i])
  //         out << "   run_numbers[" << i << "]  = " << run_numbers[i] << ";" << endl;
  //   }
  //   out << "   run_range[0]    = " << run_range[0] << ";" << endl;
  //   out << "   run_range[1]    = " << run_range[1] << ";" << endl;
  out << "   usePhysicsSelection = " << usePhysicsSelection << ";" << endl;
  out << "   useTender       = " << useTender << ";" << endl;
  out << "   useMergeViaJDL  = " << useMergeViaJDL << ";" << endl;
  out << "   useOverwriteMode  = " << useOverwriteMode << ";" << endl;
  out << "   useFastReadOption = " << useFastReadOption << ";" << endl;
  out << "   useDBG          = " << useDBG << ";" << endl;
  out << "   useMC           = " << useMC << ";" << endl;
  out << "   useTAGS         = " << useTAGS << ";" << endl;
  out << "   useKFILTER      = " << useKFILTER << ";" << endl;
  out << "   useTR           = " << useTR << ";" << endl;
  out << "   useCORRFW       = " << useCORRFW << ";" << endl;
  out << "   useAODTAGS      = " << useAODTAGS << ";" << endl;
  out << "   saveTrain       = " << "kFALSE;" << endl << endl;
  out << "   // Analysis modules" << endl;
  out << "   iAODanalysis    = " << iAODanalysis << ";" << endl;
  out << "   iAODhandler     = " << iAODhandler << ";" << endl;
  out << "   iPWG2fmd        = " << iPWG2fmd << ";" << endl;
  out << "   iPWG2femto      = " << iPWG2femto << ";" << endl;
  out << "   iPWG2spectra    = " << iPWG2spectra << ";" << endl;
  out << "}" << endl;
  ::Info("AnalysisTrainNew.C::WriteConfig", 
	 "Train configuration wrote to file %s", outName.Data());
  gSystem->ChangeDirectory(cdir);
}

//____________________________________________________________________
Bool_t LoadConfig(const char *filename)
{
  // Read train configuration from file
  if (gSystem->AccessPathName(filename)) {
    ::Error("AnalysisTrainNew.C::LoadConfig", "Config file name not found");
    return kFALSE;
  }
  gROOT->ProcessLine(Form(".x %s", filename));
  ::Info("AnalysisTrainNew.C::LoadConfig", 
	 "Train configuration loaded from file %s", filename);
  return kTRUE;
}
//____________________________________________________________________
//
// EOF
//
