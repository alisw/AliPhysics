//____________________________________________________________________
/** 
 * Run an FMD corrections job. 
 * 
 * @param runMode     Running mode (full, test, terminate, submit, offline)
 * @param anaType     What to do (background, collector, sharing)
 * @param dataDir     Input data directory 
 * @param anaSource   Analysis source (if any)
 * @param anaName     Analysis name 
 * @param colSys      Collision system (p-p, Pb-Pb, A-A)
 * @param cmsNNGeV    Center of mass energy per nucleon in GeV
 * @param bkG         Magnetic field in kilo gaus 
 * @param aliceTag    AliROOT tag to use 
 * @param rootTag     ROOT tag to use 
 * @param apiTag      API tag to use 
 */
void 
runFMDjob(const TString& runMode   = "",
	  const TString& anaType   = "background",
	  const TString& dataDir   = "",
	  const TString& anaSource = "",
	  const TString& anaName   = "", 
	  const TString& colSys    = "p-p", 
	  Float_t        cmsNNGeV  = 900, 
	  Float_t        bkG       = 5,
	  const TString& aliceTag  = "v4-21-04-AN", 
	  const TString& rootTag   = "v5-27-06b", 
	  const TString& apiTag    = "V1.1x")
{
  // --- Load libraries needed  -------------------------------------
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGLFforward");
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGLF/FORWARD/analysis");

  // --- Some initial checks and setup -------------------------------
  TString outFileName = anaName;
  outFileName.ToLower();
  outFileName += ".root";

  // --- Manager setup -----------------------------------------------
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager(anaName.Data());

  
  // --- ESD setup ---------------------------------------------------
  // Connect the EDS's to the manager and switch off unused branches
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetInactiveBranches("AliESDACORDE "
			    "AliRawDataErrorLogs "
			    "CaloClusters "
			    "Cascades "
			    "EMCALCells "
			    "EMCALTrigger "
			    "ESDfriend "
			    "Kinks "
			    "Cascades "
			    "ALIESDACORDE "
			    "MuonTracks "
			    "TrdTracks");
  mgr->SetInputEventHandler(esdH);

  // --- Task setup --------------------------------------------------
  // Configure the analysis manager to the specific task
  TString type = anaType.ToLower();
  TString addLibs;
  if (type.Contains("background")) 
    addLibs = AddBackgroundTask(outFileName);
  else if (type.Contains("collector")) 
    addLibs = AddCollectorTask(outFileName, colSys, cmsNNGeV, bkG);
  else if (type.Contains("sharing")) 
    addLibs = AddSharingTask(outFileName, colsys, cmsNNGeV, bkG);
  else {
    Error("runFMDjob", "Unknown type '%s', please fix this macro", 
	  anaType);
    return;
  }

  // --- AliEN handler setup -----------------------------------------
  // Create and configure the alien handler plugin, and connect to manager 
  if (!CreateAlienHandler(runMode, 
			  dataDir, 
			  anaSource, 
			  addLibs, 
			  anaName, 
			  outFileName,
			  aliceTag,
			  rootTag, 
			  apiTag)) {
    Error("runFMDjob", "Failed to set up alien handler");
    return;
  }
  mgr->SetGridHandler(alienHandler);
  
  // --- final job setup and execution -------------------------------
  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis()) {
    Error("runFMDjob", "Failed to initialise the train");
    return;
  }
  mgr->PrintStatus();
  mgr->StartAnalysis("grid"); 

  // --- We are done -------------------------------------------------
  Info("runFMDjob", "Job is done");
}


//____________________________________________________________________
/** 
 * Create a background correction task 
 * 
 * @param outFileName 
 *
 * @return A list of additional files that should be uploaded 
 */
void 
AddBackgroundTask(const TString& outFileName) 
{
  AliAnalysisManager* mgr = AliAnalysisManager::Instance();

  // --- Make and configure task -------------------------------------
  AliFMDAnalysisTaskGenerateCorrection *task = 
    new AliFMDAnalysisTaskGenerateCorrection("background");
  task->SetNbinsEta(200);
  mgr->AddTask(task);
  
  // --- Add the MC handler ------------------------------------------
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  
  // --- Create and connect container for input ----------------------
  AliAnalysisDataContainer *cin_esd   = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,  0, cin_esd);

  // --- Create and connect containers for output --------------------
  const char*  outList[] = { "Hits", "Primaries", "vertex", "Correction", 0 };
  const char** ptr       = outList;
  Int_t i = 1; 
  while (*ptr) { 
    AliAnalysisDataContainer* outCont = 
      mgr->CreateContainer(*ptr, TList::Class(), 
			   AliAnalysisManager::kOutputContainer, 
			   outFileName.Data());
    mgr->ConnectOutput(task, i, outCont);
    i++;
    ptr++;
  }

  return "";
}


//_______________________________________________________________
/** 
 * Create and add collector task. 
 * 
 * @param outFileName  Output file name 
 * @param col          Collision system (one of "p-p", "Pb-Pb", or "A-A")
 * @param cmsNNGeV     Center of mass energy per nucleon in GeV
 * @param bkG          Magnetic field in kilo gauss 
 *
 * @return A list of additional files that should be uploaded 
 */
void AddCollectorTask(const TString& outFileName,
		      const TString& col, 
		      Float_t        cmsNNGeV, 
		      Float_t        bkG) 
{
  AliAnalysisManager* mgr = AliAnalysisManager::Instance();

  // --- Initialize the analysis parameters --------------------------
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();  
  pars->SetCollisionSystem(col);
  pars->SetEnergy(cmsNNGeV);
  pars->SetMagField(bkG);
  pars->SetBackgroundPath("./");
  pars->Init(kTRUE,AliFMDAnaParameters::kBackgroundCorrection);
  TString bgFile = pars->GetPath(AliFMDAnaParameters::GetBackgroundID());

  // --- Create and add our task -------------------------------------
  AliFMDAnalysisTaskCollector *task = 
    new AliFMDAnalysisTaskCollector("collector");
  mgr->AddTask(task);

  // --- Add the MC handler ------------------------------------------
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  
  // --- Create containers for input/output --------------------------
  AliAnalysisDataContainer *cin_esd   = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *cout_fmd1 = 
    mgr->CreateContainer("energyDist", 
			 TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 outFileName.Data());
  
  // --- Connect input/output ----------------------------------------
  mgr->ConnectInput(task,  0, cin_esd);
  mgr->ConnectOutput(task, 1, cout_fmd1);

  return bgFile;
}

//_______________________________________________________________
/** 
 * Create and add sharing task. 
 * 
 * @param outFileName  Output file name 
 * @param col          Collision system (one of "p-p", "Pb-Pb", or "A-A")
 * @param cmsNNGeV     Center of mass energy per nucleon in GeV
 * @param bkG          Magnetic field in kilo gauss 
 *
 * @return A list of additional files that should be uploaded 
 */
TString 
AddSharingTask(const TString& outFileName,
	       const TString& col, 
	       Float_t        cmsNNGeV, 
	       Float_t        bkG) 
{
  AliAnalysisManager* mgr = AliAnalysisManager::Instance();

  // --- Initialize the analysis parameters --------------------------
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();  
  pars->SetCollisionSystem(col);
  pars->SetEnergy(cmsNNGeV);
  pars->SetMagField(bkG);
  pars->SetBackgroundPath("./");
  pars->SetEnergyPath("./");
  pars->SetEventSelectionPath("./");
  pars->Init(kTRUE,
	     AliFMDAnaParameters::kBackgroundCorrection|
	     AliFMDAnaParameters::kEnergyDistributions|
	     AliFMDAnaParameters::kEventSelectionEfficiency);
  TString files;
  files.Append(pars->GetPath(AliFMDAnaParameters::GetBackgroundID()));
  files.Append(" ");
  files.Append(pars->GetPath(AliFMDAnaParameters::GetEdistID()));
  files.Append(" ");
  files.Append(pars->GetPath(AliFMDAnaParameters::GetEventSelectionEffID()));
  
  // --- Create and add our task -------------------------------------
  AliFMDAnalysisTaskSE *task = new AliFMDAnalysisTaskSE("sharing");
  mgr->AddTask(task);

  // --- Add the MC handler ------------------------------------------
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  


  // --- Create and connect containers for output --------------------
  AliAnalysisDataContainer *cout_fmd = 
    mgr->CreateContainer("BackgroundCorrected", TList::Class(), 
                         AliAnalysisManager::kOutputContainer,outputfile); 
  
  mgr->ConnectInput(taskfmd, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskfmd, 1, cout_fmd);

  return files;
}

//____________________________________________________________________
/** 
 * Create an AliAnalysisGrid parameter object 
 * 
 * @param runMode     Running mode (full, test, terminate, submit, offline)
 * @param dataDir     Input data directory 
 * @param anaSource   Possible source to compile on worker node 
 *                    (must also be compiled and addet to train on 
 *                    submitter machine)
 * @param addLibs     Extra libraries to add 
 * @param anaName     Analysis name (i.e., script created)
 * @param outFileName Output file name 
 * 
 * @return Valid object or null
 */
Bool_t
CreateAlienHandler(const TString& runMode,
		   const TString& dataDir,
		   const TString& anaSource,
		   const TString& addLibs,
		   const TString& anaName,
		   const TString& outFileName,
		   const TString& aliceTag, 
		   const TString& rootTag, 
		   const TString& apiTag) 
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  // Overwrite all generated files, datasets and output 
  // results from a previous session
  plugin->SetOverwriteMode();

  // Set the running mode 
  plugin->SetRunMode(runMode.Data());

  // Add path to our headers 
  plugin->AddIncludePath("-I$ALICE_ROOT/PWGLF/FORWARD/analysis");

  // Set versions of used packages
  plugin->SetAPIVersion(apiTag);
  plugin->SetROOTVersion(rootTag);
  plugin->SetAliROOTVersion(aliceTag);

  // Define production directory LFN
  plugin->SetGridDataDir(dataDir.Data());

  // Set data search pattern
  plugin->SetDataPattern("AliESDs.root");
  
  // Use ESD tags (same applies for AOD's)
  //plugin->SetDataPattern("*tag.root");  
  
  // ...then add run numbers to be considered
  // If not set all runs proccessed
  //plugin->AddRunNumber(126437); 

  // Set events to run over for each file !!!
  //plugin->SetRunRange(0, 10); 
  
  // Define alien work directory where all files will be copied. 
  // Relative to alien $HOME.
  plugin->SetGridWorkingDir("work");
  
  // Declare alien output directory. Relative to working directory.
  TString outputDir = anaName;
  outputDir += "_out";
  plugin->SetGridOutputDir(outputDir.Data());
  
  // Declare the analysis source files names separated by blancs. 
  // To be compiled runtime using ACLiC on the worker nodes.
  if (!anaSource.IsNull())
    plugin->SetAnalysisSource(anaSource.Data());
  
  // Declare all libraries (other than the default ones for the framework. 
  // These will be loaded by the generated analysis macro. 
  // Add all extra files (task .cxx/.h) here.
  if (!addLibs.IsNull()) 
    plugin->SetAdditionalLibs(addLibs.Data());
  
  // No need for output file names. Procedure is automatic.
  // It's works better this way
  plugin->SetDefaultOutputs(kFALSE);
  plugin->SetOutputFiles(outFileName.Data());

  // Set a name for the generated analysis macro (default MyAnalysis.C).
  // Make this unique !!!
  TString macroName = anaName;
  macroName += "Task.C";
  plugin->SetAnalysisMacro(macroName.Data());
  
  // Optionally set maximum number of input files/subjob (default 100,
  // put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(100);

  // Optionally set number of failed jobs that will trigger killing
  // waiting sub-jobs.
  plugin->SetMaxInitFailed(5);

  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);

  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(20000);

  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");

  // Optionally modify the name of the generated JDL (default analysis.jdl)
  TString jdlName = anaName;
  jdlName += ".jdl";
  plugin->SetJDLName(jdlName.Data());
  
  // Optionally modify job price (default 1)
  plugin->SetPrice(1); 
  
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se"); 
  
  // connect to manager 
  AliAnalysisManager* mgr = AliAnalysisManager::Instance();
  mgr->SetGridHandler(alienHandler);
  
  return kTRUE
} 
//____________________________________________________________________
//
// EOF
//
