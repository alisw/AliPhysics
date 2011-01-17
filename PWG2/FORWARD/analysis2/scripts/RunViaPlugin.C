class AliAnalysisAlien;

//____________________________________________________________________
// Forward declarations 
AliAnalysisAlien*
CreateAlienHandler(const TString& runMode,
		   const TString& dataDir,
		   const TString& anaSource,
		   const TString& addLibs,
		   const TString& anaName,
		   const TString& aliceTag, 
		   const TString& rootTag, 
		   const TString& apiTag);

//____________________________________________________________________
void
RunViaPlugin(const Char_t* runMode="", 
	     const Char_t* dataDir=".", 
	     Long64_t      nEvents=-1, 
	     Bool_t        mc=false)
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", 
						    "FMD analysis train");

  // --- ESD input handler -------------------------------------------
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  esdHandler->SetInactiveBranches("AliESDACORDE "
				  "AliRawDataErrorLogs "
				  "CaloClusters "
				  "Cascades "
				  "EMCALCells "
				  "EMCALTrigger "
				  "Kinks "
				  "Cascades "
				  "MuonTracks "
				  "TrdTracks "
				  "HLTGlobalTrigger");
  mgr->SetInputEventHandler(esdHandler);      
       
  // --- Monte Carlo handler -----------------------------------------
  if (mc) {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
    mcHandler->SetReadTR(true);    
  }

  // --- AOD output handler ------------------------------------------
  AliAODHandler* aodHandler   = new AliAODHandler();
  mgr->SetOutputEventHandler(aodHandler);
  aodHandler->SetOutputFileName("AliAODs.root");

  // --- Add tasks ---------------------------------------------------
  // Physics selection 
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AddTaskPhysicsSelection(mc, kTRUE, kTRUE);

  // FMD 
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskFMD.C");
  AddTaskFMD(mc);

  // --- Create the plug-in object -----------------------------------
  TString mode(runMode);
  TString dir(dataDir);
  TString anaSource("");
  TString addLibs("");
  TString anaName("FMD");
  TString aliceTag("v4-21-04-AN");
  TString rootTag("v5-27-06b");
  TString apiTag("V1.1x");
  AliAnalysisAlien* alienHandler = CreateAlienHandler(mode,
						      dir,
						      anaSource,
						      addLibs,
						      anaName,
						      aliceTag, 
						      rootTag, 
						      apiTag);
  if (!alienHandler) { 
    Error("RunViaPlugin.C", "Failed to make plugin");
    return;
  }
  mgr->SetGridHandler(alienHandler);

  // --- final job setup and execution -------------------------------
  // Enable debug printouts
  // mgr->SetDebugLevel(2);
  if (!mgr->InitAnalysis()) {
    Error("RunViaPlugin.C", "Failed to initialise the train");
    return;
  }
  if (nEvents <= 0) nEvents = 1234567890;
  TString amode("grid");
  mode.ToLower();
  if (mode.Contains("test")) amode = local;

  mgr->PrintStatus();
  mgr->StartAnalysis(amode.Data(), nEvents);  
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
 * @param aliceTag    Tag on AliROOT
 * @param rootTag     Tag on ROOT
 * @param apiTag      AliEN tag
 * 
 * @return Valid object or null
 */
AliAnalysisAlien*
CreateAlienHandler(const TString& runMode,
		   const TString& dataDir,
		   const TString& anaSource,
		   const TString& addLibs,
		   const TString& anaName,
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
  plugin->AddIncludePath("-I$ALICE_ROOT/PWG2/FORWARD/analysis2");

  // Set versions of used packages
  plugin->SetAPIVersion(apiTag);
  plugin->SetROOTVersion(rootTag);
  plugin->SetAliROOTVersion(aliceTag);

  // Define production directory LFN
  plugin->SetGridDataDir(dataDir.Data());

  // Set data search pattern
  plugin->SetDataPattern("*ESDs*.root");
  
  // Use ESD tags (same applies for AOD's)
  //plugin->SetDataPattern("*tag.root");  

  // File used in test mode 
  plugin->SetFileForTestMode("testFiles");

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
  plugin->SetDefaultOutputs(kTRUE);

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
  // AliAnalysisManager* mgr = AliAnalysisManager::Instance();
  // mgr->SetGridHandler(plugin);
  
  return plugin;
}
     
//____________________________________________________________________
//
// EOF
//
