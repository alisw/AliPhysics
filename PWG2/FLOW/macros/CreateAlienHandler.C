AliAnalysisGrid* CreateAlienHandler() {
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init 
  // then source /tmp/gclient_env_$UID in the current shell.
  if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode("test");
  plugin->SetNtestFiles(2); // Relevant only for run mode "test" 

  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-26-00b-2");
  plugin->SetAliROOTVersion("v4-19-10-AN");

  // Declare input data to be processed - can be done in two ways:
  // METHOD 1: Create automatically XML collections using alien 'find' command.
  // ============================================================================
  //  Example 1: MC production (set in macro runFlowTask.C: DATA = kFALSE)
  // plugin->SetGridDataDir("/alice/sim/LHC10a10");
  // plugin->SetDataPattern("*AliESDs.root"); // The default data pattern, other may be "*tag.root", "*ESD.tag.root", etc
  // plugin->AddRunNumber(105054); // Alternatively use e.g. plugin->SetRunRange(105044,106044); to add more runs in one go  
  // plugin->SetOutputToRunNo();  
  // ============================================================================
  //  Example 2: Real data (set in macro runFlowTask.C: DATA = kTRUE, MCEP = kFALSE)
  plugin->SetGridDataDir("/alice/data/2009/LHC09d");
  plugin->SetDataPattern("*ESDs/pass5/*ESDs.root");
  plugin->SetRunPrefix("000");
  plugin->AddRunNumber(104044); 
  plugin->AddRunNumber(104892); // Alternatively use e.g. plugin->SetRunRange(104044,106044); to add more runs in one go 
  plugin->SetOutputToRunNo();  
  // ============================================================================
 
  // METHOD 2: Declare existing data files (raw collections, xml collections, root file)
  // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
  // XML collections added via this method can be combined with the first method if
  // the content is compatible (using or not tags)
  //plugin->AddDataFile("hijingWithoutFlow10000Evts.xml");
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir("realData/900GeV");
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes:
  // ... (if this is needed see in official tutorial example how to do it!)

  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalLibs("libPWG2flowCommon.so libPWG2flowTasks.so");
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  plugin->SetOutputFiles("AnalysisResults.root ");
  // Optionally define the files to be archived.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
  //plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  plugin->SetOutputArchive("log_archive.zip:");
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("flowAnalysis.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(100);
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  plugin->SetMaxInitFailed(5);
  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(20000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("flowAnalysis.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");

  return plugin;
}
