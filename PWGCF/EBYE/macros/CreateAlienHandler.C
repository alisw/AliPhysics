AliAnalysisGrid* CreateAlienHandler() {
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  //plugin->SetRunMode("test");
  //plugin->SetRunMode("offline");
  //plugin->SetRunMode("submit");
  plugin->SetRunMode("full");
  //plugin->SetRunMode("terminate");
  plugin->SetNtestFiles(3); // Relevant only for run mode "test" 

  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-27-06d");
  plugin->SetAliROOTVersion("v4-21-16-AN");  
  
  // Declare input data to be processed - can be done in two ways:
  // METHOD 1: Create automatically XML collections using alien 'find' command.
  // ============================================================================
  //  Example 1: MC production (set in macro runFlowTask.C: DATA = kFALSE)
  //plugin->SetGridDataDir("/alice/sim/LHC10d4");
  //plugin->SetDataPattern("*AliESDs.root"); // The default data pattern, other may be "*tag.root", "*ESD.tag.root", etc
  //plugin->AddRunNumber(119844); // Alternatively use e.g. plugin->SetRunRange(105044,106044); to add more runs in one go  
  //plugin->SetOutputToRunNo(); 
  // ============================================================================
  //  Example 2: Real data (set in macro runFlowTask.C: DATA = kTRUE, MCEP = kFALSE)
  plugin->SetGridDataDir("/alice/data/2010/LHC10h");
  plugin->SetDataPattern("*ESDs/pass1_4plus/*ESDs.root");
  plugin->SetRunPrefix("000"); // IMPORTANT!
  plugin->AddRunNumber(137161);
  plugin->AddRunNumber(137431);
  plugin->AddRunNumber(137549);
  plugin->AddRunNumber(137595);
  plugin->AddRunNumber(137638);
  plugin->AddRunNumber(137639);
  plugin->AddRunNumber(137693);
  // plugin->AddRunNumber(119844); // Alternatively use e.g. plugin->SetRunRange(104044,106044); to add more runs in one go 
  plugin->SetOutputToRunNo();  

  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir("Fluctuations/PbPb/2.76TeV/LHC10h/Pass1_4Plus/Systematics/Centrality/TPC/TPCOnly");
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  plugin->SetAnalysisSource("AliEbyEFluctuationAnalysisTask.cxx");
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalLibs("AliEbyEFluctuationAnalysisTask.h AliEbyEFluctuationAnalysisTask.cxx");
  
  // Do not specify your outputs by hand anymore:
  plugin->SetDefaultOutputs(kTRUE);
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("fluctuationsAnalysis.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(100);
  // Optionally set number of runs per masterjob:
  plugin->SetNrunsPerMaster(7);
 
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("TaskFluctuations.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");

  //Merging
  plugin->SetMergeViaJDL(kTRUE);

  return plugin;
}
