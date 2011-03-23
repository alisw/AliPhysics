AliAnalysisGrid* CreateAlienHandlerCaloEtSim(TString outputDir, TString outputName, const char * pluginRunMode)
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
  if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  // Overwrite all generated files, datasets and output results from a previous session
  plugin->SetOverwriteMode();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  // plugin->SetRunMode("full");  // VERY IMPORTANT - DECRIBED BELOW
  // plugin->SetRunMode("test");  // VERY IMPORTANT - DECRIBED BELOW
  plugin->SetRunMode(pluginRunMode);  // VERY IMPORTANT - DECRIBED BELOW
  cout<<"Running in "<<pluginRunMode<<" mode"<<endl;

  // Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-28-00a");
   plugin->SetAliROOTVersion("v4-21-17b-AN");
  // Declare input data to be processed.

  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  //   plugin->SetGridDataDir("/alice/sim/LHC10a18");
  // Set data search pattern
  //   plugin->SetDataPattern("*ESDs.root");  // simulated, tags not used
  //   plugin->SetDataPattern("*ESDs/pass4/*ESDs.root"); // real data check reco pass and data base directory
  //   plugin->SetRunPrefix("000");   // real data
  //   plugin->SetDataPattern("*tag.root");  // Use ESD tags (same applies for AOD's)
  // ...then add run numbers to be considered
  //   plugin->AddRunNumber(125020);    // simulated
  //   plugin->AddRunNumber(104065);  // real data

   plugin->SetGridDataDir("/alice/sim/LHC10d4");
   plugin->SetDataPattern("*ESDs.root");
   plugin->AddRunNumber("120741");//smallest of the above
  // Method 2: Declare existing data files (raw collections, xml collections, root file)
  // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
  // XML collections added via this method can be combined with the first method if
  // the content is compatible (using or not tags)
  //plugin->AddDataFile("tag.xml");
  // plugin->AddDataFile("wn.xml"); // test
  // file generated with:  find -x tag /alice/sim/LHC10d1/117222/* AliESDs.root > tag.xml

  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir(outputDir.Data());
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime IN THE SAME ORDER THEY ARE LISTED
  // using ACLiC on the worker nodes.
  plugin->SetAnalysisSource("AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.cxx AliAnalysisEtCommon.cxx AliAnalysisEt.cxx AliAnalysisEtMonteCarlo.cxx AliAnalysisEtMonteCarloPhos.cxx AliAnalysisEtMonteCarloEmcal.cxx AliAnalysisEtReconstructed.cxx AliAnalysisEtReconstructedPhos.cxx AliAnalysisEtReconstructedEmcal.cxx AliAnalysisEtSelectionContainer.cxx AliAnalysisEtSelectionHandler.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisTaskTotEt.cxx");
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalLibs("AliAnalysisEtCuts.h AliAnalysisEt.h AliAnalysisEtMonteCarlo.h AliAnalysisEtMonteCarloPhos.h AliAnalysisEtMonteCarloEmcal.h AliAnalysisEtReconstructed.h AliAnalysisEtReconstructedPhos.h AliAnalysisEtReconstructedEmcal.h AliAnalysisTaskTotEt.h AliAnalysisEtCuts.cxx AliAnalysisEt.cxx AliAnalysisEtMonteCarlo.cxx AliAnalysisEtMonteCarloPhos.cxx AliAnalysisEtMonteCarloEmcal.cxx AliAnalysisEtReconstructed.cxx AliAnalysisEtReconstructedPhos.cxx AliAnalysisEtReconstructedEmcal.cxx AliAnalysisTaskTotEt.cxx  AliAnalysisEtCommon.cxx  AliAnalysisEtCommon.h AliAnalysisHadEtCorrections.h AliAnalysisHadEtCorrections.cxx AliAnalysisEtSelectionContainer.h AliAnalysisEtSelectionHandler.h AliAnalysisTaskTransverseEnergy.h AliAnalysisEtSelectionContainer.cxx AliAnalysisEtSelectionHandler.cxx AliAnalysisTaskTransverseEnergy.cxx ConfigEtMonteCarlo.C ConfigEtReconstructed.C physicsSelections.root corrections.root");
  plugin->SetExecutableCommand("aliroot -b -q");
  // add extra include files/path
  plugin->AddIncludePath("-I. -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS");

  // No need for output file names. Procedure is automatic. <-- not true
  plugin->SetDefaultOutputs(kFALSE);
  //plugin->SetOutputFiles(outputName.Data());
  plugin->SetOutputFiles("Et.ESD.sim.EMCAL.root event_stat.root");
  // No need define the files to be archived. Note that this is handled automatically by the plugin.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  // Set a name for the generated analysis macro (default MyAnalysis.C) Make this unique !
  plugin->SetAnalysisMacro("DavidEtAnalysis.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore). The optimum for an analysis
  // is correlated with the run time - count few hours TTL per job, not minutes !
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
  plugin->SetJDLName("TaskEt.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1); 
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");

  return plugin;
} 
