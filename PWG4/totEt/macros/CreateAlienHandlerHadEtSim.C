 AliAnalysisGrid* CreateAlienHandlerHadEtSim()
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();

// Overwrite all generated files, datasets and output results from a previous session
   plugin->SetOverwriteMode();
// Set the run modSoon a picture of Kim Jong-un (sitting far left) flashed around the world, the first known image of him since his school dayse (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode("full");  // VERY IMPORTANT - DECRIBED BELOW
   //plugin->SetRunMode("test");  // VERY IMPORTANT - DECRIBED BELOW
   //plugin->SetCheckCopy(kFALSE);
// Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-27-06b");
   plugin->SetAliROOTVersion("v4-21-10-AN");
// Declare input data to be processed.

// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
//   plugin->SetGridDataDir("/alice/data/2010/LHC10d");
// Set data search pattern
//   plugin->SetDataPattern("*ESDs.root");  // simulated, tags not used
//   plugin->SetDataPattern("*ESDs/pass2/*ESDs.root"); // real data check reco pass and data base directory
//   plugin->SetRunPrefix("000");   // real data
//   plugin->SetDataPattern("*tag.root");  // Use ESD tags (same applies for AOD's)
// ...then add run numbers to be considered
//   plugin->AddRunNumber(125020);    // simulated
//   plugin->AddRunNumber(126403);  // real data
//   plugin->AddRunNumber(126404);  // real data
//   plugin->AddRunNumber(126405);  // real data

   plugin->SetGridDataDir("/alice/sim/LHC10d4");
   plugin->SetDataPattern("*ESDs.root");
   plugin->AddRunNumber(120741);//smallest of the above
//    plugin->AddRunNumber(120750);
//    plugin->AddRunNumber(120758);
//    plugin->AddRunNumber(120820);
//    plugin->AddRunNumber(120821);
//    plugin->AddRunNumber(120822);
//    plugin->AddRunNumber(120823);
//    plugin->AddRunNumber(120824);
//    plugin->AddRunNumber(120825);
//    plugin->AddRunNumber(120829);

// Method 2: Declare existing data files (raw collections, xml collections, root file)
// If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
// XML collections added via this method can be combined with the first method if
// the content is compatible (using or not tags)
//   plugin->AddDataFile("tag.xml");
//   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");

// Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("et");
// Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
// Declare the analysis source files names separated by blancs. To be compiled runtime IN THE SAME ORDER THEY ARE LISTED
// using ACLiC on the worker nodes.
   //plugin->SetAnalysisSource("AliAnalysisTaskHadEt.cxx");
   //plugin->SetAnalysisSource("AliAnalysisEt.cxx AliAnalysisEtMonteCarlo.cxx AliAnalysisEtMonteCarloPhos.cxx AliAnalysisEtReconstructed.cxx AliAnalysisEtReconstructedPhos.cxx AliAnalysisHadEt.cxx AliAnalysisHadEtMonteCarlo.cxx AliAnalysisHadEtReconstructed.cxx AliAnalysisTaskHadEt.cxx AliAnalysisTaskTotEt.cxx");
   plugin->SetAnalysisSource("AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.cxx AliAnalysisEtCommon.cxx AliAnalysisHadEt.cxx AliAnalysisHadEtMonteCarlo.cxx AliAnalysisHadEtReconstructed.cxx AliAnalysisEtSelectionContainer.cxx AliAnalysisEtSelectionHandler.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisTaskHadEt.cxx");
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->SetAdditionalLibs("AliAnalysisEtCuts.h AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.h AliAnalysisHadEtCorrections.cxx AliAnalysisEtCommon.h AliAnalysisEtCommon.cxx AliAnalysisHadEt.cxx AliAnalysisHadEtMonteCarlo.cxx AliAnalysisHadEtReconstructed.cxx  AliAnalysisEtSelectionContainer.cxx AliAnalysisEtSelectionHandler.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisTaskHadEt.cxx AliAnalysisHadEt.h AliAnalysisHadEtMonteCarlo.h AliAnalysisHadEtReconstructed.h  AliAnalysisEtSelectionContainer.h AliAnalysisEtSelectionHandler.h AliAnalysisTaskTransverseEnergy.h AliAnalysisTaskHadEt.h physicsSelections.root corrections.root ConfigHadEtAnalysis.C ConfigHadEtMonteCarlo.C ConfigHadEtReconstructed.C");
// No need for output file names. Procedure is automatic. <-- not true
   plugin->SetDefaultOutputs(kFALSE);
   plugin->SetOutputFiles("Et.ESD.new.sim.root event_stat.root");
// No need define the files to be archived. Note that this is handled automatically by the plugin.
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
// Set a name for the generated analysis macro (default MyAnalysis.C) Make this unique !
   plugin->SetAnalysisMacro("ChristinesEtAnalysis.C");
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
