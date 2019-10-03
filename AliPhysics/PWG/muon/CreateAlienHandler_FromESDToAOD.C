AliAnalysisGrid* CreateAlienHandler_FromESDToAOD()
{
//========================================================================
// Macro to configure the GRID plugin
// (see Alice offline web pages for definitions)
//========================================================================

// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell
//=====================================================================

   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();

// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
//=====================================================================
   plugin->SetRunMode("full");  // VERY IMPORTANT 

// Set versions of used packages
//=====================================================================
   plugin->SetAPIVersion("V2.4");
   plugin->SetROOTVersion("v5-23-04");
   plugin->SetAliROOTVersion("v4-17-03");

// Declare input data to be processed.
//=====================================================================

// Method 1: Create automatically XML collections using alien 'find' command.
//===========
// Define production directory LFN
//   plugin->SetGridDataDir("/alice/cern.ch/user/a/arnaldi/FromESDToAOD/JPSI_generation/1001");
// Set data search pattern
//   plugin->SetDataPattern("*tag.root");
// ...then add run numbers to be considered
//   plugin->AddRunNumber(300000);
//   plugin->AddRunNumber(1001);

// Method 2: Declare existing data files (raw collections, xml collections, root file)
//===========
// If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
// XML collections added via this method can be combined with the first method if
// the content is compatible (using or not tags)
   plugin->AddDataFile("/alice/cern.ch/user/a/arnaldi/FromESDToAOD/Plugin/FromESDToAOD.xml");

// Define alien work directory where all files will be copied. Relative to alien $HOME.
//=====================================================================
   plugin->SetGridWorkingDir("FromESDToAOD/Plugin");
   
// Declare alien output directory. Relative to working directory.
//=====================================================================
   plugin->SetGridOutputDir("outputPlugin"); // In this case will be $HOME/work/output
   
// Declare the analysis source files names separated by blancs. 
// Declare all libraries (other than the default ones for the framework). These will be
// loaded by the generated analysis macro and compiled runtime.
// Add par files, if needed.
// Add all extra files (task .cxx/.h/.C) here.
// (AddTaskTagCreation.C can be removed from SetAdditionalLibs, if available in the grid aliroot version)
//=====================================================================
   plugin->SetAdditionalLibs("libPWGHFbase.so libPWGmuon.so AddTaskTagCreation.C");

// Declare the output file names separated by blancs.
//=====================================================================
// (can be like: file.root or file.root@ALICE::Niham::File)
   plugin->SetOutputFiles("AliAODs.root AOD.tag.root");

// Optionally define the files to be archived.
//=====================================================================
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");

// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
//=====================================================================
   plugin->SetAnalysisMacro("analysisFromESDToAOD_Plugin.C");

// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
//=====================================================================
   plugin->SetSplitMaxInputFileNumber(0);

// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//=====================================================================
   plugin->SetMaxInitFailed(5);

// Optionally resubmit threshold.
//=====================================================================
   plugin->SetMasterResubmitThreshold(90);

// Optionally set time to live (default 30000 sec)
//=====================================================================
   plugin->SetTTL(20000);

// Optionally set input format (default xml-single)
//=====================================================================
   plugin->SetInputFormat("xml-single");

// Optionally modify the name of the generated JDL (default analysis.jdl)
//=====================================================================
   plugin->SetJDLName("analysisFromESDToAOD_Plugin.jdl");

// Optionally modify job price (default 1)
//=====================================================================
   plugin->SetPrice(1); 

// Optionally modify split mode (default 'se')    
//=====================================================================
   plugin->SetSplitMode("se");

// Optionally define preferred SE
//=====================================================================
   plugin->SetPreferedSE("ALICE::Torino::DPM");
   return plugin;
} 
