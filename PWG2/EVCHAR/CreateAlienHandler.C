AliAnalysisGrid* CreateAlienHandler(const char * filename)
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
//   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   plugin->SetOverwriteMode();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode("full");
   plugin->SetNtestFiles(3);
// Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
//    plugin->SetROOTVersion("v5-27-05");
//    plugin->SetAliROOTVersion("v4-20-08-AN");
   //change: 28/09/10
   //VO_ALICE@AliRoot::v4-20-09-AN-build1 VO_ALICE@ROOT::v5-27-05-build3,VO_ALICE@GEANT3::v1-11-14 
   plugin->SetROOTVersion("v5-27-05-build3");
   plugin->SetAliROOTVersion("v4-20-09-AN-build1");

// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
// LHC10f
// /alice/sim/LHC10f4a/198000  
//   plugin->SetGridDataDir("/alice/sim/LHC10f4a");
   TString filenameStr = Form("%s",filename);
   filenameStr = TString("/alice/sim/")+filenameStr;
   const char * inputfilename = filenameStr.Data();
   plugin->SetGridDataDir(inputfilename);
// Set data search pattern
   plugin->SetDataPattern("*ESDs.root");
// Data pattern for reconstructed data
//   plugin->SetDataPattern("AliESDs.root");
// ...then add run numbers to be considered
   plugin->AddRunNumber(130844);
   //   plugin->AddRunNumber(198000);
// Method 2: Declare existing data files (raw collections, xml collections, root file)
// If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
// XML collections added via this method can be combined with the first method if
// the content is compatible (using or not tags)
   // plugin->AddDataFile("tag.xml");
   // plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
// Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir(Form("work2/%s",filename));
   //plugin->SetGridWorkingDir("work");
// Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); 

// Declare the analysis source files names separated by blancs. To be compiled runtime
// using ACLiC on the worker nodes.
   plugin->SetAnalysisSource("AliAnalysisTaskCentralityTreeMaker.cxx");
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->SetAdditionalLibs("AliAnalysisTaskCentralityTreeMaker.h AliAnalysisTaskCentralityTreeMaker.cxx");

// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
   plugin->SetDefaultOutputs();
// Optionally define the files to be archived.
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");

// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("CentralityAnalysis.C");

// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(100);
// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable("CentralityAnalysis.sh");
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
   plugin->SetMaxInitFailed(5);
// Optionally resubmit threshold.
   plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(30000);
// Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
// Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("CentralityAnalysis.jdl");
// Optionally modify job price (default 1)
   plugin->SetPrice(1);      
// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   return plugin;
}
