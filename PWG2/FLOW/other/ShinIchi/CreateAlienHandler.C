// 115393 115401 115413 115414 115514
// 117048 117050 117052 117053 117054
// 117059 117060 117063 117065 117077
// 117082 117086 117092
AliAnalysisGrid* CreateAlienHandler()
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
//   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
// plugin->SetRunMode("test");
   plugin->SetRunMode("full");
// plugin->SetRunMode("terminate");
// Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
// plugin->SetROOTVersion("v5-26-00b");
   plugin->SetROOTVersion("v5-27-06b");
// plugin->SetAliROOTVersion("v4-19-19-AN");
// plugin->SetAliROOTVersion("v4-19-Rev-25");
// plugin->SetAliROOTVersion("v4-21-05-AN");
   plugin->SetAliROOTVersion("v4-21-06-AN");
// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
//   plugin->SetGridDataDir("/alice/sim/LHC10a6");
   // On real reconstructed data:
   plugin->SetGridDataDir("/alice/data/2010/LHC10h");
// Set data search pattern
   // plugin->SetDataPattern("*ESDs.root");
// Data pattern for reconstructed data
// plugin->SetDataPattern("*ESDs/pass1_plusplusplus/*ESDs.root");
   plugin->SetDataPattern("*ESDs/pass1_4plus/*ESDs.root");
   plugin->SetRunPrefix("000");   // real data
// ...then add run numbers to be considered
//   plugin->AddRunNumber(125020);
// plugin->AddRunNumber(137161); // 137161 137162 real data
   plugin->AddRunNumber(137609); // 137161 *137431 *137549 *137595 *137609 137638 137639 137693
//   plugin->SetOutputSingleFolder("output");
//   plugin->SetOutputToRunNo();
// Method 2: Declare existing data files (raw collections, xml collections, root file)
// If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
// XML collections added via this method can be combined with the first method if
// the content is compatible (using or not tags)
//   plugin->AddDataFile("tag.xml");
//   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
// Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("work_137609");
// Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
// Declare the analysis source files names separated by blancs. To be compiled runtime
// using ACLiC on the worker nodes.
   plugin->SetAnalysisSource("AliAnalysisTaskPt.cxx");
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->SetAdditionalLibs("AliAnalysisTaskPt.h AliAnalysisTaskPt.cxx");
// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
// plugin->SetDefaultOutputs(kTRUE);
   plugin->SetDefaultOutputs(kFALSE);
   plugin->SetOutputFiles("EveMyTree.root");
// Optionally define the files to be archived.
// plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
// plugin->SetOutputArchive("log_archive.zip:mydata.dat");
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("TaskPt.C");
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(100);
// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutableCommand("aliroot -b -q");
   plugin->SetExecutable("TaskPt.sh");
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
// Optionally resubmit threshold.
//   plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(30000);
// Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
// Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("TaskPt.jdl");
// Optionally modify job price (default 1)
   plugin->SetPrice(1);      
// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   return plugin;
}
