
AliAnalysisGrid* CreateAlienHandler()
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
//   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   plugin->SetOverwriteMode();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   //  plugin->SetRunMode("terminate"); //RR
   plugin->SetRunMode("full"); // RR
   // plugin->SetRunMode("test"); // RR
   // Set versions of used packages   
   plugin->SetAPIVersion("V1.1x");
   plugin->SetAliROOTVersion("vAN-20140619");
   plugin->SetROOTVersion("v5-34-08-6");
   //   plugin->SetAliROOTVersion("vAN-20150513");
   //   plugin->SetAliROOTVersion("v5-06-17-30");
   // The aliroot versions compatible with Yosemite with ROOT versions v5-34-30 are:
   //   VO_ALICE@AliRoot::v5-06-15-30
   //   VO_ALICE@AliRoot::v5-06-17-30
   //   plugin->SetROOTVersion("v5-34-30");
   plugin->SetNtestFiles(2);

// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
   // For DATA
//   plugin->SetGridDataDir("/alice/data/2010/LHC10h"); //RR
   // For MonteCarlo
   plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis"); //RR
   //Set data search pattern
   //   plugin->SetDataPattern("*ESDs/pass2/AOD160/*AOD.root");
   //SetData Search pattern for MonteCarlo
   plugin->SetDataPattern("*AOD162/*AOD.root");
   //   plugin->SetDataPattern("*ESDs/pass2/AOD156/*AOD.root");
   //   plugin->SetDataPattern("*ESDs/pass2/AOD160/*AOD.root");
   //   plugin->SetDataPattern("*ESDs.root");
// Data pattern for reconstructed data
//   plugin->SetDataPattern("*ESDs/pass4/*ESDs.root");
//   plugin->SetRunPrefix("000");   // real data
// ...then add run numbers to be considered

   //   plugin->AddRunNumber(137161); // Golden RUN

   //   plugin->AddRunNumber(139438);

   //139510, 139507, 139505, 139503, 139465,
   //   plugin->AddRunNumber(139507);
   //   137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595,
   //     137549, 137546, 137544, 137541,137539, 137531, 137530, 
   //   138396, 138364, 138275, 138225, 138201,138197, 138192, 138190, 137848, 137844,
   //     137752, 137751, 137724,137722, 137718, 137704,
   
   //   139105,138872, 138871, 138870,
   //138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638,
   //138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469,
   //138442, 138439, 138438, 
   plugin->AddRunNumber(139105);
   plugin->AddRunNumber(138872);
   plugin->AddRunNumber(138871);
   plugin->AddRunNumber(138870);
   plugin->AddRunNumber(138837);
   plugin->AddRunNumber(138732);
   plugin->AddRunNumber(138730);
   plugin->AddRunNumber(138666);
   plugin->AddRunNumber(138662);
   plugin->AddRunNumber(138653);
   plugin->AddRunNumber(138652);
   plugin->AddRunNumber(138638);
   plugin->AddRunNumber(138624);
   plugin->AddRunNumber(138621);
   plugin->AddRunNumber(138583);
   plugin->AddRunNumber(138582);
   plugin->AddRunNumber(138579);
   plugin->AddRunNumber(138578);
   plugin->AddRunNumber(138534);
   plugin->AddRunNumber(138469);
   plugin->AddRunNumber(138442);
   plugin->AddRunNumber(138439);
   plugin->AddRunNumber(138438);

	 
//   plugin->SetOutputSingleFolder("output");
   plugin->SetOutputToRunNo();

// Method 2: Declare existing data files (raw collections, xml collections, root file)
// If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
// XML collections added via this method can be combined with the first method if
// the content is compatible (using or not tags)
//   plugin->AddDataFile("tag.xml");
//   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
   
// Define alien work directory where all files will be copied. Relative to alien $HOME.
   // plugin->SetGridWorkingDir("workDATA/20150520a");
    plugin->SetGridWorkingDir("workMCData/20150528a");
// Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
// Declare the analysis source files names separzated by blancs. To be compiled runtime
// using ACLiC on the worker nodes.
   plugin->SetAnalysisSource("AliAnalysisFBMultFluct.cxx");
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->SetAdditionalLibs("AliAnalysisFBMultFluct.h AliAnalysisFBMultFluct.cxx");
// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
   //  plugin->SetOutputFiles("Pt.ESD.1.root");
   plugin->SetDefaultOutputs(kFALSE);
   plugin->SetOutputFiles("AOD.DATA.root");
// Optionally define the files to be archived.
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("TaskFBMultFromAOD.C");
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(100);
// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable("TaskFBMultFromAOD.sh");
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
   plugin->SetFileForTestMode("files.txt");
// Optionally resubmit threshold.
//   plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(30000);
// Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
// Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("TaskFBMultFromAOD.jdl");
// Optionally modify job price (default 1)
   plugin->SetPrice(1);      
// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   
//Merging
   plugin->SetMergeViaJDL(kTRUE);
   return plugin;
}
