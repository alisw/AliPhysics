AliAnalysisGrid* CreateAlienHandler()
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
  //   if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  //plugin->SetRunMode("test");
  plugin->SetRunMode("terminate");
  //plugin->SetRunMode("full");
  //plugin->SetFileForTestMode("TestFiles.txt");
  // Set versions of used packages
  //plugin->SetAPIVersion("V1.1x");

  //  plugin->SetROOTVersion("v5-34-08-6");
  //  plugin->SetAliROOTVersion("v5-06-30");

  //  plugin->SetROOTVersion("v5-34-26");
  //plugin->SetAliROOTVersion("v5-06-Rev-14");
  //plugin->SetAliPhysicsVersion("vAN-20171130-1");
  //plugin->SetAliPhysicsVersion("vAN-20180115-1");
  plugin->SetAliPhysicsVersion("vAN-20180314-1"); 
  
  plugin->SetNtestFiles(1);
  // Declare input data to be processed.
  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  plugin->SetGridDataDir("/alice/sim/2017//LHC17j7_ZDCfix_extra"); //"/alice/sim/2017//LHC17j7_ZDCfix_extra" for MC
                                                                   //"/alice/data/2017/LHC17n/" for data 
  // Set data search pattern
  plugin->SetDataPattern("*ESDs.root"); // for running on MC
  // Data pattern for reconstructed data
  //plugin->SetDataPattern("/pass1/*ESDs.root");
  //plugin->SetRunPrefix("000");   // real data
  // ...then add run numbers to be considered
  //
  // REMOVE THESE NUMBERS IF ONLY GOLDEN RUN NEEDS TO BE PROCESSED
  //

  plugin->AddRunNumber(280234);
  plugin->AddRunNumber(280235);
  //TString date = gSystem->GetFromPipe("date +%F");
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  //plugin->SetGridWorkingDir(Form("XeXeTpcPikPSpectra_%s/", date.Data()));
  plugin->SetGridWorkingDir("XeXe_MC_2018_03_18/");

  // Declare alien output directory. Relative to working directory.
  //plugin->SetGridOutputDir(Form("XeXeTpcPikPSpectraOutput_%sdata", date.Data())); // In this case will be $HOME/work/output
  plugin->SetGridOutputDir("Output_2018_03_18_MC");
  
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  plugin->SetAnalysisSource("AliAnalysisTaskTpcSpectra.cxx");
  
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalLibs("libOADB.so AliAnalysisTaskTpcSpectra.h AliAnalysisTaskTpcSpectra.cxx libpythia6_4_28.so");
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  //   plugin->SetOutputFiles("Pt.ESD.1.root");
  plugin->SetDefaultOutputs();
  // Optionally define the files to be archived.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  //  plugin->SetAnalysisMacro("TaskSpectra.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(200); //10 for Data
  // Optionally modify the executable name (default analysis.sh)
  //  plugin->SetExecutable("TaskPt.sh");
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  //   plugin->SetMaxInitFailed(5);
  // Optionally resubmit threshold.
  //   plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("TaskSpectra.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");
  
  // merging detauls
  plugin->SetMergeViaJDL(kTRUE);
  // plugin->SetMergeViaJDL(kFALSE);
  
  plugin->SetMaxMergeFiles(10000);
  plugin->SetMaxMergeStages(2); //this number you can change depending how many files you have
  
  return plugin;
}
