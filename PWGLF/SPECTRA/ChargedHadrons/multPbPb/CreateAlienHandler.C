
AliAnalysisGrid* CreateAlienHandler(TString dataset, TList * listToLoad, const char * mode = "full", Bool_t isMC = 0)
{
  

// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
//   if (!AliAnalysisGrid::CreateToken()) return NULL;
  TGrid::Connect("alien:");

   AliAnalysisAlien *plugin = new AliAnalysisAlien();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(mode);
   //plugin->SetRunMode("test");
// Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-27-06-2");
   plugin->SetAliROOTVersion("v4-21-03-AN");
// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
// LHC09d
// /alice/data/2009/LHC09d/000104892/ESDs/pass6/
   plugin->SetGridDataDir(dataset);
// Set data search pattern
   plugin->SetDataPattern("AliESDs.root");
// ...then add run numbers to be considered
//   plugin->AddRunNumber(104892);
//   plugin->AddRunNumber(300001);
// Method 2: Declare existing data files (raw collections, xml collections, root file)
// If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
// XML collections added via this method can be combined with the first method if
// the content is compatible (using or not tags)
   // plugin->AddDataFile("tag.xml");
   // plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
// Define alien work directory where all files will be copied. Relative to alien $HOME.
   TString output = "MultPb/Output_";
   dataset.ReplaceAll("/","_");
   output += dataset;
   plugin->SetGridWorkingDir(output.Data());
// Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("out"); 
// Declare the analysis source files names separated by blancs. To be compiled runtime
// using ACLiC on the worker nodes.
   TIterator * iter = listToLoad->MakeIterator();
   TObjString * name = 0;
   TString sources = "";
   // while (name = (TObjString *)iter->Next()) {
   //   gSystem->ExpandPathName(name->String());
   //   name->String().ReplaceAll("+","");
   //   sources = sources + name->String() + " ";
   // }
   while (name = (TObjString *)iter->Next()) {
     gSystem->ExpandPathName(name->String());
     name->String().ReplaceAll("+","");     
     sources = sources + gSystem->BaseName(name->String().Data()) + " ";
   }
   plugin->SetAnalysisSource(sources.Data());
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
//   plugin->SetAdditionalLibs("AliCollisionsNormalization.cxx AliCollisionNormalizationTask.cxx AliPhysicsSelection.cxx AliCollisionsNormalization.h AliCollisionNormalizationTask.h AliPhysicsSelection.h");
   // iter = listToLoad->MakeIterator();
   // name = 0;
   // sources = "";
   // while (name = (TObjString *)iter->Next()) {
   //   gSystem->ExpandPathName(name->String());
   //   name->String().ReplaceAll("+","");     
   //   sources = sources + gSystem->BaseName(name->String().Data()) + " ";
   // }
   // while (name = (TObjString *)iter->Next()) {
   //   gSystem->ExpandPathName(name->String());
   //   name->String().ReplaceAll("+","");
   //   sources = sources + name->String() + " ";
   // }
   plugin->SetAdditionalLibs(sources.Data());
   // I'm using a modified ANALYSISalice package, so I need to load par files for everything
   // plugin->EnablePackage("STEERBase");
   // plugin->EnablePackage("ESD");
   // plugin->EnablePackage("AOD");
   // plugin->EnablePackage("CORRFW");
   // plugin->EnablePackage("ANALYSIS");
   // plugin->EnablePackage("ANALYSISalice");

// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
//   plugin->SetDefaultOutputs(kFALSE);
   //   plugin->SetOutputFiles(Form("EventStat_temp.root %s",outfilename));
   //   plugin->SetOutputFiles("EventStat_temp.root multPbPbtracks.root");
// Optionally define the files to be archived.
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("AnalysisMultPb.C");
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(100);
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
   plugin->SetMaxInitFailed(5);
// Optionally resubmit threshold.
   plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(30000);
// Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
// Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("TaskMultPb.jdl");
// Optionally modify job price (default 1)
   plugin->SetPrice(1);      
// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   return plugin;
}
