
AliAnalysisGrid* CreateAlienHandler(const char * runlist, TList * listCode, const char * mode, Bool_t isMC)
{
  // FIXME: different paths
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
//  if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien(); 
  plugin->SetOverwriteMode();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(mode);
  // Set versions of used packages
  // FIXME: PAR FILES OPTIONAL?
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-28-00a");
    plugin->SetAliROOTVersion("v4-21-20-AN");
// PAR files: I'm using a modified ANALYSISalice package, so I need to load par files for everything
/*   plugin->EnablePackage("STEERBase");
   plugin->EnablePackage("ESD");
   plugin->EnablePackage("AOD");
   plugin->EnablePackage("CORRFW");
   plugin->EnablePackage("ANALYSIS");
   plugin->EnablePackage("ANALYSISalice");
   plugin->EnablePackage("OADB");
*/

  // Declare input data to be processed.
  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
   if(!isMC) 
     //plugin->SetGridDataDir("/alice/data/2010/LHC10h/");
     plugin->SetGridDataDir("/alice/sim/LHC11a10a/");
   else
     //plugin->SetGridDataDir("/alice/sim/LHC10h8/");
     plugin->SetGridDataDir("/alice/sim/LHC11a10a/");

   // Set data search pattern
   if(!isMC)
     //plugin->SetDataPattern("*/pass1/*AliESDs.root");
     //plugin->SetDataPattern("*/pass2/*AliESDs.root");
     plugin->SetDataPattern("*AliESDs.root");
   else
     plugin->SetDataPattern("AliESDs.root");

   plugin->SetNrunsPerMaster(1); // One run per master job

   // ...then add run numbers to be considered   
   plugin->AddRunNumber(runlist);


   //plugin->SetGridWorkingDir("LambdaK0/");
   plugin->SetGridWorkingDir("LambdaK0MC/");
   plugin->SetGridOutputDir("out");

   // plugin->SetDefaultOutputs(kFALSE);
   // plugin->SetOutputFiles("EventStat_temp.root lambdak0.root");
   // plugin->SetOutputArchive("log_archive.zip:stdout,stderr");

 
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.

   // Load additional stuff
   TString cxxToLoad, hToLoad;
   TIterator * iter = listToLoad->MakeIterator();
   TObjString * str = 0;
   Bool_t first = kTRUE;
   while (str = (TObjString*) iter->Next()) {
     cxxToLoad = cxxToLoad + (first ? "" : " ") + str->String();
     str->String().ReplaceAll("cxx","h");
     hToLoad = hToLoad + (first ? "" : " ") + str->String();
     first = kFALSE;
   }
   plugin->SetAnalysisSource(cxxToLoad.Data());
   plugin->SetAdditionalLibs((hToLoad + " " + cxxToLoad).Data());



  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("AnalysisStrange.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  //plugin->SetSplitMaxInputFileNumber(100);
  plugin->SetSplitMaxInputFileNumber(50);
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  plugin->SetMaxInitFailed(5);
  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(90000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("TaskStrange.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");
  
  return plugin;
}
