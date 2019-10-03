AliAnalysisGrid* CreateAlienHandler()
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
  //   if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();

  //========================================
  // SET RUN MODE :  "full",   "est",   "offline",   "ubmit" or "terminate"
  //========================================
  plugin->SetRunMode("test"); 

  plugin->SetNtestFiles(1); // Relevant only for run mode "test" 

  //========================================
  // Set versions of used packages
  //========================================
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-30-06-1");plugin->SetAliROOTVersion("v5-03-04-AN");

  //========================================
  // Declare input data to be processed.
  //========================================
  // /alice/data/2010/LHC10h/000137366/ESDs/pass1/10000137366005.2140/AliESDs.root
  plugin->SetGridDataDir("/alice/data/2011/LHC11h_2/");
  //plugin->SetDataPattern("/pass2/11000170308079.20/*AOD.root"); // TEST MODE

  plugin->SetDataPattern("*/pass2/*AOD.root");
  //plugin->SetDataPattern("*/pass2/*AOD.root");
  plugin->SetRunPrefix("000");   // real data
  plugin->AddRunNumber(170308);

  plugin->SetAdditionalLibs("libANALYSIS.so libANALYSISalice.so libEMCALUtils.so libPHOSUtils.so libGui.so libCDB.so libRAWDatabase.so libRAWDatarec.so libProof.so libSTEER.so libTOFbase.so libTOFrec.so libMinuit.so libRAWDatabase.so libRAWDatarec.so libAOD.so libCORRFW.so libPWGCFJCORRAN.so");

//  plugin->SetAdditionalLibs("PWG4JCORRAN.par libEMCALUtils.so libPHOSUtils.so");//AliCentralityBy1D_137161_GLAU.root AliCentralitySelectionTask.cxx");

  //========================================
  // Set Ouput Information  
  //========================================
  plugin->SetGridWorkingDir("PWG_CF/test_120320");
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  plugin->SetOutputToRunNo();
  plugin->SetDefaultOutputs(kFALSE);
  //plugin->SetPreferedSE("ALICE::NDGF::DCACHE");
  //plugin->SetOutputArchive("log_archive.zip:std*@disk=1 root_archive.zip:jcorran.root,AnalysisResults.root,EventStat_temp.root,*.stat@disk=2");
  plugin->SetOutputArchive("log_archive.zip:std*@disk=1 root_archive.zip:*.root,*stat,*.xml@disk=2");
  plugin->SetKeepLogs(kTRUE);
  plugin->SetOutputFiles("jcorran.root");
  plugin->SetTerminateFiles("event_stat.root");
  // plugin->SetOutputSingleFolder("output");
   
  //========================================
  // Optional
  //========================================
  //FK//   plugin->SetAnalysisSource("AliJCORRANTask.cxx");
  plugin->SetAnalysisMacro("TaskJC.C");
  plugin->SetSplitMaxInputFileNumber(20);
  plugin->SetExecutable("TaskJC.sh");
//  plugin->SetExecutableCommand("export MALLOC_CHECK_=0 ; root -b -q");
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  //   plugin->SetMaxInitFailed(5);
  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(95);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("TaskJC.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");


  return plugin;
}
