
AliAnalysisGrid* CreateAlienHandlerAODpp()
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
//   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   plugin->SetOverwriteMode();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode("test");
//    plugin->SetRunMode("full");
// plugin->SetRunMode("terminate");
//   plugin->SetRunMode("offline");
//   plugin->SetRunMode("submit");
// Set versions of used packages
   plugin->SetNtestFiles(1);

   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-33-02b");
   plugin->SetAliROOTVersion("v5-03-23-AN");
// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
//   plugin->SetGridDataDir("/alice/cern.ch/user/k/kleinb/analysisESD/LHC10b/output_pwg4train_LHC10b_pass2_100910a");
// On real reconstructed data:
   plugin->SetGridDataDir("/alice/data/2010/LHC10d/");
   // Set data search pattern
 //plugin->SetDataPattern("*/pass2/AOD043/*/AliAOD.root");
 //plugin->SetDataPattern("*/pass2/AOD039/*/AliAOD.root");
 plugin->SetDataPattern("*/pass2/AOD057/*/AliAOD.root");
// Data pattern for reconstructed data
//   plugin->SetDataPattern("*AliAOD.root"); //esta linea sirve para pruebas
   plugin->SetRunPrefix("000");   // real data
// ...then add run numbers to be considered

//Use Evas RCT-Good-RUN-list splitter:
   //LHC10eAOD43
   //TString runs="130850 130848 130844 130834 130803 130799 130795 130793 130704 130696 130628 130623 130621 130601 130526 130524 130520 130519 130517 130481 130480 130375 130369 130358 130356 130343 130178 130158 130157 130151 130149 130148 129983 129966 129962 129961 129738 129736 129735 129729 129667 129654 129641 129540 129536 129520 129042 128912 128855 128853 128850 128843 128836 128833 128824 128820 128777 128621 128615 128609 128596 128592 128504 128503 128498 128486 128452 128366 128260 128257 128189 128186 128185 128182 128180 128175 127942 127941 127937 127936 127935 127933 127932 127931 127822 127817 127815 127814 127729 127724 127719 127718 127714";	
  // TString runs="130848 130844 130834 130803 130799 130795 130793 130704 130696 130628 130623 130621 130601 130526 130524 130520 130519 130517 130481 130480 130375 130369 130358 130356 130343 130178 130158 130157 130151 130149 130148 129983 129966 129962 129961 129738 129736 129735 129729 129667 129654 129641 129540 129536 129520 129042 128912 128855 128853 128850 128843 128836 128833 128824 128820 128777 128621 128615 128609 128596 128592 128504 128503 128498 128486 128452 128366 128260 128257 128189 128186 128185 128182 128180 128175 127942 127941 127937 127936 127935 127933 127932 127931 127822 127817 127815 127814 127729 127724 127719 127718 127714";	
   
   //LHC10dAOD39:
   //TString runs="126437 126432 126425 126424 126422 126409 126408 126407 126406 126404 126359 126352 126351 126284 126283 126168 126160 126158 126097 126090 126088 126082 126081 126078 126073 126008 126007 125855 125851 125850 125849 125842 125632 125630 125628 125186 125140 125139 125134 125101 125097 125085 122374";
   //LHC10dAOD57:
   //TString runs="122374 124358 124360 124371 124374 124378 124380 124381 124383 124385 124388 124600 124606 125085 125097 125101 125134 125139 125140 125186 125295 125632 125842 125849 125850 125851 125855 126007 126073 126078 126081 126082 126088 126090 126097 126158 126160 126167 126168 126177 126283 126284 126351 126352 126359 126404 126406 126407 126408 126409 126422 126424 126425 126432 126437";
  //   TString runs="       124358 124360 124371 124374 124378 124380 124381 124383 124385        124600 124606 125085 125097 125101 125134 125139 125140 125186 125295 125632 125842 125849 125850 125851 125855 126007 126073 126078 126081 126082 126088 126090 126097 126158 126160 126167 126168 126177 126283 126284 126351 126352 126359 126404 126406 126407 126408 126409 126422        126425 126432 126437";
   //Taken out: 125296 
/*
 TObjArray* array = runs.Tokenize ( " " );
   TObjString *str;
   TString strr;
   for ( Int_t i = 0;i < array->GetEntriesFast();i++ ) {
	  str = ( TObjString * ) array->At ( i );
	  strr = str->GetString();
	  if ( !strr.IsNull() ) {
		 plugin->AddRunNumber(strr.Atoi());
	  }
   }  

*/

   
   plugin->AddRunNumber(126432); //sim data
  
   plugin->SetOutputToRunNo();
// Method 2: Declare existing data files (raw collections, xml collections, root file)
// If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
// XML collections added via this method can be combined with the first method if
// the content is compatible (using or not tags)
//   plugin->AddDataFile("tag.xml");
//   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
// Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("work");
// Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
// Declare the analysis source files names separated by blancs. To be compiled runtime
// using ACLiC on the worker nodes.
//   plugin->SetAnalysisSource("AliAnalysisHelperJetTasks.cxx AliAnalysisTaskPartonDisc.cxx");
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/Tender -I$ALICE_ROOT/TenderSupplies -I$ALICE_ROOT/PWGHF/ -I$ALICE_ROOT/PWGHF/hfe/macros -I$ALICE_ROOT/PWGHF/hfe -I$ALICE_ROOT/PWG/ -I$ALICE_ROOT/PWG/FLOW -I$ALICE_ROOT/PWG/Base -I$ALICE_ROOT/PWG/Tasks");
   plugin->SetAdditionalLibs("libGui.so  libXMLParser.so libSTEERBase.so libAOD.so libCDB.so libANALYSIS.so libANALYSISalice.so libCORRFW.so  libTender.so libProof.so libRAWDatabase.so libSTEER.so libTOFbase.so libPWGflowBase.so libPWGflowTasks.so libPWGHFhfe.so");// ConfigLowMassDiE.C")
   
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
//   plugin->SetAdditionalLibs("AliAnalysisHelperJetTasks.h AliAnalysisHelperJetTasks.cxx AliAnalysisTaskPartonDisc.h AliAnalysisTaskPartonDisc.cxx");
// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
//   plugin->SetOutputFiles("Output.root");

   // plugin->SetDefaultOutputs(); 
   plugin->SetDefaultOutputs(kFALSE);
   plugin->SetOutputFiles("rbailhac_tpctofppaod_out.root");


// Optionally define the files to be archived.
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@disk=2 root_archive.zip:*.root@disk=2");
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("Tpctofaodpp.C");
// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable("Tpctofaodpp.sh");
   plugin->SetExecutableCommand("aliroot -b -q");
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
// Optionally resubmit threshold.
//   plugin->SetMasterResubmitThreshold(90);

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
//   plugin->SetJDLName("DielAnalysis.jdl");
// Optionally modify job price (default 1)
   plugin->SetPrice(1);      
// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");

// Enable Merging on the GRID:
//   plugin->SetMergeViaJDL(kTRUE);



   return plugin;
}
