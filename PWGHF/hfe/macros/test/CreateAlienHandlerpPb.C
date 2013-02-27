//Macro to test Analysis Macros on the GRID 
//please check settings for output files
//for local test use 'test' mode

AliAnalysisGrid* CreateAlienHandlerpPb(bool isAOD = kFALSE)
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
//   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   plugin->SetOverwriteMode();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode("test");
//   plugin->SetRunMode("offline");
//   plugin->SetRunMode("full");
//   plugin->SetRunMode("terminate");
   plugin->SetNtestFiles(1);
// Set versions of used packages

   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-34-02-1");
   plugin->SetAliROOTVersion("v5-04-32-AN");
// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
// On real reconstructed data:
   plugin->SetGridDataDir("/alice/data/2013/LHC13b");
// Set data search pattern
   if(isAOD)
      plugin->SetDataPattern("*/pass2/AOD/*/AliAOD.root");
   else
      plugin->SetDataPattern("*/pass2/*/AliESDs.root");

//same for pp MC:
//   plugin->SetGridDataDir("/alice/sim/LHC10f6a");
//  plugin->SetDataPattern("*/*/AliESDs.root");
// Data pattern for reconstructed data
//   plugin->SetDataPattern("*AliAOD.root"); //esta linea sirve para pruebas

   plugin->SetRunPrefix("000");   // real data

// ...then add run numbers to be considered
//   plugin->SetRunRange(122374,126437); //sim data
//10d
// plugin->AddRunNumber(126437); //sim data
//11h.pass2
plugin->AddRunNumber(195351); 
//   plugin->SetOutputSingleFolder("output");
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
   plugin->AddIncludePath("-I. .I$ALIEN_ROOT/api/lib -I$ROOTSYS/lib -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/PWGHF/ -I$ALICE_ROOT/PWGHF/hfe/macros -I$ALICE_ROOT/PWGHF/hfe -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/ANALYSIS/Tender -I$ALICE_ROOT/ANALYSIS/TenderSupplies -I$ALICE_ROOT/PWG/ -I$ALICE_ROOT/PWG/FLOW -I$ALICE_ROOT/PWG/Base -I$ALICE_ROOT/PWG/Tasks");
 //  plugin->SetAdditionalLibs("libGui.so  libXMLParser.so libSTEERBase.so libESD.so libAOD.so libCDB.so libANALYSIS.so libANALYSISalice.so libCORRFW.so  libTENDER.so libProof.so libRAWDatabase.so libSTEER.so libTOFbase.so libPWGDQdielectron.so");// ConfigLowMassDiE.C")
   plugin->SetAdditionalLibs("libGui.so  libXMLParser.so libSTEERBase.so libESD.so libAOD.so libCDB.so libANALYSIS.so libANALYSISalice.so libCORRFW.so  libPWGflowBase.so libPWGflowTasks.so libPWGHFhfe.so libTENDER.so libProof.so libRAWDatabase.so libSTEER.so libTOFbase.so");
   
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
//   plugin->SetAdditionalLibs("AliAnalysisHelperJetTasks.h AliAnalysisHelperJetTasks.cxx AliAnalysisTaskPartonDisc.h AliAnalysisTaskPartonDisc.cxx");
// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
//   plugin->SetOutputFiles("Output.root");
   //plugin->SetDefaultOutputs(); 
   plugin->SetDefaultOutputs(kFALSE);
   plugin->SetOutputFiles("AnalysisResults.root"); 
     //plugin->SetOutputFiles("cbaumann_LMEEpp2010_out.root"); 
//   plugin->SetOutputFiles("cbaumann_lowmass_out.root cbaumann_lowmass_CF.root");
// Optionally define the files to be archived.
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@disk=2 root_archive.zip:*.root@disk=2");
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("TPCTOFanalysispPb.C");
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
//   plugin->SetSplitMaxInputFileNumber(2);
// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable("TPCTOFanalysispPb.sh");
   plugin->SetExecutableCommand("aliroot -b -q");
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
// Optionally resubmit threshold.
//   plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(30000);
// Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
// Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("TPCTOFanalysispPb.jdl");
// Optionally modify job price (default 1)
   plugin->SetPrice(1);      
// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   return plugin;
}
