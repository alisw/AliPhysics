//Macro to test Analysis Macros on the GRID 
//please check settings for output files
//for local test use 'test' mode
/*#include <fstream>
using namespace std;

#include <TString.h>

#include "AliAnalysisAlien.h"

void AddRunNumbers(AliAnalysisAlien* plugin, const Char_t* filename);
*/

AliAnalysisGrid* CreateAlienHandlerPbPb(const Char_t* inputRunList)
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
//   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   plugin->SetOverwriteMode();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
//   plugin->SetRunMode("test");
//   plugin->SetRunMode("offline");

   plugin->SetRunMode("full");

//   plugin->SetRunMode("terminate");
   plugin->SetNtestFiles(1);
// Set versions of used packages

   plugin->SetAPIVersion("V1.1x");
   //   plugin->SetROOTVersion("v5-33-02b");
   plugin->SetROOTVersion("v5-34-30-alice5-1");
   plugin->SetAliROOTVersion("v5-08-15-1");
   plugin->SetAliPhysicsVersion("vAN-20160821-1");
// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
// On real reconstructed data:
   plugin->SetGridDataDir("/alice/data/2011/LHC11h_2");
// Set data search pattern
   plugin->SetDataPattern("*/pass2/*/AliESDs.root");

//same for pp MC:
//   plugin->SetGridDataDir("/alice/sim/LHC10f6a");
//  plugin->SetDataPattern("*/*/AliESDs.root");
// Data pattern for reconstructed data
//   plugin->SetDataPattern("*AliAOD.root"); //esta linea sirve para pruebas
   
   plugin->SetRunPrefix("000");   // real data
   AddRunNumbers(plugin, inputRunList);
   
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
   plugin->SetGridOutputDir("outputDst");
// Declare the analysis source files names separated by blancs. To be compiled runtime
// using ACLiC on the worker nodes.
   //plugin->SetAnalysisSource("alien:///alice/cern.ch/user/i/iarsene/work/AliReducedEvent.cxx alien:///alice/cern.ch/user/i/iarsene/work/AliAnalysisTaskReducedTree.cxx");

   plugin->AddIncludePath("-I. -I$ALIEN_ROOT/api/lib -I$ROOTSYS/lib -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
   //plugin->SetAdditionalLibs("");
   
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   //plugin->SetAdditionalLibs("AliReducedEvent.cxx");


// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
//   plugin->SetOutputFiles("Output.root");
   //plugin->SetDefaultOutputs(); 
   //or specify files:
   plugin->SetDefaultOutputs(kFALSE);
   plugin->SetOutputFiles("dstTree.root");    // TODO: adapt depending on what is attached to the train
     
// Optionally define the files to be archived.
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@disk=2 root_archive.zip:*.root@disk=2");
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("ReducedEventAnalysis.C");
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(10);
// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable("ReducedEventAnalysis.sh");
   plugin->SetExecutableCommand("aliroot -b -q");
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
// Optionally resubmit threshold.
//   plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(20000);
// Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
// Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("ReducedEventAnalysis.jdl");
// Optionally modify job price (default 1)
   plugin->SetPrice(1);      
// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   return plugin;
}

//_______________________________________________________________________________
void AddRunNumbers(AliAnalysisAlien* plugin, const Char_t* filename) {
   //
   // Add all the runs listed in filename to the plugin
   //
   ifstream in;
   in.open(filename);
   
   TString line;
   //loop over file                                                                                                                                                             
   while(in.good()) {
      in >> line;
      if (!line.IsNull()) {
         cout << "Adding run: " << line <<endl;
         plugin->AddRunNumber(line.Atoi());
      }
   }
}

