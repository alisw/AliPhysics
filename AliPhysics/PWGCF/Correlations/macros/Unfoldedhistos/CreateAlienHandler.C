#ifdef __ECLIPSE_IDE
//  few includes and external declarations just for the IDE
#include "AliAnalysisAlien.h"
#endif // ifdef __ECLIPSE_IDE declaration and includes for the ECLIPSE IDE

#include "runCorrelationsStudiesConfigMacro.H"

AliAnalysisGrid* CreateAlienHandler(const char *runMode,Bool_t gridMerge)
{
    // Check if user has a valid token, otherwise make one. This has limitations.
    // One can always follow the standard procedure of calling alien-token-init then
    //   source /tmp/gclient_env_$UID in the current shell.

  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  //Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(runMode);

  // Set MC chain if MC truth as input
  plugin->SetUseMCchain(bMC && bMConlyTruth);

  plugin->SetNtestFiles(nNoOfTestFiles); // num of test files in "test" mode

  if (TString(runMode).EqualTo("test"))
    plugin->SetCheckCopy(kFALSE);


  // Set versions of used packages
  plugin->SetAliPhysicsVersion(szAliPhysicsVersion.Data());

  plugin->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  /////////////////////////////////////////////////////////////////////////
  plugin->SetDataPattern(szDataPattern.Data());
  plugin->SetGridDataDir(szDataDir.Data()); // Data

  /* CHECK: for MC seems not needed */
  if (!szRunPrefix.IsWhitespace())
    plugin->SetRunPrefix(szRunPrefix.Data());

  /* run numbers to analyse */
  for (Int_t i=0;i<listOfActiveRuns.GetEntriesFast();i++) {
    plugin->AddRunNumber(((TObjString*) listOfActiveRuns.At(i))->GetString().Data());
  }

  /* alternatively provide run number */
//  plugin->AddRunNumber( 139505 );
// Alternatively use run range
//  plugin->SetRunRange(138653, 138666);

  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir(szGridWorkingDir.Data());  // NOTE: Change name here every new run!!!eclare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  plugin->SetOutputToRunNo(); // we want the run number as output subdirectory
  plugin->SetDefaultOutputs(kTRUE);

//  plugin->SetMergeExcludes("Viscosity.root EventStat_temp.root");
  plugin->SetMergeViaJDL(gridMerge);
  plugin->SetMaxMergeFiles(30);

  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  /* plugin->SetAnalysisSource("AliCSAnalysisCutsBase.cxx "
      "AliCSTrackCutsBase.cxx "
      "AliCSTrackMaps.cxx "
      "AliCSEventCuts.cxx "
      "AliCSTrackCuts.cxx "
      "AliCSPIDCuts.cxx "
      "AliCSTrackSelection.cxx "
      "AliDptDptCorrelations.cxx "
      "AliCSPairAnalysis.cxx "
      "AliAnalysisTaskCorrelationsStudies.cxx"); */

  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  /* plugin->SetAdditionalLibs("AliCSAnalysisCutsBase.cxx AliCSAnalysisCutsBase.h "
      "AliCSTrackCutsBase.cxx AliCSTrackCutsBase.h "
      "AliCSTrackMaps.cxx AliCSTrackMaps.h "
      "AliCSEventCuts.cxx AliCSEventCuts.h "
      "AliCSTrackCuts.cxx AliCSTrackCuts.h "
      "AliCSPIDCuts.cxx AliCSPIDCuts.h "
      "AliCSTrackSelection.cxx AliCSTrackSelection.h "
      "AliDptDptCorrelations.cxx AliDptDptCorrelations.h "
      "AliCSPairAnalysis.cxx AliCSPairAnalysis.h "
      "AliAnalysisTaskCorrelationsStudies.h AliAnalysisTaskCorrelationsStudies.cxx"); */

  // alternatively pass a par file with source information
//  plugin->EnablePackage("PWGCFCorrelationsDPhi.par");

// Declare the output file names separated by blanks.
// (can be like: file.root or file.root@ALICE::Niham::File)
  plugin->SetOverwriteMode(kFALSE);

// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  if (nNoOfInputFiles != 0)
    plugin->SetSplitMaxInputFileNumber(nNoOfInputFiles);
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(15);
// Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
// Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  return plugin;
}
