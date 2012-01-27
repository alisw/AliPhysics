//
// This is an example steering macro for running RSN analysis task
// with the AliEn plugin to launch a multiple analysis.
//
// Inputs:
//
//   - runMode     = AliEn plugin run mode
//   - suffix      = eventual suffix which is attached to the automatic name created for outputs
//   - partname    = a string which defines the resonance name (used for output)
//   - runList     = configuration file whith gives list of runs, pattern, prefix, options and AliEn path
//   - split       = corresponding JDL value
//   - nmerge      = number of outputs to be merged per stage
//
//   - taskList    = a string containin all the 'add-task' macros to be used
//   - taskPath    = the uniqu pathe where all 'add-task' macros (and all their needs) are stored
//
//   - workDirBase = path of the working directory (starting from AliEn home)
//   - outDir      = path of the output directory (w.r. to workDir)
//
// Notes:
//
//   - in case the source is an ESD, and if inputs are a MC production
//     the MC input handler is created by default
//
//
void runPlugin
(
   const char *runMode,
   const char *options,
   const char *alienDir,
   const char *addLibs,
   const char *addPars,
   
   const char *alienJobName,
   const char *alienRunPath,
   const char *alienRunPrefix,
   const char *alienRunPattern,
   const char *alienRunList,
   
   Int_t       nmix,
   Int_t       split,
   Int_t       nmerge,
   const char *macroPath,
   Int_t       runsPerMaster
)
{
   //
   // === PREPARATION ==============================================================================
   //
   
   // execute the general setup from the apposite macro
   // it returns also a TString value with the input tree name
   gROOT->LoadMacro("../AnalysisSetupRsnMini.C");
   TString out = Setup(nmix, options, "analysis.root", macroPath);
   if (out.Length() < 1) return;
   
   //
   // === PLUGIN CONFIGURATION =====================================================================
   //

   // load macro for plugin setup
   gSystem->AddIncludePath("-I$ALICE_ROOT/include");
   gROOT->LoadMacro("SetupPlugin.C++g");
   
   // define inputs
   PluginSetup::alienInputRuns = kTRUE;
   PluginSetup::nRunsPerMaster = runsPerMaster;
   PluginSetup::runPath        = alienRunPath;
   PluginSetup::runPrefix      = alienRunPrefix;
   PluginSetup::runPattern     = alienRunPattern;
   PluginSetup::runList        = alienRunList;
   PluginSetup::split          = split;
   PluginSetup::maxMergeFiles  = nmerge;

   // define output path
   PluginSetup::workDir  = "analysis/RSNpackage/mini/";
   PluginSetup::workDir += alienDir;
   
   // define common root for all generated files
   PluginSetup::jobName = alienJobName;
   
   // define additional libraries
   PluginSetup::addLibs = addLibs;
   PluginSetup::addPar  = addPars;
   
   // additional modalities
   PluginSetup::rootVersion = "v5-28-00d";
   PluginSetup::aliVersion = "v4-21-25-AN";
   
   // create the plugin (need to know if we want tender)
   if (!PluginSetup::CreatePlugin()) return;
   PluginSetup::plugin->SetRunMode(runMode);

   //
   // === ANALYSIS EXECUTION =======================================================================
   //
   
   // add plugin to analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) return;
   mgr->SetGridHandler(PluginSetup::plugin);
   
   // initialize and start analysis
   if (!mgr->InitAnalysis()) {
      ::Error("runPlugin.C", "Failed to init analysis");
      return;
   }
   mgr->PrintStatus();
   mgr->StartAnalysis("grid");
}
