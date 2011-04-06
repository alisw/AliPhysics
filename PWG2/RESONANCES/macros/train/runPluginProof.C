//
// This is an example steering macro for running RSN analysis task
// with the AliEn plugin on PROOF.
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
// In principle, the user should never modify this macro.
//
void runPluginProof
(
   const char *runMode         = "full",
   Int_t       nmix            = 0,
   const char *dataSet         = "/alice/data/LHC10h_000137162#esdTree",
   const char *clusterName     = "alice-caf.cern.ch",
   const char *testFile        = "pbpb_data.txt",
   const char *runOptions      = "esd_data_phys_cent",
   const char *analysisOptions = "tpcpid_tofpid",
   const char *outName         = "test_proof.root",
   const char *taskList        = "AddRsnAnalysisTask.C",
   //const char *taskPath        = "$(ALICE_ROOT)/PWG2/RESONANCES/macros/train"
   const char *taskPath        = "$(HOME)/code/resonances/alice-rsn-package/PWG2resonances/RESONANCES/macros/test/pulvir"
)
{
   //
   // === PREPARATION ==============================================================================
   //
   
   // this option is not needed when using plugin
   // gEnv->SetValue("XSec.GSI.DelegProxy","2");
   
   // some options
   TString opt(runOptions);
   opt.ToUpper();
   Bool_t useTender = opt.Contains("TENDER");
   Bool_t isMC      = opt.Contains("MC") || (!opt.Contains("DATA"));
   ::Info("runPlugin.C", "useTender = %d", useTender);
   ::Info("runPlugin.C", "isMC      = %d", isMC     );
   ::Info("runPlugin.C", "runOpts   = %s", runOptions);
   ::Info("runPlugin.C", "anaOpts   = %s", analysisOptions);
   
   // basic configurations
   gROOT->LoadMacro(Form("%s/AnalysisSetup.C", taskPath));
   AnalysisSetup(isMC, nmix, runOptions, outName, taskPath);

   //
   // === PLUGIN CONFIGURATION =====================================================================
   //

   // load and execute plugin configuration macro
   // pass to the macro, as FIRST argument, the common name
   // which is used for the output, since it must be equal
   // to the one defined here for the common output (for merging)
   gROOT->LoadMacro("PluginByRunProof.C");

   // build the list of arguments (for them see the 'PluginByRun.C' code)
   TString args("PluginByRunProof(");
   args += Form("\"%s\", ", dataSet);
   args += Form("\"%s\", ", testFile);
   args += Form("\"%s\") ", clusterName);

   // create the plugin
   AliAnalysisAlien *plugin = (AliAnalysisAlien*)gROOT->ProcessLine(args.Data());

   // set run mode
   plugin->SetRunMode(runMode);
   
   // add to manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) return;
   mgr->SetGridHandler(plugin);
   
   //
   // === ANALYSIS EXECUTION =======================================================================
   //

   gROOT->LoadMacro(Form("%s/AddRsnAnalysisTask.C", taskPath));
   AddRsnAnalysisTask(isMC, kFALSE, runOptions, analysisOptions);
   
   // initialize and start analysis
   if (!mgr->InitAnalysis()) {
      ::Error("runPlugin.C", "Failed to init analysis");
      return;
   }
   mgr->PrintStatus();
   mgr->StartAnalysis("proof");
}
