void runPluginProof
(
   const char *runMode     = "full",
   
   Int_t       nmix        = 50,
   
   //const char *clusterName = "alice-caf.cern.ch",
   //const char *dataSet     = "/alice/sim/LHC10d1_000117112",
   //const char *options     = "esd_mc",
   //const char *dataSet     = "/alice/sim/LHC11a10b_000139507",
   //const char *options     = "esd_mc_pbpb",	
   //const char *dataSet     = "/PWG3/rbala/LHC11a10b_000137748_AOD048",
   //const char *options     = "aod_mc_pbpb",	
   
   const char *clusterName = "skaf.saske.sk",
   const char *dataSet     = "/alice/sim/LHC10d2_117099",
   const char *options     = "esd_mc",
   //const char *dataSet     = "/alice/sim/LHC11a10b_000139314_AOD048",
   //const char *options     = "aod_mc_pbpb",
   
   const char *outName     = "proof.root",
   const char *macroPath   = "..",
   const char *testFile    = "",
   const char *addLibs     = "libEventMixing.so PWG2resonances.par",
   const char *addPars     = ""
)
{
   //
   // === PREPARATION ==============================================================================
   //
   
   // this option is not needed when using plugin
   // gEnv->SetValue("XSec.GSI.DelegProxy","2");
   
   // execute the general setup from the apposite macro
   // it returns also a TString value with the input tree name
   gROOT->LoadMacro("../AnalysisSetupRsnMini.C");
   TString out = Setup(nmix, options, outName, macroPath);
   if (out.Length() < 1) return;

   //
   // === PLUGIN CONFIGURATION =====================================================================
   //

   // load macro for plugin setup
   gSystem->AddIncludePath("-I$ALICE_ROOT/include");
   gROOT->LoadMacro("SetupPlugin.C++g");
   
   // set run modalities
   PluginSetup::runMode      = runMode;
   PluginSetup::targetAlien  = kFALSE;

   // define inputs
   PluginSetup::dataSet      = dataSet;
   PluginSetup::proofTest    = testFile;
   PluginSetup::proofCluster = clusterName;

   // define additional libraries
   PluginSetup::addLibs = addLibs;
   PluginSetup::addPar  = addPars;
   
   // additional modalities
   PluginSetup::aliVersion = "v4-21-26-AN";
   
   // create the plugin (need to know if we want tender)
   if (!PluginSetup::CreatePlugin()) return;
   
   PluginSetup::plugin->SetRunMode(runMode);
   PluginSetup::plugin->SetOverwriteMode();
   PluginSetup::plugin->SetProofParameter("PROOF_UseMergers", "-1");
   if (nmix > 0) {::Info("run", "Forcing local"); PluginSetup::plugin->SetProofParameter("PROOF_ForceLocal", "1");}
   PluginSetup::plugin->SetRootVersionForProof("current");
   PluginSetup::plugin->SetAliRootMode("default");
   PluginSetup::plugin->SetClearPackages(kFALSE);
   
   //
   // === ANALYSIS EXECUTION =======================================================================
   //

   // add plugin to analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   mgr->SetGridHandler(PluginSetup::plugin);
   
   // initialize and start analysis
   if (!mgr->InitAnalysis()) {
      ::Error("runPlugin.C", "Failed to init analysis");
      return;
   }
   mgr->PrintStatus();
   mgr->StartAnalysis("proof");
}
