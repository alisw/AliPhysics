// ======================================
// ===== ALIEN PLUGIN CONFIGURATION =====
// ======================================
//
// This macro configures an AliEn plugin in order to launch a job
// which runs a task from the resonance package on a list of runs
// which are processed separately.
//
// All the possible configuration parameters are arguments
// of the macro function, even if most of them have default
// values which the user will rarely change.
//
// The macro tries to synchronize some output names, using
// a unique name ('analysisName') to define all files that
// describe the output, the analysis macros/executables/JDL.
//
// Since the run mode can be more variable than the config
// it is not set here, but it is required in the run macro
// which uses the plugin.
//
// Considered that the arguments are many, they are explained
// inside the list of arguments in the macro definition.
// In ALL cases where a list of strings must be provided, its
// elements must be separated by a blank character.
//
AliAnalysisAlien* PluginByRunProof
(
   // all parameters which could often be customized
   // are placed at the beginning of the macro, while
   // all others can stay there with their default values:
   // -- analysisName   --> common name used for naming all analysis related files
   // -- dataset        --> dataset to be processed
   // -- testFile       --> used for test mode only
   // -- clusterName    --> PROOF cluster to be used
   const char *dataSet,
   const char *testFile,
   const char *clusterName,
   
   // -- proofReset     --> require or not the PROOF reset
   // -- alirootMode    --> the way AliROOT library are loaded
   // -- clearPack      --> to clear or not the PAR libraries
   Bool_t      proofReset  = kFALSE,
   const char *alirootMode = "default",
   Bool_t      clearPack   = kFALSE,

   // standard package versions
   const char *rootVersion    = "v5-28-00a",
   const char *aliVersion     = "v4-21-17a-AN"
)
{
   // create plugin object
   // with specifications which apply to a run-by-run execution
   // this creates by default also the job structure for merging
   AliAnalysisAlien *plugin = new AliAnalysisAlien;

   // package versions
   plugin->SetROOTVersion(rootVersion);
   plugin->SetAliROOTVersion(aliVersion);

   // additional libraries/includes
   //plugin->SetAdditionalLibs("libTENDER.so TENDERSupplies.par libEventMixing.so libPWG2resonances.so");
   plugin->SetAdditionalLibs("libEventMixing.so PWG2resonances.par");
   
   // PROOF-specific settings
   plugin->SetProofCluster(clusterName);
   plugin->SetProofDataSet(dataSet);
   plugin->SetProofReset(proofReset);
   plugin->SetProofConnectGrid(kTRUE);
   plugin->SetAliRootMode(alirootMode);
   plugin->SetClearPackages(clearPack);
   plugin->SetFileForTestMode(testFile);

   // the end!
   return plugin;
}
