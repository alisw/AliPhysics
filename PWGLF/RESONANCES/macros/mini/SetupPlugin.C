#include <Riostream.h>
#include <TError.h>
#include <TObjArray.h>
#include <TObjString.h>
#include "AliAnalysisAlien.h"

namespace PluginSetup
{
   //
   // === DEFINITIONS ==============================================================================
   //
   
   // output object
   AliAnalysisAlien *plugin   = 0x0;                 // --> output of this code
   
   // basic setup
   Bool_t  targetAlien        = kTRUE;               // [*all*] kTRUE --> alien, kFALSE --> proof

   // executables
   TString aliVersion         = "v4-21-19-AN";       // [*all*] tag used for AliRoot
   TString rootVersion        = "v5-28-00a";         // [*all*] tag used for ROOT
   TString apiVersion         = "V1.1x";             // [AliEn] tag used for AliEn API
                                                       
   // file naming 
   TString jobName            = "";                  // [AliEn] basic name for all files                                   
   TString fileJDL            = "";                  // [AliEn] JDL filename
   TString fileC              = "";                  // [AliEn] run macro name
   TString fileSH             = "";                  // [AliEn] executable shell script name
                              
   // input definitions 
   Bool_t  alienInputRuns     = kTRUE;               // [AliEn] choose if inputs are runs (kTRUE) or lists of XML (kFALSE)      
   TString runList            = "";                  // [AliEn] list of run numbers from which XML collections are created
   TString runPath            = "";                  // [AliEn] common path where the above runs are stored
   TString runPrefix          = "";                  // [AliEn] a common string prefix to all run numbers (usually, "" for MC, "000" for data)
   TString runPattern         = "";                  // [AliEn] pattern to be used for the 'find -x' command
   TString xmlList            = "";                  // [AliEn] alternative to runs, pass a list of already-made XML files (NOTE: cannot do "test")
   TString dataSet            = "";                  // [PROOF] dataset name
   TString proofTest          = "";                  // [PROOF] test TXT collection of local files for running in "test" mode for PROOF
   TString proofCluster       = "alice-caf.cern.ch"; // [PROOF] name of used PROOF cluster
   Bool_t  proofReset         = kFALSE;              // [PROOF] tell the cluster if compiled libraries must be cleared
   TString alirootMode        = "default";           // [PROOF] modality for AliRoot loading
   Bool_t  proofClear         = kFALSE;              // [PROOF] tell the cluster if compiled libraries must be cleared
                              
   // output paths            
   TString workDir            = "";                  // [AliEn] working dir w.r. to user $ALIEN_HOME
   TString outDir             = "out";               // [AliEn] output dir w.r. to 'workDir'

   // additional libraries
   TString addTaskRuntime     = "";                  // [*all*] list of tasks compiled at runtime (just the name)
   TString addLibs            = "";                  // [*all*] list of all libraries to be added (order matters)
   TString addIncludes        = "";                  // [*all*] list of additional include paths to use
   TString addPar             = "";                  // [*all*] list of additional PARs to use
                                                   
   // AliEn JDL parameters                         
   Int_t   split              = 100;                 // [AliEn] split parameter in JDL
   Int_t   maxMergeFiles      = 50;                  // [AliEn] how many files to be merged per chunk
   Int_t   nRunsPerMaster     = 1;                   // [AliEn] how many runs are used together in a single master job (for ESD: 1, for AOD: even all)
   Int_t   maxMergeStages     = 2;                   // [AliEn] how many maximum merge stages to have
   Int_t   maxInitFailed      = 0;                   // [AliEn] how many failed initializations are tolerated (?)
   Int_t   resubmitThr        = 0;                   // [AliEn] threshold that triggers automatic job resubmit (?)
   Int_t   TTL                = 30000;               // [AliEn] time-to-live (in sec)
   TString inputFormat        = "xml-single";        // [AliEn] input format in JDL
   TString splitMode          = "se";                // [AliEn] splitting modality
   Int_t   price              = 1;                   // [AliEn] job price
   TString jobTag             = "";                  // [AliEn] tag assigned to job
   Int_t   nTestFiles         = 1;                   // [AliEn] number of input files used locally in test mode
   
   //
   // === FUNCTIONS ================================================================================
   //
   
   Bool_t AssignNames();
   Bool_t CreatePlugin();
   Bool_t SetupForAlien();
   Bool_t SetupForProof();
   
   //_______________________________________________________________________________________________
   //
   // Define all 'automatic' names from the job name
   //
   Bool_t AssignNames()
   {
      if (jobName.Length() < 1) {
         ::Error("SetupPlugin::AssignNames()", "Job name not defined");
         return kFALSE;
      }
         
      
      fileC = jobName;
      fileSH = jobName;
      fileJDL = jobName;
      
      fileC.Append(".C");
      fileSH.Append(".sh");
      fileJDL.Append(".jdl");
      
      ::Info("SetupPlugin::AssignNames()", "macro  name: \"%s\"", fileC.Data());
      ::Info("SetupPlugin::AssignNames()", "JDL    name: \"%s\"", fileJDL.Data());
      ::Info("SetupPlugin::AssignNames()", "script name: \"%s\"", fileSH.Data());
      
      return kTRUE;
   }
   
   //_______________________________________________________________________________________________
   //
   // Initialize plugin with all its common values
   //
   Bool_t CreatePlugin()
   {
      // create object
      plugin = new AliAnalysisAlien;
      
      // framework version
      plugin->SetROOTVersion(rootVersion.Data());
      plugin->SetAliROOTVersion(aliVersion.Data());
      
      // additional libraries/includes
      if (addLibs.Length() > 0) plugin->SetAdditionalLibs(addLibs.Data());
      if (addIncludes.Length() > 0) plugin->AddIncludePath(addIncludes.Data());
      
      // runtime tasks
      if (addTaskRuntime.Length() > 0) {
         TObjArray *list = addTaskRuntime.Tokenize(" ");
         TObjArrayIter next(list);
         TObjString *os = 0x0;
         TString sources("");
         while ( (os = (TObjString*)next()) ) {
            const char *taskName = os->GetString().Data();
            addLibs.Append(Form("%s.h %s.cxx", taskName, taskName));
            sources.Append(Form("%s.cxx", taskName));
         }
         plugin->SetAnalysisSource(sources.Data());
      }
      
      // specific setups
      if (targetAlien) {
         ::Info("SetupPlugin::CreatePlugin()", "Setting up for ALIEN");
         return SetupForAlien();
      } else {
         ::Info("SetupPlugin::CreatePlugin()", "Setting up for PROOF");
         return SetupForProof();
      }
   }
   
   
   //_______________________________________________________________________________________________
   //
   // Initialize plugin with all needed for AliEn
   //
   Bool_t SetupForAlien()
   {
      if (!plugin) {
         ::Error("SetupPlugin::SetupForAlien()", "Initialize plugin first");
         return kFALSE;
      }
      
      // create names
      if (!AssignNames()) {
         ::Error("SetupPlugin::SetupForAlien()", "Failed name initializations");
         return kFALSE;
      }
      
      // API version
      plugin->SetAPIVersion(apiVersion.Data());
      
      // merging detauls
      plugin->SetMergeViaJDL();
      plugin->SetMaxMergeFiles(maxMergeFiles);
      plugin->SetMaxMergeStages(maxMergeStages);
      
      // output paths
      plugin->SetGridWorkingDir(workDir.Data());
      plugin->SetGridOutputDir(outDir.Data());
      plugin->SetDefaultOutputs(kTRUE);
      
      // excutable
      plugin->SetExecutableCommand("aliroot -q -b");
      plugin->SetExecutableArgs(">& std.log");
      plugin->SetExecutable(fileSH.Data());
      
      // automatically created files
      plugin->SetAnalysisMacro(fileC.Data());
      plugin->SetJDLName(fileJDL.Data());
      
      // JDL parameters
      plugin->SetSplitMaxInputFileNumber(split);
      plugin->SetMaxInitFailed(maxInitFailed);
      plugin->SetMasterResubmitThreshold(resubmitThr);
      plugin->SetTTL(TTL);
      plugin->SetPrice(price);
      plugin->SetInputFormat(inputFormat.Data());
      plugin->SetSplitMode(splitMode.Data());
      if (jobTag.Length() > 0) plugin->SetJobTag(jobTag.Data());
   
      // input definition
      if (alienInputRuns) {
         plugin->SetOutputToRunNo(kTRUE);
         plugin->SetNtestFiles(nTestFiles);
         plugin->SetNrunsPerMaster(nRunsPerMaster);
         plugin->SetRunPrefix(runPrefix.Data());
         plugin->SetGridDataDir(runPath.Data());
         plugin->SetDataPattern(runPattern.Data());
         plugin->AddRunList(runList.Data());
      } else {
         TObjArray *list = xmlList.Tokenize(" ");
         TObjArrayIter next(list);
         TObjString *os = 0x0;
         while ( (os = (TObjString*)next()) ) {
            plugin->AddDataFile(os->GetString().Data());
         }
         plugin->SetOutputToRunNo(kFALSE);
      }
      
      return kTRUE;
   }
   
   //_______________________________________________________________________________________________
   //
   // Initialize plugin with all needed for PROOF
   //
   Bool_t SetupForProof()
   {
      if (!plugin) {
         ::Error("SetupPlugin::SetupForAlien()", "Initialize plugin first");
         return kFALSE;
      }
      
      plugin->SetProofCluster(proofCluster.Data());
      plugin->SetProofDataSet(dataSet.Data());
      plugin->SetProofReset(proofReset);
      plugin->SetProofConnectGrid(kTRUE);
      plugin->SetAliRootMode(alirootMode.Data());
      plugin->SetClearPackages(proofClear);
      plugin->SetFileForTestMode(proofTest.Data());
      
      return kTRUE;
   }
};
