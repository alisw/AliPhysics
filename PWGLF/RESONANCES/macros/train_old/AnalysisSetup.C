//
// This macro sets all the aspects of configuration of an Analysis Train run
// which are always the same for all kinds of analysis (local, PROOF, AliEn)
//
// Inputs:
//
//   - taskList       = a string containin all the 'add-task' macros to be used
//   - options        = a set of keywords which drive some configurations
//   - outputFileName = name of file produced by train
//   - configPath     = a path where all required config macros are stored
//
// Notes:
//
//   - in case the source is an ESD, and if inputs are a MC production
//     the MC input handler is created by default
//
Bool_t AnalysisSetup
(
   Bool_t      isMC,
   const char *options,
   const char *outputFileName,
   const char *configPath
)
{
   //
   // === EXAMINE OPTIONS ==========================================================================
   //

   // this is done using the utility 'RsnOptions.C'
   // which provides a unique way to interpret them
   
   TString opt(options);
   opt.ToUpper();
   opt.ReplaceAll(" ", "");
   
   Bool_t isMix      = opt.Contains("MIX");
   Bool_t isESD      = opt.Contains("ESD");
   Bool_t isAOD      = opt.Contains("AOD");
   Bool_t central    = opt.Contains("CEN");
   Bool_t peripheral = opt.Contains("PER");
   Bool_t useTender  = opt.Contains("Tender");
   Bool_t usePhysSel = opt.Contains("PHYS");
   Bool_t noV0       = opt.Contains("NOV0");
   
   //
   // === LOAD LIBRARIES ===========================================================================
   //

   gSystem->Load("libVMC");
   gSystem->Load("libTree");
   gSystem->Load("libPhysics");
   gSystem->Load("libMatrix");
   gSystem->Load("libMinuit");
   gSystem->Load("libXMLParser");
   gSystem->Load("libGui");
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libEventMixing");
   gSystem->Load("libCORRFW");
   
   if (useTender) {
      ::Info("AnalysisSetup", "Loading tender libraries");
      gSystem->Load("libTender");
      gSystem->Load("libTenderSupplies");
   }
   
   if (!AliAnalysisAlien::SetupPar("PWG2resonances.par")) return kFALSE;

   //
   // === ANALYSIS MANAGER CONFIGURATION ===========================================================
   //

   // create analysis manager
   AliAnalysisManager *mgr = new AliAnalysisManager("RsnAnalysisManager");
   mgr->SetCommonFileName(outputFileName);
   ::Info("AnalysisSetup", "Common file name: %s", outputFileName);

   //
   // === INPUT / OUTPUT HANDLER CONFIGURATION =====================================================
   //

   // create input handler
   // since there is an exit point above if the job
   // isn't either ESD or AOD, here we don't recheck that
   if (isESD) {
      ::Info("AnalysisSetup", "Configuring for ESD");
      AliESDInputHandler *esdHandler = new AliESDInputHandler();
      mgr->SetInputEventHandler(esdHandler);
      // if possible, create also MC handler
      if (isMC) {
         ::Info("AnalysisSetup", "Creating MC handler");
         AliMCEventHandler *mcHandler  = new AliMCEventHandler();
         mgr->SetMCtruthEventHandler(mcHandler);
      }
   } else if (isAOD) {
      ::Info("AnalysisSetup", "Configuring for AOD");
      AliAODInputHandler *aodHandler = new AliAODInputHandler();
      mgr->SetInputEventHandler(aodHandler);
   } else {
      ::Error("AnalysisSetup", "Require ESD or AOD");
      return kFALSE;
   }

   //
   // === CONFIGURE AND INSERT PHYSICS SELECTION & Tender SUPPLY ===================================
   //

   // add event selection for data if running ESD
   if (isESD) {
      // add tender supply for TOF
      if (useTender) {
         ::Info("AnalysisSetup", "options '%s' require to add tender", options);
         gROOT->LoadMacro(Form("%s/AddTenderSupplies.C", configPath));
         AddTenderSupplies(100.0, kTRUE, kFALSE);
      }

      // add event selection for data
      // and swtich off VZERO if tender is not used
      if (usePhysSel) {
         ::Info("AnalysisSetup", "options '%s' require to add physics selection", options);
         gROOT->LoadMacro("$(ALICE_ROOT)/OADB/macros/AddTaskPhysicsSelection.C");
         AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);
         if (noV0) {
            ::Info("AnalysisSetup", "options '%s' require to skip V0 info", options);
            physSelTask->GetPhysicsSelection()->SetSkipV0(kTRUE);
         }
      }
   }
   
   ::Info("AnalysisSetup", "Setup successful");
   return kTRUE;
}
