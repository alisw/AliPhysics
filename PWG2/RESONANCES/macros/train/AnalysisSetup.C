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
   Int_t       nmix,
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
   
   Bool_t isMix      = (nmix > 0);
   Bool_t isESD      = opt.Contains("ESD");
   Bool_t isAOD      = opt.Contains("AOD");
   Bool_t useTender  = opt.Contains("TENDER");
   Bool_t usePhysSel = opt.Contains("PHYS");
   Bool_t noV0       = opt.Contains("NOV0");
   Bool_t useCent    = opt.Contains("CENT");
   
   //
   // === LOAD LIBRARIES ===========================================================================
   //

   gSystem->Load("libVMC.so");
   gSystem->Load("libTree.so");
   gSystem->Load("libPhysics.so");
   gSystem->Load("libMatrix.so");
   gSystem->Load("libMinuit.so");
   gSystem->Load("libXMLParser.so");
   gSystem->Load("libGui.so");
   gSystem->Load("libSTEERBase.so");
   gSystem->Load("libESD.so");
   gSystem->Load("libAOD.so");
   gSystem->Load("libANALYSIS.so");
   gSystem->Load("libANALYSISalice.so");
   gSystem->Load("libEventMixing.so");
   gSystem->Load("libCORRFW.so");
   
   if (useTender) {
      ::Info("AnalysisSetup", "Loading tender libraries");
      gSystem->Load("libTENDER.so");
      gSystem->Load("libTENDERSupplies.so");
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

   ::Info("AnalysisSetup", "Using Multi Handler");
   // create multi input event handler
   AliMultiInputEventHandler *multiHandler = new AliMultiInputEventHandler();

   // since there is an exit point above if the job
   // isn't either ESD or AOD, here we don't recheck that
   if (isESD) {
      ::Info("AnalysisSetup", "Creating ESD handler");
      AliESDInputHandler *esdHandler = new AliESDInputHandler();
      multiHandler->AddInputEventHandler(esdHandler);

      // if possible, create also MC handler
      if (isMC) {
         ::Info("AnalysisSetup", "Creating MC handler");
         AliMCEventHandler *mcHandler  = new AliMCEventHandler();
         multiHandler->AddInputEventHandler(mcHandler);
      }
   } else {
      ::Info("AnalysisSetup", "Creating AOD handler");
      AliAODInputHandler *aodHandler = new AliAODInputHandler();
      multiHandler->AddInputEventHandler(aodHandler);
   }

   // add RSN input handler
   ::Info("AnalysisSetup", "Adding RSN input handler");
   gROOT->LoadMacro(Form("%s/AddRsnInputHandler.C", configPath));
   AddRsnInputHandler(isMC, multiHandler);
   
   // set the input event handler for manager
   mgr->SetInputEventHandler(multiHandler);
   
   // set event mixing properties, if required
   if (isMix) {
      ::Info("AnalysisSetup", "Setting mixing to n = %d", nmix);
      // define mixing handler
      AliMixInputEventHandler *mixHandler = new AliMixInputEventHandler(1, nmix);
      mixHandler->SetInputHandlerForMixing(multiHandler);

      // define binnings
      AliMixEventPool   *evPool  = new AliMixEventPool();
      AliMixEventCutObj *multip  = new AliMixEventCutObj(AliMixEventCutObj::kMultiplicity, 1., 10000., 1000.);
      AliMixEventCutObj *zvertex = new AliMixEventCutObj(AliMixEventCutObj::kZVertex, -10., 10., 4.);
      AliMixEventCutObj *cent    = new AliMixEventCutObj(AliMixEventCutObj::kCentrality, 0.0, 100.0, 10.0, "V0M");

      // add cuts to handler and handler to manager
      evPool->AddCut(zvertex);
      if (useCent) evPool->AddCut(cent); else evPool->AddCut(multip);
      mixHandler->SetEventPool(evPool);
      multiHandler->AddInputEventHandler(mixHandler);
   }

   //
   // === CONFIGURE AND INSERT PHYSICS SELECTION & TENDER SUPPLY ===================================
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
         gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPhysicsSelection.C");
         AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);
         if (noV0) {
            ::Info("AnalysisSetup", "options '%s' require to skip V0 info", options);
            physSelTask->GetPhysicsSelection()->SetSkipV0(kTRUE);
         }
         if (multiHandler) {
            AliESDInputHandler *esdHandler = dynamic_cast<AliESDInputHandler*>(multiHandler->GetFirstInputEventHandler());
            esdHandler->SetEventSelection(multiHandler->GetEventSelection());
         }
      }
   }
   
   // add PID response task
   //gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPIDResponse.C");
   //AddTaskPIDResponse(isMC);
   
   ::Info("AnalysisSetup", "Setup successful");
   return kTRUE;
}
