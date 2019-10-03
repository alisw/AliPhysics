/* Macro for setup of Lambda* analysis in PPb collisions.
// Created by Sarita Sahoo, 30 Jan. 2014
//
//
// This macro sets all the aspects of configuration of an Analysis Train run
// which are always the same for all kinds of analysis (local, PROOF, AliEn)
//
// Inputs:
//
//   - nmix           = number of mixings to do (if > 0, initialize mixing stuff)
//   - options        = a set of keywords which drive some configurations
//   - outputFileName = name of file produced by train
//   - configPath     = a path where all required config macros are stored
//
// Notes:
//
//   - in case the source is an ESD, and if inputs are a MC production
//     the MC input handler is created by default
//
// Returns:
//   
//   - if successful: the name of the expected input TTree (esdTree or aodTree)
//   - if failed    : NULL
*/

TString Setup
(
   Int_t       nmix,
   const char *options,
   const char *outputFileName,
   const char *macroPath = "."
)
{
  
  
  // prepare output
   TString out("");
   
   //
   // === EXAMINE OPTIONS ==========================================================================
   //
   
   // this is done using the utility 'RsnOptions.C'
   // which provides a unique way to interpret them
   
   TString opt(options);
   opt.ToUpper();
   
   Bool_t isMC      = opt.Contains("MC") || (!opt.Contains("DATA"));
   Bool_t isPP      = opt.Contains("PP") || (!opt.Contains("PBPB"));
   Bool_t isESD     = opt.Contains("ESD");
   Bool_t useTender = opt.Contains("Tender");
   Bool_t noV0      = opt.Contains("NOV0");
   
   //
   // === LOAD LIBRARIES ===========================================================================
   //

   // load analysis libraries
   gSystem->Load("libCore");
   gSystem->Load("libGeom");
   gSystem->Load("libVMC");
   gSystem->Load("libMinuit");
   gSystem->Load("libPhysics");
   gSystem->Load("libTree");
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libEventMixing");
   gSystem->Load("libCORRFW");
   gSystem->Load("libPWGLFresonances");


   // tender-related libraries
   if (isESD && useTender) {
      ::Info("AnalysisSetup", "Loading tender libraries");
      gSystem->Load("libTender");
      gSystem->Load("libTenderSupplies");
   } else if (!isESD) {
      useTender = kFALSE;
   }
   
   // load development RSN library
   if (!AliAnalysisAlien::SetupPar("PWGLFresonances.par")) return "";

   //
   // === CREATE ANALYSIS MANAGER ==================================================================
   //

   AliAnalysisManager *mgr = new AliAnalysisManager("RsnAnalysisManager");
   mgr->SetCommonFileName(outputFileName);
   ::Info("AnalysisSetup", "Common file name: %s", outputFileName);

   //
   // === INPUT / OUTPUT HANDLER CONFIGURATION =====================================================
   //

   if (isESD) {
      out = "esdTree";
      ::Info("AnalysisSetup", "Creating ESD handler");
      AliESDInputHandler *esdHandler = new AliESDInputHandler();
      mgr->SetInputEventHandler(esdHandler);
      if (isMC) {
         ::Info("AnalysisSetup", "Creating MC handler");
         AliMCEventHandler *mcHandler  = new AliMCEventHandler();
         mgr->SetMCtruthEventHandler(mcHandler);
      }
   } else {
      out = "aodTree";
      ::Info("AnalysisSetup", "Creating AOD handler");
      AliAODInputHandler *aodHandler = new AliAODInputHandler();
      mgr->SetInputEventHandler(aodHandler);
   }
   
   //
   // === Tender TASK (ESD only -- optional) =======================================================
   //

   if (isESD && useTender) {
      ::Info("AnalysisSetup", "Adding tender (and then accepting V0 info)", options);
      gROOT->LoadMacro(Form("%s/AddTaskTender.C", macroPath));
      AddTaskTender();
      noV0 = kFALSE;
   }

   //
   // === PHYSICS SELECTION (ESD only) =============================================================
   //

   if (isESD) {
      ::Info("AnalysisSetup", "Add physics selection by default on ESD analysis");
      gROOT->LoadMacro("$(ALICE_PHYSICS)/OADB/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);
      if (noV0) {
         ::Info("AnalysisSetup", "Skip of V0 info is required");
         physSelTask->GetPhysicsSelection()->SetSkipV0(kTRUE);
      }
   }
   
   //
   // === CENTRALITY/PLANE (ESD only) ==============================================================
   //
   if (isESD && !isPP) {
     ::Info("AnalysisSetup", "Add centrality and event plane computation tasks");
      gROOT->LoadMacro("$(ALICE_PHYSICS)/OADB/macros/AddTaskCentrality.C");
      gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskEventplane.C");
      AliCentralitySelectionTask* taskCentrality = (AliCentralitySelectionTask*)AddTaskCentrality();
      if (isMC) {
	::Info("AnalysisSetup", "Setting centrality computation for MC");
	taskCentrality->SetMCInput();
      }
      AddTaskEventplane();
   }
   //
   // === PID RESPONSE =============================================================================
   //
   
   gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPIDResponse.C");
   AddTaskPIDResponse(isMC,kTRUE,kTRUE);

   //gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPIDqa.C ");
   //AddTaskPIDqa();
   
   //
   // === OTHER TASKS ==============================================================================
   //
   
   // add RSN task
   gROOT->LoadMacro(Form("%s/AddTaskLStarPPb.C", macroPath));
   if (!AddTaskLStarPPb(isMC, isPP, macroPath, nmix)) return "";
   
   ::Info("AnalysisSetup", "Setup successful");
   return out;
}
