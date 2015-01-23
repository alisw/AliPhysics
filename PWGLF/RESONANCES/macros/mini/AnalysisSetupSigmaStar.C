// Created by Enrico Fragiacomo - 13/02/14
// Based on AnalysisSetupRsnMini
//
// This macro sets all those aspects of configuration of an Analysis Train run
// which are always the same for all kinds of analysis (local, PROOF, AliEn)
//
// Inputs:
//
//   - options        = a set of keywords which drive some configurations
//   - outputFileName = name of file produced by train
//   - configPath     = a path where all required config macros are stored
//
// Returns:
//   
//   - if successful: the name of the expected input TTre (esdTree or aodTree)
//   - if failed    : NULL
//

enum ERsnCollType_t { kPP=0,
		      kPPb,
		      kPbPb};

TString Setup
(
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

   TString opt(options);
   opt.ToUpper();
   
   Bool_t isMC      = opt.Contains("MC") || (!opt.Contains("DATA"));
   Bool_t isESD     = opt.Contains("ESD");
   Bool_t useTender = opt.Contains("Tender");
   Bool_t noV0      = opt.Contains("NOV0");

   Int_t collSyst;
   if( opt.Contains("PP") && (!opt.Contains("PB")) ) collSyst=kPP;
   else if( opt.Contains("PPB") ) collSyst=kPPb;
   else collSyst=kPbPb;

   //
   // === LOAD LIBRARIES ===========================================================================
   //

   // load analysis libraries
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libEventMixing");
   gSystem->Load("libCORRFW");
   
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

   if (isESD && !(collSyst==kPP) ) {
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
   AddTaskPIDResponse(isMC);
      
   //
   // === MAIN TASK FOR THE SIGMASTAR ANALYSIS ==============================================================================
   //
   
   gROOT->LoadMacro(Form("%s/AddTaskSigmaStar.C", macroPath));
   if (!AddTaskSigmaStar(isMC, collSyst, 10.0 ,0, 0, 0, 5, 3.0, 3.0, 7.0, 0.01, 0.3, 0.97, 0.5, 70, 1.0, 10.0, 20.0, 68, "Sigma1385")) return ""; 

   ::Info("AnalysisSetup", "Setup successful");
   return out;
}
