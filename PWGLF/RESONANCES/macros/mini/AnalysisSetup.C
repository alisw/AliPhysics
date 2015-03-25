//fbellini

TString AnalysisSetup(Int_t       nmix,
		      const char *options,
		      const char *outputFileName,
		      const char *macroPath,
		      Bool_t      isAOD49=kFALSE,
		      Bool_t      enableMon=kTRUE,
		      Bool_t      runMonOnly=kFALSE,
		      Bool_t      isMcTrueOnly=kFALSE,
		      TString     monitorOpt = "NoSIGN")
{
   // prepare output
   TString out("");
   
   // === EXAMINE OPTIONS ==========================================================================
   TString opt(options);
   opt.ToUpper();
   
   Bool_t isMC      = opt.Contains("MC") || (!opt.Contains("DATA"));
   Bool_t isPP      = opt.Contains("PP") || (!opt.Contains("PBPB"));
   Bool_t isESD     = opt.Contains("ESD");
   Bool_t useTender = opt.Contains("TENDER");
   Bool_t noV0      = opt.Contains("NOV0");
   
   //
   // === LOAD LIBRARIES ===========================================================================
   // uncomment only if not already done in the steering macro, eg. RunGrid.C
   // load analysis libraries
   gSystem->AddIncludePath("-I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/ANALYSISalice");
   gSystem->Load("libTree.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libPhysics.so");
   gSystem->Load("libSTEERBase.so");
   gSystem->Load("libESD.so");
   gSystem->Load("libAOD.so");
   gSystem->Load("libANALYSIS.so");
   gSystem->Load("libOADB");
   gSystem->Load("libANALYSISalice.so");
   gSystem->Load("libEventMixing.so");
   gSystem->Load("libPWGLFresonances.so");
   
   // tender-related libraries
   if (isESD && useTender) {
      ::Info("AnalysisSetup", "Loading tender libraries");
      gSystem->Load("libTENDER.so");
      gSystem->Load("libTENDERSupplies.so");
   } else if (!isESD) {
      useTender = kFALSE;
   }
   
   // load development RSN library
   // uncomment only for par files usage (not recommended)
   // if (!AliAnalysisAlien::SetupPar("PWGLFresonances.par")) return "";

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
      esdHandler->SetNeedField();
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
   // === PHYSICS SELECTION =============================================================
   // Check that the correct trigger class is selected!!! USER DEFINED!!!

   if (isESD) {
     ::Info("AnalysisSetup", "Add physics selection");
     gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
     AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);
     physSelTask->SelectCollisionCandidates(AliVEvent::kMB); 
     if (noV0) {
       ::Info("AnalysisSetup", "Skip of V0 info is required");
       physSelTask->GetPhysicsSelection()->SetSkipV0(kTRUE);
     }
   }
   
   //
   // === CENTRALITY/PLANE ==============================================================
   //
   if (isESD & !isPP) {
     ::Info("AnalysisSetup", "Add centrality and event plane computation tasks");
     gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
     AliCentralitySelectionTask* taskCentrality = (AliCentralitySelectionTask*)AddTaskCentrality();
     if (isMC) {
       ::Info("AnalysisSetup", "Setting centrality computation for MC");
       taskCentrality->SetMCInput();
     } 
   }
   
   //
   // === PID RESPONSE =============================================================================
   //   
   ::Info("AnalysisSetup", "Add task for PID response");
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   AddTaskPIDResponse(isMC, 
		      isMC, //autoMCesd
		      isMC, //tuneOnData
		      2, //recopass
                      kFALSE, //cachePID, default
		      "",//detResponse, default
                      kTRUE,//useTPCEtaCorrection,/*Please use default value! Otherwise splines can be off*/
                      kFALSE, //useTPCMultiplicityCorrection --> default was kTRUE, but not avail for LHC13b2_efix_p1 
		      -1); //recodatapass, default=-1
   
   //
   // === PID QA ===================================================================================
   //   
   gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPIDqa.C ");
   AddTaskPIDqa();
   
   //
   // === RSN TASKS ==============================================================================
   //    
   ::Info("AnalysisSetup", "Add task for KStar");
   gROOT->LoadMacro("${ALICE_PHYSICS}/PWGLF/RESONANCES/macros/mini/AddTaskRhoPP7TeV.C");
   AddTaskRhoPP7TeV(isMC, isPP);

     ::Info("AnalysisSetup", "Setup successful");
   return out;
}

  
  
