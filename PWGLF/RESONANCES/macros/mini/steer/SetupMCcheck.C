//fbellini
enum ERsnSpecie_t{kRho=0,
		  kKstar,
		  kAKstar,
		  kKstarP,
		  kKstarM,
		  kPhi,
		  kDelta0,
		  kDeltaPP,
		  kADeltaMM,
		  kLstar,
		  kALstar,
		  kF0,
		  kNrsn};

TString SetupMCcheck(Int_t       nmix,
		      const char *options,
		      const char *outputFileName,
		      const char *macroPath,
		      Bool_t      enableMon=kTRUE,
		      Bool_t      runMonOnly=kFALSE,
		      Bool_t      isMcTrueOnly=kTRUE,
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
  gSystem->Load("libEventMixing.so");
  gSystem->Load("libPWGLFresonances.so");
   
  // tender-related libraries
  if (isESD && useTender) {
    ::Info("AnalysisSetup", "Loading tender libraries");
    gSystem->AddIncludePath("-I$ALICE_PHYSICS/TENDER/Tender -I$ALICE_PHYSICS/TENDER/Tender -g");
    gSystem->Load("libTender.so");
    gSystem->Load("libTenderSupplies.so");
  } else if (!isESD) {
    useTender = kFALSE;
  }

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
  // ::Info("AnalysisSetup", "Add task for PID response");
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
  TString partname [12] = {"Rho","Kstar","aKstar","Kstarp","Kstarm", "Phi", "Delta0","Deltapp","Deltamm", "Lstar", "aLstar", "f0"};
  Int_t    pdgCode [12] = {113, 313, -313, 323, -323, 333, 2114, 2224, -2224, 3124, -3124, 9010221};
  Float_t mass [12] = {0.770, 0.89495, 0.89445, 0.89166, 0.89166, 1.019461, 1.232, 1.2349, 1.2349, 1.51954, 1.51954, 0.990};
  Float_t masslow [12] = {0.5, 0.7, 0.7, 0.7, 0.7, 0.9, 0.8, 0.8, 0.8, 1.3, 1.3, 0.8};
  Float_t massup  [12] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.1, 1.4, 1.4, 1.4, 1.7, 1.7, 1.2};
  Int_t   nbins   [12] = {250, 300, 300, 300, 300, 200, 300, 300, 300, 200, 200, 200};
  Char_t  charge1 [12] = {'+', '+','-','0','0','+','+','+','-','+','-','+'};
  Char_t  charge2 [12] = {'-', '-','+','+','-','-','-','+','-','-','+','-'};
  RSNPID d1[12] = { AliRsnDaughter::kPion,  AliRsnDaughter::kKaon,  AliRsnDaughter::kKaon, 
		    AliRsnDaughter::kKaon0, AliRsnDaughter::kKaon0, AliRsnDaughter::kKaon,
		    AliRsnDaughter::kProton,AliRsnDaughter::kProton,AliRsnDaughter::kProton,
		    AliRsnDaughter::kProton, AliRsnDaughter::kProton, AliRsnDaughter::kPion};
  RSNPID d2[12] = { AliRsnDaughter::kPion,  AliRsnDaughter::kPion,  AliRsnDaughter::kPion,
		    AliRsnDaughter::kPion,  AliRsnDaughter::kPion, AliRsnDaughter::kKaon,
		    AliRsnDaughter::kPion, AliRsnDaughter::kPion, AliRsnDaughter::kPion,
		    AliRsnDaughter::kKaon, AliRsnDaughter::kKaon, AliRsnDaughter::kPion};

  gROOT->LoadMacro("${ALICE_PHYSICS}/PWGLF/RESONANCES/macros/mini/steer/AddTaskMC.C");
  for (Int_t i=0;i<ERsnSpecie_t::kNrsn;i++){
    AddTaskMC(isPP, Form("rsn%i",i), 0, 0, 5, partname[i].Data(), pdgCode[i], mass[i], masslow[i], massup[i], nbins[i], charge1[i], charge2[i], d1[i],d2[i], (i==0), "NoSIGN");
  }  
  ::Info("AnalysisSetup", "Setup successful");
  return out;
}

  
  
