//
// This function configures the entire task for all resonances the user is interested in.
// This is done by creating all configuration objects which are defined in the package.
//
// Generally speaking, one has to define the following objects for each resonance:
//
//  1 - an AliRsnPairDef to define the resonance decay channel to be studied
//  2 - an AliRsnPair{Ntuple|Functions} where the output is stored
//  3 - one or more AliRsnCut objects to define track selections
//      which will have then to be organized into AliRsnCutSet objects
//  4 - an AliRsnCutManager to include all cuts to be applied (see point 3)
//  5 - definitions to build the TNtuple or histograms which are returned
//
// The return value is used to know if the configuration was successful
//
Bool_t RsnConfig
(
  const char *taskName, 
  const char *options,
  const char *config
)
{
  // load useful macros
  gROOT->LoadMacro("$(ALICE_ROOT)/PWG2/RESONANCES/macros/train/LHC2010-7TeV-phi/ConfigESDCutsITS.C");
  gROOT->LoadMacro("$(ALICE_ROOT)/PWG2/RESONANCES/macros/train/LHC2010-7TeV-phi/ConfigESDCutsTPC.C");
  
  // interpret the useful information from second argument
  TString opt(options);
  Bool_t isSim  = opt.Contains("sim");
  Bool_t isData = opt.Contains("data");
  Bool_t isESD  = opt.Contains("ESD");
  Bool_t isAOD  = opt.Contains("AOD");
  
  // interpret the specific info from third argument
  // which should be fixed in the various calls to this function
  TString opt(config);
  Bool_t realPID     = opt.Contains("realistic");
  Bool_t perfPID     = opt.Contains("perfect");
  Bool_t addITSSA    = opt.Contains("its");
  Bool_t dipAngleCut = opt.Contains("dip");
  Int_t  typePID     = 0;
  if (realPID) typePID = 1;
  else if (perfPID) typePID = 2;
      
  // info
  const Char_t *pid[3] = {"nopid", "realistic", "perfect"};
  Info("RsnConfig2010PhiFcn", "=== Specific configuration: %s ===", config);
  Info("RsnConfig2010PhiFcn", "--> PID           : %s", pid[typePID]);
  Info("RsnConfig2010PhiFcn", "--> ITS standalone: %s", (addITSSA ? "INCLUDED" : "EXCLUDED"));
  Info("RsnConfig2010PhiFcn", "--> dip-angle cut : %s", (dipAngleCut ? "INCLUDED" : "EXCLUDED"));
  
  // generate a common suffix depending on chosen options
  TString suffix(pid[typePID]);
  if (addITSSA) suffix += "_sa"; else suffix += "_nosa";
  if (dipAngleCut) suffix += "_dip";
  Info("RsnConfig2010PhiFcn", "--> suffix used   : %s", suffix.Data());

  // retrieve analysis manager & task
  AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
  AliRsnAnalysisSE   *task = (AliRsnAnalysisSE*)mgr->GetTask(taskName);

  // for safety, return if no task is passed
  if (!task)
  {
    Error("RsnConfig2010PhiFcn", "Task not found");
    return kFALSE;
  }

  //
  // -- Setup pairs ---------------------------------------------------------------------------------
  //

  // decay channels
  AliRsnPairDef *pairDefPM = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455);
  AliRsnPairDef *pairDefPP = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '+', 333, 1.019455);
  AliRsnPairDef *pairDefMM = new AliRsnPairDef(AliPID::kKaon, '-', AliPID::kKaon, '-', 333, 1.019455);

  // computation objects
  AliRsnPairFunctions *pairPM = new AliRsnPairFunctions(Form("PairPM_%s", suffix.Data()), pairDefPM);
  AliRsnPairFunctions *truePM = new AliRsnPairFunctions(Form("TruePM_%s", suffix.Data()), pairDefPM);
  AliRsnPairFunctions *pairPP = new AliRsnPairFunctions(Form("PairPP_%s", suffix.Data()), pairDefPP);
  AliRsnPairFunctions *pairMM = new AliRsnPairFunctions(Form("PairMM_%s", suffix.Data()), pairDefMM);

  //
  // -- Setup cuts ----------------------------------------------------------------------------------
  //

  // track cut -----------------------
  // --> global cuts for 2010 analysis
  // --> most options are set to right values by default
  AliRsnCutESD2010 *cuts2010_esd = new AliRsnCutESD2010(Form("cuts2010_esd_%s", suffix.Data()));
  AliRsnCutAOD2010 *cuts2010_aod = new AliRsnCutAOD2010(Form("cuts2010_aod_%s", suffix.Data()));
  // ----> set the flag for sim/data management (which sets some other options)
  cuts2010_esd->SetMC(isSim);
  cuts2010_aod->SetMC(isSim);
  // ----> include or not the ITS standalone (TPC is always in)
  cuts2010_esd->SetUseITSTPC(kTRUE);
  cuts2010_esd->SetUseITSSA (addITSSA);
  //cuts2010_aod->SetUseITSTPC(kTRUE);
  //cuts2010_aod->SetUseITSSA (addITSSA);
  // ----> require to check PID or not, depending on the label
  if (realPID)
  {
    // if doing realistic PID, it must be activated
    cuts2010_esd->SetCheckITS (kTRUE);
    cuts2010_esd->SetCheckTPC (kTRUE);
    cuts2010_esd->SetCheckTOF (kTRUE);
    cuts2010_aod->SetCheckITS (kTRUE);
    cuts2010_aod->SetCheckTPC (kTRUE);
    cuts2010_aod->SetCheckTOF (kTRUE);
  }
  else
  {
    // otherwise (both for no pid and perfect PID)
    // the PID cuts are deactivated
    cuts2010_esd->SetCheckITS (kFALSE);
    cuts2010_esd->SetCheckTPC (kFALSE);
    cuts2010_esd->SetCheckTOF (kFALSE);
    cuts2010_aod->SetCheckITS (kFALSE);
    cuts2010_aod->SetCheckTPC (kFALSE);
    cuts2010_aod->SetCheckTOF (kFALSE);
  }
  // ----> set all other defaults
  ConfigESDCutsTPC(cuts2010_esd->GetCutsTPC());
  ConfigESDCutsITS(cuts2010_esd->GetCutsITS());
  
  // track cut -----------------------------
  // --> perfect PID for check of PID issues
  AliRsnCutPID *cutPID = new AliRsnCutPID("cutPID", AliPID::kKaon, 0.0, kTRUE);
  
  // pair cut ----------------------
  // --> dip angle between daughters
  AliRsnCutValue *cutDip = new AliRsnCutValue("cutDip", AliRsnValue::kPairDipAngle, 0.03, 1.01);

  // cut set for tracks------------------------
  // --> only common cuts for tracks are needed
  // --> standard 2010 cuts are applied always
  TString cutSchemeTrack;
  if (isESD)
  {
    pairPM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010_esd);
    truePM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010_esd);
    pairPP->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010_esd);
    pairMM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010_esd);
    
    cutSchemeTrack += cuts2010_esd->GetName();
  }
  else if (isAOD)
  {
    pairPM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010_aod);
    truePM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010_aod);
    pairPP->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010_aod);
    pairMM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010_aod);
    
    cutSchemeTrack += cuts2010_aod->GetName();
  }
  else
  {
    Error("Required ESD or AOD");
    return kFALSE;
  }
  if (perfPID)
  {
    pairPM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cutPID);
    truePM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cutPID);
    pairPP->GetCutManager()->GetCommonDaughterCuts()->AddCut(cutPID);
    pairMM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cutPID);
    
    cutSchemeTrack += "&cutPID";
  }
  pairPM->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cutSchemeTrack.Data());
  truePM->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cutSchemeTrack.Data());
  pairPP->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cutSchemeTrack.Data());
  pairMM->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cutSchemeTrack.Data());
  
  // cut set for pairs---------------------------------------
  // --> add dip angle cut (but then include only if required)
  AliRsnCutSet *cutSetPair = new AliRsnCutSet("cutsPair", AliRsnCut::kMother);
  cutSetPair->AddCut(cutDip);
  cutSetPair->SetCutScheme(cutDip->GetName());
  if (dipAngleCut)
  {
    pairPM->GetCutManager()->GetMotherCuts()->AddCut(cutDip);
    truePM->GetCutManager()->GetMotherCuts()->AddCut(cutDip);
    pairPP->GetCutManager()->GetMotherCuts()->AddCut(cutDip);
    pairMM->GetCutManager()->GetMotherCuts()->AddCut(cutDip);
    
    pairPM->GetCutManager()->GetMotherCuts()->SetCutScheme(cutDip->GetName());
    truePM->GetCutManager()->GetMotherCuts()->SetCutScheme(cutDip->GetName());
    pairPP->GetCutManager()->GetMotherCuts()->SetCutScheme(cutDip->GetName());
    pairMM->GetCutManager()->GetMotherCuts()->SetCutScheme(cutDip->GetName());
  }
  
  // set additional option for true pairs when needed
  truePM->SetOnlyTrue  (kTRUE);
  truePM->SetCheckDecay(kTRUE);

  //
  // -- Setup functions -----------------------------------------------------------------------------
  //

  // axis definition
  // 0) invariant mass
  // 1) transverse momentum
  // 2) rapidity
  // 3) multiplicity
  Double_t mult[] = {0., 6., 10., 15., 23., 1E6};
  Int_t    nmult  = sizeof(mult) / sizeof(mult[0]);
  AliRsnValue *axisIM   = new AliRsnValue("IM"  , AliRsnValue::kPairInvMass     , 0.9,  1.4, 0.001);
  AliRsnValue *axisPt   = new AliRsnValue("PT"  , AliRsnValue::kPairPt          , 0.0, 10.0, 0.100);
  AliRsnValue *axisY    = new AliRsnValue("Y"   , AliRsnValue::kPairY           ,-1.2,  1.2, 0.100);
  AliRsnValue *axisMult = new AliRsnValue("Mult", AliRsnValue::kEventMultESDCuts, nmult, mult);
  
  // add the support cut to the value which computes the multiplicity
  AliESDtrackCuts *cuts = new AliESDtrackCuts;
  ConfigESDCutsTPC(cuts);
  axisMult->SetSupportObject(cuts);

  // create function and add axes
  AliRsnFunction *fcnImPtY = new AliRsnFunction;
  if ( !fcnImPtY->AddAxis(axisIM  ) ) return kFALSE;
  if ( !fcnImPtY->AddAxis(axisPt  ) ) return kFALSE;
  if ( !fcnImPtY->AddAxis(axisY   ) ) return kFALSE;
  if ( !fcnImPtY->AddAxis(axisMult) ) return kFALSE;

  // add functions to pairs
  pairPM->AddFunction(fcnImPtY);
  truePM->AddFunction(fcnImPtY);
  pairPP->AddFunction(fcnImPtY);
  pairMM->AddFunction(fcnImPtY);

  //
  // -- Conclusion ----------------------------------------------------------------------------------
  //

  // add all created AliRsnPair objects to the AliRsnAnalysisManager in the task
  task->GetAnalysisManager()->Add(pairPM);
  task->GetAnalysisManager()->Add(pairPP);
  task->GetAnalysisManager()->Add(pairMM);
  if (isSim) task->GetAnalysisManager()->Add(truePM);

  return kTRUE;
}
