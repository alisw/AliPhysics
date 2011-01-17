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
Bool_t RsnConfigMC(const char *taskName, const char *options)
{
  // info
  Info("RsnConfig2010PhiFcnMC", "Starting configuration");

  // retrieve analysis manager & task
  AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
  AliRsnAnalysisSE   *task = (AliRsnAnalysisSE*)mgr->GetTask(taskName);

  // for safety, return if no task is passed
  if (!task)
  {
    Error("ConfigTaskRsn", "Task not found");
    return kFALSE;
  }

  // interpret the useful information from second argument
  TString opt(options);
  Bool_t isMC    = opt.Contains("MC");
  Bool_t isSim   = opt.Contains("sim");
  Bool_t isData  = opt.Contains("data");
  Bool_t isPass1 = opt.Contains("pass1");
  Bool_t isPass2 = opt.Contains("pass2");
  if (!isMC)
  {
    Info("RsnConfig2010PhiFcnMC", "Config skipped for not pure MonteCarlo samples");
    return kTRUE;
  }

  //
  // -- Setup pairs ---------------------------------------------------------------------------------
  //

  // decay channels
  AliRsnPairDef *pairDefPM = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455);
  AliRsnPairDef *pairDefPP = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '+', 333, 1.019455);
  AliRsnPairDef *pairDefMM = new AliRsnPairDef(AliPID::kKaon, '-', AliPID::kKaon, '-', 333, 1.019455);

  // computation objects
  AliRsnPairFunctions *pairPM = new AliRsnPairFunctions("PairPM_mc", pairDefPM);
  AliRsnPairFunctions *truePM = new AliRsnPairFunctions("TruePM_mc", pairDefPM);
  AliRsnPairFunctions *pairPP = new AliRsnPairFunctions("PairPP_mc", pairDefPP);
  AliRsnPairFunctions *pairMM = new AliRsnPairFunctions("PairMM_mc", pairDefMM);

  //
  // -- Setup cuts ----------------------------------------------------------------------------------
  //

  // track cut -----------------------------
  // --> perfect PID for check of PID issues
  AliRsnCutPID *cutPID = new AliRsnCutPID("cutPID", AliPID::kKaon, 0.0, kTRUE);

  // cut sets ---------------------------------
  // --> only common cuts for tracks are needed
  // --> standard 2010 cuts are applied only when working on ESD
  AliRsnCutSet *cutSetDaughterCommon = new AliRsnCutSet("commonDaughterCuts", AliRsnCut::kDaughter);
  cutSetDaughterCommon->AddCut(cutPID);
  cutSetDaughterCommon->SetCutScheme("cutPID");
  cout << "Cut scheme: " << cutSetDaughterCommon->GetCutScheme() << endl;

  // configure cut managers -------------------
  pairPM->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);
  truePM->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);
  pairPP->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);
  pairMM->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);

  // set additional option for true pairs when needed
  truePM->SetOnlyTrue  (kTRUE);
  truePM->SetCheckDecay(kTRUE);

  //
  // -- Setup functions -----------------------------------------------------------------------------
  //

  // function axes
  Double_t y[] = {-0.8, -0.7, -0.6, -0.5, 0.5, 0.6, 0.7, 0.8};
  Int_t    ny  = sizeof(y) / sizeof(y[0]);
  AliRsnValue *axisIM = new AliRsnValue("IM", AliRsnValue::kPairInvMass, 2000, 0.9,  2.9);
  AliRsnValue *axisPt = new AliRsnValue("PT", AliRsnValue::kPairPt,       100, 0.0, 10.0);
  AliRsnValue *axisY  = new AliRsnValue("Y" , AliRsnValue::kPairY,         ny, y);

  // create function and add axes
  AliRsnFunction *fcnImPtY = new AliRsnFunction;
  fcnImPtY->AddAxis(axisIM);
  fcnImPtY->AddAxis(axisPt);
  fcnImPtY->AddAxis(axisY);

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
  task->GetAnalysisManager()->Add(truePM);

  return kTRUE;
}
