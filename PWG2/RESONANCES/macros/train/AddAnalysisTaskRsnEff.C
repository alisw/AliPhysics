//
// This macro add an analysis task for computing efficiency.
// It will have as output an AliCFContainer with several steps:
//
//  0) all resonances in MC which decay in the pair specified
//  1) subset of (0) whose daughters are in acceptance
//  2) subset of (1) whose daughters satisfy quality track cuts (covariance, chi square && nTPCclusters)
//  3) subset of (2) whose daughters satisfy primary track cuts (nsigma to vertex, no kink daughters)
//  4) subset of (3) whose daughters satisty the BB TPC compatibility cut at 3 sigma
//
Bool_t AddAnalysisTaskRsnEff
(
  Bool_t      useBB    = kFALSE,
  Double_t    sigmaTPC = 0.065,
  const char *outFile  = "eff"
)
{
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // common suffixes
  TString suf[2];
  suf[0] = "nopid";
  suf[1] = "pid";

  // create task
  AliRsnAnalysisEffSE *task[2];
  task[0] = new AliRsnAnalysisEffSE("EffNoPID");
  task[1] = new AliRsnAnalysisEffSE("EffPID");

  // set prior probabilities for PID
  for (Int_t i = 0; i < 2; i++)
  {
    task[i]->SetPriorProbability(AliPID::kElectron, 0.02);
    task[i]->SetPriorProbability(AliPID::kMuon,     0.02);
    task[i]->SetPriorProbability(AliPID::kPion,     0.83);
    task[i]->SetPriorProbability(AliPID::kKaon,     0.07);
    task[i]->SetPriorProbability(AliPID::kProton,   0.06);
    task[i]->DumpPriors();
  }

  // pair definitions:
  // phi   --> K+ K-
  // kstar --> K+ pi- & K- pi+
  AliRsnPairDef *pairPhi    = new AliRsnPairDef('+', AliPID::kKaon, '-', AliPID::kKaon, 333);
  AliRsnPairDef *pairKStar1 = new AliRsnPairDef('-', AliPID::kKaon, '+', AliPID::kPion, 313);
  AliRsnPairDef *pairKStar2 = new AliRsnPairDef('+', AliPID::kKaon, '-', AliPID::kPion, 313);
  for (Int_t i = 0; i < 2; i++)
  {
    task[i]->AddPairDef(pairPhi);
    task[i]->AddPairDef(pairKStar1);
    task[i]->AddPairDef(pairKStar2);
  }

  // axis definition
  //  0) transverse momentum
  //  1) pseudo-rapidity
  //  2) multiplicity (estimated with SPD tracklets - uncorrected)
  AliRsnFunctionAxis *axisPt   = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairPt,       50,  0.0,  10.0);
  AliRsnFunctionAxis *axisEta  = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairEta,      20, -1.5,   1.5);
  AliRsnFunctionAxis *axisMult = new AliRsnFunctionAxis(AliRsnFunctionAxis::kEventMult,     8,  0.0, 200.0);
  for (Int_t i = 0; i < 2; i++)
  {
    task[i]->AddAxis(axisMult);
    task[i]->AddAxis(axisPt);
    task[i]->AddAxis(axisEta);
  }

  // define cuts for event selection:
  // this will determine the filling of bins in the "info" histograms
  // and should be computed as additional correction factor in efficiency
  AliRsnCutPrimaryVertex *cutVertex   = new AliRsnCutPrimaryVertex("cutVertex", 3);
  AliRsnCutSet           *cutSetEvent = new AliRsnCutSet("eventCuts");
  cutSetEvent->AddCut(cutVertex);
  cutSetEvent->SetCutScheme("cutVertex");
  for (Int_t i = 0; i < 2; i++)
  {
    task[i]->SetEventCuts(cutSetEvent);
  }

  //
  // *** STEP 0 - All resonances which decay in the specified pairs
  //
  // This step does not need any kind of definition, since
  // its requirement is automatically checked during execution,
  // but to avoid segfaults, it is better to initialize a cut manager.
  //
  AliRsnCutMgr *cutMgrMC_step0 = new AliRsnCutMgr("mc_step0", "");

  //
  // *** STEP 1 - Acceptance
  //
  // Here we add a cut on the pseudorapidity for both tracks
  //
  AliRsnCutStd *cutEta = new AliRsnCutStd("cutEta", AliRsnCutStd::kEta, -0.9, 0.9);

  AliRsnCutMgr *cutMgrMC_step1    = new AliRsnCutMgr("mc_step1", "");
  AliRsnCutSet *cutSetTrack_step1 = new AliRsnCutSet("mc_step1_tracks");
  
  cutSetTrack_step1->AddCut(cutEta);
  cutSetTrack_step1->SetCutScheme("cutEta");
  cutMgrMC_step1   ->SetCutSet(AliRsnCut::kParticle, cutSetTrack_step1);

  //
  // *** STEP 2 - Reconstruction & track quality
  //
  // Use the interface to AliESDtrackCuts
  // and set only the cuts we are interested in
  AliRsnCutESDPrimary *cutQuality = new AliRsnCutESDPrimary("cutCov");
  cutQuality->GetCuts()->SetMaxCovDiagonalElements(2.0, 2.0, 0.5, 0.5, 2.0);
  cutQuality->GetCuts()->SetRequireSigmaToVertex(kTRUE);
  cutQuality->GetCuts()->SetMaxNsigmaToVertex(10000.0);
  cutQuality->GetCuts()->SetRequireTPCRefit(kTRUE);
  cutQuality->GetCuts()->SetAcceptKinkDaughters(kTRUE);
  cutQuality->GetCuts()->SetMinNClustersTPC(50);
  cutQuality->GetCuts()->SetMaxChi2PerClusterTPC(3.5);

  AliRsnCutMgr *cutMgrESD_step2   = new AliRsnCutMgr("esd_step2", "");
  AliRsnCutSet *cutSetTrack_step2 = new AliRsnCutSet("esd_step2_tracks");

  cutSetTrack_step2->AddCut(cutQuality);
  cutSetTrack_step2->SetCutScheme("cutQuality");
  cutMgrESD_step2  ->SetCutSet(AliRsnCut::kParticle, cutSetTrack_step2);

  //
  // *** STEP 3 - Primary tracks
  //
  // Use the interface to AliESDtrackCuts
  // and set only the cuts we are interested in
  // we also disable the cuts we applied before, for clarity
  AliRsnCutESDPrimary *cutESDPrimary = new AliRsnCutESDPrimary("cutESDPrimary");
  cutESDPrimary->GetCuts()->SetMaxCovDiagonalElements(100.0, 100.0, 100.0, 100.0, 100.0);
  cutESDPrimary->GetCuts()->SetRequireSigmaToVertex(kTRUE);
  cutESDPrimary->GetCuts()->SetMaxNsigmaToVertex(3.0);
  cutESDPrimary->GetCuts()->SetRequireTPCRefit(kFALSE);
  cutESDPrimary->GetCuts()->SetAcceptKinkDaughters(kFALSE);
  cutESDPrimary->GetCuts()->SetMinNClustersTPC(0);
  cutESDPrimary->GetCuts()->SetMaxChi2PerClusterTPC(100000000.0);

  AliRsnCutMgr *cutMgrESD_step3   = new AliRsnCutMgr("esd_step3", "");
  AliRsnCutSet *cutSetTrack_step3 = new AliRsnCutSet("esd_step3_tracks");

  cutSetTrack_step3->AddCut(cutESDPrimary);
  cutSetTrack_step3->SetCutScheme("cutESDPrimary");
  cutMgrESD_step3  ->SetCutSet(AliRsnCut::kParticle, cutSetTrack_step3);

  //
  // *** STEP 4 - Two possibilities (depend on the first macro argument)
  //
  // option 1 = Bethe-Bloch cut in 3 sigma (the sigma is one argument)
  // option 2 = realistic Bayesian PID with all detectors
  //
  AliRsnCutMgr *cutMgrESD_step4[2];
  AliRsnCutSet *cutSetTrack_step4[2];
  
  cutMgrESD_step4[0] = new AliRsnCutMgr("esd_step4_nopid", "");
  cutMgrESD_step4[1] = new AliRsnCutMgr("esd_step4_pid", "");
  cutSetTrack_step4[0] = new AliRsnCutSet("esd_step4_tracks_nopid");
  cutSetTrack_step4[1] = new AliRsnCutSet("esd_step4_tracks_pid");
  
  // Bethe-Bloch with kaon mass hypothesis
  AliRsnCutBetheBloch *cutKaonBB = new AliRsnCutBetheBloch("cutKaonBB", 3.0 * sigmaTPC, AliPID::kKaon);
  cutKaonBB->SetCalibConstant(0, 0.76176e-1);
  cutKaonBB->SetCalibConstant(1, 10.632);
  cutKaonBB->SetCalibConstant(2, 0.13279e-4);
  cutKaonBB->SetCalibConstant(3, 1.8631);
  cutKaonBB->SetCalibConstant(4, 1.9479);

  // cuts for realistic PID match
  AliRsnCutStd *cutRealisticPID = new AliRsnCutStd("cutKaonPID", AliRsnCutStd::kRealisticPIDMatch, 0);

  cutSetTrack_step4[0]->AddCut(cutKaonBB);
  cutSetTrack_step4[0]->SetCutScheme("cutKaonBB");

  cutSetTrack_step4[1]->AddCut(cutRealisticPID);
  cutSetTrack_step4[1]->SetCutScheme("cutRealisticPID");

  cutMgrESD_step4[0]->SetCutSet(AliRsnCut::kParticle, cutSetTrack_step4[0]);
  cutMgrESD_step4[1]->SetCutSet(AliRsnCut::kParticle, cutSetTrack_step4[1]);

  // add all steps to the task
  for (Int_t i = 0; i < 2; i++)
  {
    task[i]->AddStepMC (cutMgrMC_step0);
    task[i]->AddStepMC (cutMgrMC_step1);
    task[i]->AddStepESD(cutMgrESD_step2);
    task[i]->AddStepESD(cutMgrESD_step3);
    task[i]->AddStepESD(cutMgrESD_step4[i]);
  }

  // connect containers and finalize
  for (Int_t i = 0; i < 2; i++)
  {
    mgr->AddTask(task[i]);
    mgr->ConnectInput(task[i], 0, mgr->GetCommonInputContainer());

    // create paths for the output in the common file
    Char_t infoPath[500], effPath[500];
    sprintf(infoPath , "%s:PWG2RSNINFO" , AliAnalysisManager::GetCommonFileName());
    sprintf(effPath  , "%s:PWG2RSNEFF%s", AliAnalysisManager::GetCommonFileName(), suf[i].Data());

    // initialize and connect container for the output
    AliAnalysisDataContainer *info = 0x0, *out = 0x0;
    info = mgr->CreateContainer(Form("EffInfo_%s", suf[i].Data()), TList::Class(), AliAnalysisManager::kOutputContainer, infoPath);
    out  = mgr->CreateContainer(Form("EFF_%s", suf[i].Data()), TList::Class(), AliAnalysisManager::kOutputContainer, effPath);

    mgr->ConnectOutput(task[i], 1, info);
    mgr->ConnectOutput(task[i], 2, out);
  }

  return kTRUE;
}
