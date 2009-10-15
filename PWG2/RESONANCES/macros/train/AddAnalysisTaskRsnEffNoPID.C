Bool_t AddAnalysisTaskRsnEffNoPID
(
  const char *outFile = "eff_nopid.root",    // output file name
)
{
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // create task
  AliRsnAnalysisEffSE *task = new AliRsnAnalysisEffSE("EffNoPID");

  // set prior probabilities for PID
  task->SetPriorProbability(AliPID::kElectron, 0.02);
  task->SetPriorProbability(AliPID::kMuon,     0.02);
  task->SetPriorProbability(AliPID::kPion,     0.83);
  task->SetPriorProbability(AliPID::kKaon,     0.07);
  task->SetPriorProbability(AliPID::kProton,   0.06);
  task->DumpPriors();

  // pair definition
  AliRsnPairDef *pairDef1 = new AliRsnPairDef('+', AliPID::kKaon, '-', AliPID::kKaon, 333);
  AliRsnPairDef *pairDef2 = new AliRsnPairDef('-', AliPID::kKaon, '+', AliPID::kPion, 313);
  AliRsnPairDef *pairDef3 = new AliRsnPairDef('+', AliPID::kKaon, '-', AliPID::kPion, 313);
  task->AddPairDef(pairDef1);
  task->AddPairDef(pairDef2);
  task->AddPairDef(pairDef3);

  // axis definition
  AliRsnFunctionAxis *axisPt   = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairPt,      100,  0.0,  10.0);
  AliRsnFunctionAxis *axisEta  = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairEta,      10, -1.0,   1.0);
  AliRsnFunctionAxis *axisMult = new AliRsnFunctionAxis(AliRsnFunctionAxis::kEventMult,   500,  0.0, 500.0);
  task->AddAxis(axisMult);
  task->AddAxis(axisPt);
  task->AddAxis(axisEta);

  // setup cuts for events (good primary vertex)
  AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 3);
  AliRsnCutSet *cutSetEvent = new AliRsnCutSet("eventCuts");
  cutSetEvent->AddCut(cutVertex);
  cutSetEvent->SetCutScheme("cutVertex");
  task->SetEventCuts(cutSetEvent);

  // *** STEP 0 - All resonances

  AliRsnCutMgr *cutMgr0      = new AliRsnCutMgr("step0", "");
  AliRsnCutSet *cutSetTrack0 = new AliRsnCutSet("step0_tracks");
  AliRsnCutSet *cutSetPair0  = new AliRsnCutSet("step0_pairs");

  cutMgr0->SetCutSet(AliRsnCut::kParticle, cutSetTrack0);
  cutMgr0->SetCutSet(AliRsnCut::kPair    , cutSetPair0 );

  task->AddStepMC(cutMgr0);

  // *** STEP 1 - Acceptance

  AliRsnCutMgr *cutMgr1      = new AliRsnCutMgr("step1", "");
  AliRsnCutSet *cutSetTrack1 = new AliRsnCutSet("step1_tracks");
  AliRsnCutSet *cutSetPair1  = new AliRsnCutSet("step1_pairs");

  AliRsnCutStd *cutEta = new AliRsnCutStd("cutEta", AliRsnCutStd::kEta, -0.9, 0.9);

  cutSetTrack1->AddCut(cutEta);
  cutSetTrack1->SetCutScheme("cutEta");

  cutMgr1->SetCutSet(AliRsnCut::kParticle, cutSetTrack1);
  cutMgr1->SetCutSet(AliRsnCut::kPair    , cutSetPair1 );

  task->AddStepMC(cutMgr1);

  // *** STEP 2 - Reconstruction & quality

  AliRsnCutMgr *cutMgr2      = new AliRsnCutMgr("step2", "");
  AliRsnCutSet *cutSetTrack2 = new AliRsnCutSet("step2_tracks");
  AliRsnCutSet *cutSetPair2  = new AliRsnCutSet("step2_pairs");

  // cuts for tracks:
  // -- primary track quality
  AliRsnCutESDPrimary *cutESDPrimary = new AliRsnCutESDPrimary("cutESDPrimary");
  cutESDPrimary->GetCuts()->SetMaxCovDiagonalElements(2.0, 2.0, 0.5, 0.5, 2.0);
  cutESDPrimary->GetCuts()->SetRequireSigmaToVertex(kTRUE);
  cutESDPrimary->GetCuts()->SetMaxNsigmaToVertex(4.0);
  cutESDPrimary->GetCuts()->SetRequireTPCRefit(kTRUE);
  cutESDPrimary->GetCuts()->SetAcceptKinkDaughters(kFALSE);
  cutESDPrimary->GetCuts()->SetMinNClustersTPC(50);
  cutESDPrimary->GetCuts()->SetMaxChi2PerClusterTPC(3.5);

  cutSetTrack2->AddCut(cutESDPrimary);
  cutSetTrack2->SetCutScheme("cutESDPrimary");

  cutMgr2->SetCutSet(AliRsnCut::kParticle, cutSetTrack2);
  cutMgr2->SetCutSet(AliRsnCut::kPair    , cutSetPair2 );

  task->AddStepESD(cutMgr2);

  // *** STEP 3 - Bethe-Bloch cut (0.2)

  AliRsnCutMgr *cutMgr3      = new AliRsnCutMgr("step3", "");
  AliRsnCutSet *cutSetTrack3 = new AliRsnCutSet("step3_tracks");
  AliRsnCutSet *cutSetPair3  = new AliRsnCutSet("step3_pairs");

  // cuts for tracks:
  // -- Bethe-Bloch with kaon mass hypothesis
  AliRsnCutBetheBloch *cutKaonBB = new AliRsnCutBetheBloch("cutKaonBB", 0.2, AliPID::kKaon);
  cutKaonBB->SetCalibConstant(0, 0.76176e-1);
  cutKaonBB->SetCalibConstant(1, 10.632);
  cutKaonBB->SetCalibConstant(2, 0.13279e-4);
  cutKaonBB->SetCalibConstant(3, 1.8631);
  cutKaonBB->SetCalibConstant(4, 1.9479);

  cutSetTrack3->AddCut(cutKaonBB);
  cutSetTrack3->SetCutScheme("cutKaonBB");
  //AliLog::SetClassDebugLevel("AliRsnCut", AliLog::kDebug+3);
  //AliLog::SetClassDebugLevel("AliRsnCutStd", AliLog::kDebug+3);
  //AliLog::SetClassDebugLevel("AliRsnCutBetheBloch", AliLog::kDebug+3);

  cutMgr3->SetCutSet(AliRsnCut::kParticle, cutSetTrack3);
  cutMgr3->SetCutSet(AliRsnCut::kPair    , cutSetPair3 );

  task->AddStepESD(cutMgr3);

  // add the task to manager
  mgr->AddTask(task);

  // connect input container according to source choice
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // initialize and connect container for the output
  AliAnalysisDataContainer *outputInfo = mgr->CreateContainer("EffNoPIDInfo", TList::Class(), AliAnalysisManager::kOutputContainer, "info.root");
  AliAnalysisDataContainer *out = mgr->CreateContainer("EFF_NOPID", TList::Class(), AliAnalysisManager::kOutputContainer, outFile);
  mgr->ConnectOutput(task, 1, outputInfo);
  mgr->ConnectOutput(task, 2, out);

  return kTRUE;
}
