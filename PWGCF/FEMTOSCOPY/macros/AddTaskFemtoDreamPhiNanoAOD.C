AliAnalysisTaskSE *AddTaskFemtoDreamPhiNanoAOD(bool isMC = false,
                                               TString CentEst = "kInt7",
                                               const char *cutVariation = "0") {
  TString suffix = TString::Format("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);
  evtCuts->SetSphericityCuts(0.7, 1);

  // Proton cuts
  const float ProtonPtlow = 0.425;
  const float ProtonPtup = 0.575;
  const float ProtonEtaLow = 0.75;
  const float ProtonEtaUp = 0.85;
  const float ProtonNsigmaLow = 2.5;
  const float ProtonNsigmaUp = 3.5;
  const float ProtonNClsLow = 70;
  const float ProtonNClsUp = 90;

  AliFemtoDreamTrackCuts *TrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetCutCharge(-1);

  if (suffix != "0" && suffix != "999") {
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
  }





  // Kaon cuts
  const float KaonPtlow = 0.075;
  const float KaonPtup = 0.225;
  const float KaonEtaLow = 0.75;
  const float KaonEtaUp = 0.85;
  const float KaonNsigmaLow = 4.25;
  const float KaonNsigmaUp = 5.75;
  const float KaonNClsLow = 70;
  const float KaonNClsUp = 90;
  const float SpheriLow = 0.665;
  const float SpheriUp = 0.735;

  AliFemtoDreamTrackCuts *TrackPosKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackPosKaonCuts->SetCutCharge(1);
  TrackPosKaonCuts->SetFilterBit(128);

  AliFemtoDreamTrackCuts *TrackNegKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackNegKaonCuts->SetCutCharge(-1);
  TrackNegKaonCuts->SetFilterBit(128);

  TrackPosKaonCuts->SetDCAVtxZ(0.4);
  TrackNegKaonCuts->SetDCAVtxZ(0.4);
  TrackPosKaonCuts->SetDCAVtxXY(0.8);
  TrackNegKaonCuts->SetDCAVtxXY(0.8);

  if (suffix != "0" && suffix != "999") {
    TrackPosKaonCuts->SetMinimalBooking(true);
    TrackNegKaonCuts->SetMinimalBooking(true);
  }



  AliFemtoDreamv0Cuts *TrackCutsPhi = new AliFemtoDreamv0Cuts();
  TrackCutsPhi->SetIsMonteCarlo(isMC);
  TrackCutsPhi->SetAxisInvMassPlots(400, 0.95, 2);
  TrackCutsPhi->SetCutInvMass(0.008);
  AliFemtoDreamTrackCuts *dummyCutsPos = new AliFemtoDreamTrackCuts();
  dummyCutsPos->SetIsMonteCarlo(isMC);
  AliFemtoDreamTrackCuts *dummyCutsNeg = new AliFemtoDreamTrackCuts();
  dummyCutsNeg->SetIsMonteCarlo(isMC);
  TrackCutsPhi->SetPosDaugterTrackCuts(dummyCutsPos);
  TrackCutsPhi->SetNegDaugterTrackCuts(dummyCutsNeg);
  TrackCutsPhi->SetPDGCodePosDaug(321);
  TrackCutsPhi->SetPDGCodeNegDaug(321);
  TrackCutsPhi->SetPDGCodev0(333);

  double Phimass = TDatabasePDG::Instance()->GetParticle(333)->Mass();

  if (suffix != "0") {
    TrackCutsPhi->SetMinimalBooking(true);
  }


  // Now we define stuff we want for our Particle collection
  // Thanks, CINT - will not compile due to an illegal constructor
  // std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  // First we need to tell him about the particles we mix, from the
  // PDG code the mass is obtained.
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);  // 0 protons REC MC
  PDGParticles.push_back(2212);  // 1 antiprot REC MC
  PDGParticles.push_back(333);   // 2 v0 particle REC MC
  PDGParticles.push_back(2212);  // 3 proton truth MC TRUTH
  PDGParticles.push_back(2212);  // 4 antiprot truth MC TRUTH
  PDGParticles.push_back(333);   // 5 True Phi inv mass MC TRUTH
  PDGParticles.push_back(333);   // 6 Fake Phi inv mass (combi. BG) MC TRUTHn
  PDGParticles.push_back(321);   // 7 Kaon
  PDGParticles.push_back(321);   // 8 AntiKaon
  PDGParticles.push_back(321);   // 9 Kaon inv mass phi
  PDGParticles.push_back(321);   // 10 AngiKaon inv mass phi


  // We need to set the ZVtx bins
  std::vector<float> ZVtxBins;
  ZVtxBins.push_back(-10);
  ZVtxBins.push_back(-8);
  ZVtxBins.push_back(-6);
  ZVtxBins.push_back(-4);
  ZVtxBins.push_back(-2);
  ZVtxBins.push_back(0);
  ZVtxBins.push_back(2);
  ZVtxBins.push_back(4);
  ZVtxBins.push_back(6);
  ZVtxBins.push_back(8);
  ZVtxBins.push_back(10);
  // The Multiplicity bins are set here
  std::vector<int> MultBins;
  MultBins.push_back(0);
  MultBins.push_back(4);
  MultBins.push_back(8);
  MultBins.push_back(12);
  MultBins.push_back(16);
  MultBins.push_back(20);
  MultBins.push_back(24);
  MultBins.push_back(28);
  MultBins.push_back(32);
  MultBins.push_back(36);
  MultBins.push_back(40);
  MultBins.push_back(44);
  MultBins.push_back(48);
  MultBins.push_back(52);
  MultBins.push_back(56);
  MultBins.push_back(60);
  MultBins.push_back(64);
  MultBins.push_back(68);
  MultBins.push_back(72);
  MultBins.push_back(76);
  MultBins.push_back(80);
  MultBins.push_back(84);
  MultBins.push_back(88);
  MultBins.push_back(92);
  MultBins.push_back(96);
  MultBins.push_back(100);

  // Number of bins
  std::vector<int> NBins;
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  std::vector<float> kMin;
  // minimum k* value
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  // maximum k* value
  std::vector<float> kMax;
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  // pair QA extended
  std::vector<int> pairQA;
  //  pairQA.push_back(11); //pp
  //  pairQA.push_back(0); //pap
  //  pairQA.push_back(12); //pphi
  //  pairQA.push_back(11); //apap
  //  pairQA.push_back(12); //apphi
  //  pairQA.push_back(22); //phiphi

  for (int i = 0; i < (1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10 + 11); i++) {
    NBins.push_back(750);
    kMin.push_back(0.);
    kMax.push_back(3.);
    pairQA.push_back(0);
  }

  pairQA[0] = 11;   // pp
  pairQA[1] = 11;   // pap
  pairQA[2] = 12;   // pphi
  pairQA[3] = 0;
  pairQA[4] = 0;
  pairQA[5] = 0;
  pairQA[6] = 0;
  pairQA[7] = 11;   // KK
  pairQA[8] = 11;   // KaK
  pairQA[9] = 11;
  pairQA[10] = 11;


  pairQA[11] = 11;   // apap
  pairQA[12] = 12;   // apphi
  pairQA[13] = 0;
  pairQA[14] = 0;
  pairQA[15] = 0;
  pairQA[16] = 0;
  pairQA[17] = 11;   // KK
  pairQA[18] = 11;   // KaK
  pairQA[19] = 11;
  pairQA[20] = 11;

  pairQA[21] = 0;  // phiphi

  if (isMC) {
   for (int i=22;i<66;i++)  {
       pairQA[i] = 0;
      }
  }

  AliFemtoDreamCollConfig *config =
      new AliFemtoDreamCollConfig("Femto", "Femto");
//  config->SetPtQA(true);
//  config->SetMassQA(true);
//  config->SetmTBinning(true);
//  config->SetdPhidEtaPlots(true);
  config->SetExtendedQAPairs(pairQA);
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetMultBinning(true);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetUseEventMixing(true);
  config->SetMixingDepth(30);
//  config->SetPhiEtaBinnign(true);
  //config->SetMomentumResolution(true);

//  if (suffix == "0") {
//      config->SetMassQA(true);
//  }
//  if (suffix == "1") {
//      config->SetMassQA(true);
//  }
//  if (suffix == "2") {
//      config->SetMassQA(true);
//  }


  //-------MIXED EVENTS---------------------------
  // config->SetUseEventMixing(true);
  // config->SetMixingDepth(10);

  //-------ROTATED SAMPLE-------------------------
  //  config->SetUseEventMixing(false);
  //  config->SetUsePhiSpinning(true);
  //  config->SetControlMethod(AliFemtoDreamCollConfig::kPhiSpin);
  //  config->SetSpinningDepth(1);

  //-------COMMON ANCESTORS-------------------------

//  if (isMC) {
//    config->SetAncestors(true);
//    config->GetDoAncestorsPlots();
//  }

  if (isMC) {
    config->SetMomentumResolution(true);
  } else {
    std::cout << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
  }

  /*
  //This is just to show off what would be possible in case you are interested,
  don't be confused by this at the beginning
  //you can just ignore it!
  if (false) {
    config->SetkTBinning(false);
    config->SetmTBinning(false);
    config->SetkTCentralityBinning(false);
    std::vector<float> centBins;
    centBins.push_back(20);
    centBins.push_back(40);
    centBins.push_back(90);
    config->SetCentBins(centBins);
    config->SetZBins(ZVtxBins);

    if (isMC) {
      config->SetMomentumResolution(true);
    } else {
      std::cout << "You are trying to request the Momentum Resolution without MC
  Info; fix it wont work! \n";
    }
    if (isMC) {
      config->SetPhiEtaBinnign(true);
    } else {
      std::cout << "You are trying to request the Eta Phi Plots without MC Info;
  fix it wont work! \n";
    }
  }
   */
  // now we create the task
  AliAnalysisTaskNanoAODFemtoDreamPhi *task =
      new AliAnalysisTaskNanoAODFemtoDreamPhi(
          "AliAnalysisTaskNanoAODFemtoDreamPhi", isMC);
  // THIS IS VERY IMPORTANT ELSE YOU DONT PROCESS ANY EVENTS
  // kINT7 == Minimum bias
  // kHighMultV0 high multiplicity triggered by the V0 detector
  if (CentEst == "kInt7") {
    task->SetTrigger(AliVEvent::kINT7);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (CentEst == "kHM") {
    task->SetTrigger(AliVEvent::kHighMultV0);
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
  } else {
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "Centrality Estimator not set, fix it else your Results will "
                 "be empty!"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
  }

  // Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetPosKaonCuts(TrackPosKaonCuts);
  task->SetNegKaonCuts(TrackNegKaonCuts);
  task->SetCollectionConfig(config);
  task->SetPhiCuts(TrackCutsPhi);
  task->SetOEventMixing(false);
  task->SetMCTruth(false);

//  task->SetCutWindowMCTRUTH(Phimass - 0.008, Phimass + 0.008);
//  task->SetmixREC(true);
//  task->SetMCTruth(true);
//  task->SetmixTRUTHREAL(true);
//  task->SetmixTRUTFAKE(true);
//  task->SetpTprotonMCTRUTH(0.5);
//  task->SetetaprotonMCTRUTH(0.8);
//  task->SetpTkaonMCTRUTH(0.15);
//  task->SetetakaonMCTRUTH(0.8);


//  if (suffix == "0") {
//      task->SetCutWindowMCTRUTH(Phimass - 0.008, Phimass + 0.008);
//      config->SetPtQA(true);
//      config->SetMassQA(true);
//      config->SetmTBinning(true);
//      config->SetdPhidEtaPlots(true);
//      config->SetExtendedQAPairs(pairQA);
//      task->SetmixREC(true);
//      task->SetmixTRUTHREAL(true);
//      task->SetmixTRUTFAKE(true);
//  }
//    if (suffix == "1") {
//        task->SetCutWindowMCTRUTH(0.995, 1.011);
//        TrackCutsPhi->SetMinimalBooking(true);
//        task->SetmixREC(true);
//        task->SetmixTRUTHREAL(true);
//        task->SetmixTRUTFAKE(true);
//    }
//    if (suffix == "2") {
//        task->SetCutWindowMCTRUTH(1.028, 1.044);
//        TrackCutsPhi->SetMinimalBooking(true);
//        task->SetmixREC(true);
//        task->SetmixTRUTHREAL(true);
//        task->SetmixTRUTFAKE(true);
//    }



  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);

  AliAnalysisDataContainer *coutputQA;
  TString addon = "";
  if (CentEst == "kInt7") {
    addon += "MB";
  } else if (CentEst == "kHM") {
    addon += "HM";
  }
  TString QAName = Form("%sPhiResults%s", addon.Data(), suffix.Data());
  coutputQA = mgr->CreateContainer(QAName.Data(), TList::Class(),
                                   AliAnalysisManager::kOutputContainer,
                                   Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}
