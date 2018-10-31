AliAnalysisTaskSE *AddTaskSigma0Femto(bool isMC = false,
                                      bool isHeavyIon = false,
                                      bool MomRes = false,
                                      bool etaPhiPlotsAtTPCRadii = false,
                                      TString trigger = "kINT7",
                                      const char *cutVariation = "0") {
  TString suffix;
  suffix.Form("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskSigma0Run2()", "No analysis manager found.");
    return 0x0;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton;
  cutnumberPhoton = "00200008400000002280920000";
  TString cutnumberEvent = "00000000";
  TString periodNameV0Reader = "";
  Bool_t enableV0findingEffi = kFALSE;
  Bool_t runLightOutput = kFALSE;
  if (suffix != "0" && suffix != "999") {
    runLightOutput = kTRUE;
  }
  if (suffix == "28") {
    // eta < 0.8
    cutnumberPhoton = "0d200008400000002280920000";
  } else if (suffix == "29") {
    // eta < 0.75
    cutnumberPhoton = "04200008400000002280920000";
  } else if (suffix == "30") {
    // r in 0 - 180 cm
    cutnumberPhoton = "00000008400000002280920000";
  } else if (suffix == "31") {
    // r in 10 - 180 cm
    cutnumberPhoton = "00500008400000002280920000";
  } else if (suffix == "32") {
    // single pT > 0, gammapT > 0
    cutnumberPhoton = "00200078400000002280920000";
  } else if (suffix == "33") {
    // single pT > 0.150, gammapT > 0.02
    cutnumberPhoton = "00200028400000002280920000";
  } else if (suffix == "34") {
    // single pT > 0.050, gammapT > 0.150
    cutnumberPhoton = "002000a8400000002280920000";
  } else if (suffix == "35") {
    // TPC cluster, findable > 0.6
    cutnumberPhoton = "00200009400000002280920000";
  } else if (suffix == "36") {
    // TPC PID -10,10
    cutnumberPhoton = "00200008000000002280920000";
  } else if (suffix == "37") {
    // TPC PID -3,3
    cutnumberPhoton = "00200008a00000002280920000";
  } else if (suffix == "38") {
    // TPC pi, p, k rejection -1,1
    cutnumberPhoton = "00200008400020002280920000";
  } else if (suffix == "39") {
    // 1-D Qt cut, qt < 0.1
    cutnumberPhoton = "00200008400000001280920000";
  } else if (suffix == "40") {
    // 2-D Qt cut, qt < 0.02
    cutnumberPhoton = "00200008400000006280920000";
  } else if (suffix == "41") {
    // chi2 < 100
    cutnumberPhoton = "00200008400000002080920000";
  } else if (suffix == "42") {
    // chi2 < 10
    cutnumberPhoton = "00200008400000002780920000";
  } else if (suffix == "43") {
    // psiPair < 0.2, 1-D
    cutnumberPhoton = "00200008400000002240920000";
  } else if (suffix == "44") {
    // psiPair < 0.1, 2-D
    cutnumberPhoton = "00200008400000002250920000";
  } else if (suffix == "45") {
    // cosPA < 0.98
    cutnumberPhoton = "00200008400000002280820000";
  } else if (suffix == "46") {
    // cosPA < 0.995
    cutnumberPhoton = "00200008400000002280a20000";
  } else if (suffix == "47") {
    // shared electron cut, no photon QA selection
    cutnumberPhoton = "00200008400000002280910000";
  } else if (suffix == "48") {
    // DCA_R < 5
    cutnumberPhoton = "00200008400000002280920200";
  } else if (suffix == "49") {
    // DCA_Z < 5
    cutnumberPhoton = "00200008400000002280920020";
  }
  // 50 - pile up cut, done in the task

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName =
      Form("V0ReaderV1_%s_%s", cutnumberEvent.Data(), cutnumberPhoton.Data());
  AliConvEventCuts *fEventCuts = NULL;

  if (!(AliV0ReaderV1 *)mgr->GetTask(V0ReaderName.Data())) {
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    if (periodNameV0Reader.CompareTo("") != 0)
      fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);  // AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return nullptr;
    }

    if (cutnumberEvent != "") {
      fEventCuts =
          new AliConvEventCuts(cutnumberEvent.Data(), cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      fEventCuts->SetLightOutput(runLightOutput);
      if (periodNameV0Reader.CompareTo("") != 0)
        fEventCuts->SetPeriodEnum(periodNameV0Reader);
      fV0ReaderV1->SetEventCuts(fEventCuts);
      fEventCuts->SetFillCutHistograms("", kTRUE);
    }

    // Set AnalysisCut Number
    AliConversionPhotonCuts *fCuts = NULL;
    if (cutnumberPhoton != "") {
      fCuts = new AliConversionPhotonCuts(cutnumberPhoton.Data(),
                                          cutnumberPhoton.Data());
      fCuts->SetPreSelectionCutFlag(kTRUE);
      fCuts->SetIsHeavyIon(isHeavyIon);
      fCuts->SetV0ReaderName(V0ReaderName);
      fCuts->SetLightOutput(runLightOutput);
      if (fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())) {
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("", kTRUE);
      }
    }

    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kFatal);

    // connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1, 0, cinput);
  }

  //========= Init subtasks and start analyis ============================
  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  TrackCuts->SetCheckFilterBit(false);
  TrackCuts->SetCheckESDFiltering(true);
  TrackCuts->SetDCAReCalculation(false);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetCheckFilterBit(false);
  AntiTrackCuts->SetCheckESDFiltering(true);
  AntiTrackCuts->SetDCAReCalculation(false);
  AntiTrackCuts->SetCutCharge(-1);

  if (suffix != "0" && suffix != "999") {
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
  }
  if (suffix == "51") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);
  } else if (suffix == "52") {
    TrackCuts->SetPtRange(0.5, 4.05);
    AntiTrackCuts->SetPtRange(0.5, 4.05);
  } else if (suffix == "53") {
    TrackCuts->SetEtaRange(-0.7, 0.7);
    AntiTrackCuts->SetEtaRange(-0.7, 0.7);
  } else if (suffix == "54") {
    TrackCuts->SetEtaRange(-0.9, 0.9);
    AntiTrackCuts->SetEtaRange(-0.9, 0.9);
  } else if (suffix == "55") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, 2);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2);
  } else if (suffix == "56") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, 5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 5);
  } else if (suffix == "57") {
    TrackCuts->SetFilterBit(96);
    AntiTrackCuts->SetFilterBit(96);
  } else if (suffix == "58") {
    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);
  } else if (suffix == "59") {
    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);
  }

  AliSigma0V0Cuts *v0Cuts = AliSigma0V0Cuts::LambdaCuts();
  v0Cuts->SetIsMC(isMC);
  v0Cuts->SetPID(3122);
  v0Cuts->SetPosPID(AliPID::kProton, 2212);
  v0Cuts->SetNegPID(AliPID::kPion, -211);
  AliSigma0V0Cuts *antiv0Cuts = AliSigma0V0Cuts::LambdaCuts();
  antiv0Cuts->SetIsMC(isMC);
  antiv0Cuts->SetPID(-3122);
  antiv0Cuts->SetPosPID(AliPID::kPion, 211);
  antiv0Cuts->SetNegPID(AliPID::kProton, -2212);

  if (suffix != "0") {
    v0Cuts->SetLightweight(true);
    antiv0Cuts->SetLightweight(true);
  }
  if (suffix == "1") {
    v0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    antiv0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
  } else if (suffix == "2") {
    v0Cuts->SetV0PtMin(0.24);
    antiv0Cuts->SetV0PtMin(0.24);
  } else if (suffix == "3") {
    v0Cuts->SetV0PtMin(0.36);
    antiv0Cuts->SetV0PtMin(0.36);
  } else if (suffix == "4") {
    v0Cuts->SetV0CosPAMin(0.995);
    antiv0Cuts->SetV0CosPAMin(0.995);
  } else if (suffix == "5") {
    v0Cuts->SetV0CosPAMin(0.98);
    antiv0Cuts->SetV0CosPAMin(0.98);
  } else if (suffix == "6") {
    v0Cuts->SetPIDnSigma(3);
    antiv0Cuts->SetPIDnSigma(3);
  } else if (suffix == "7") {
    v0Cuts->SetPIDnSigma(3);
    antiv0Cuts->SetPIDnSigma(3);
    v0Cuts->SetArmenterosCut(0, 1, -1, 1);
    antiv0Cuts->SetArmenterosCut(0, 1, -1, 1);
  } else if (suffix == "8") {
    v0Cuts->SetPIDnSigma(6);
    antiv0Cuts->SetPIDnSigma(6);
  } else if (suffix == "9") {
    v0Cuts->SetPIDnSigma(6);
    antiv0Cuts->SetPIDnSigma(6);
    v0Cuts->SetArmenterosCut(0, 1, -1, 1);
    antiv0Cuts->SetArmenterosCut(0, 1, -1, 1);
  } else if (suffix == "10") {
    v0Cuts->SetArmenterosCut(0, 1, -1, 1);
    antiv0Cuts->SetArmenterosCut(0, 1, -1, 1);
  } else if (suffix == "11") {
    v0Cuts->SetTPCclusterMin(80);
    antiv0Cuts->SetTPCclusterMin(80);
  } else if (suffix == "12") {
    v0Cuts->SetTPCclusterMin(60);
    antiv0Cuts->SetTPCclusterMin(60);
  } else if (suffix == "13") {
    v0Cuts->SetEtaMax(0.8);
    antiv0Cuts->SetEtaMax(0.8);
  } else if (suffix == "14") {
    v0Cuts->SetEtaMax(0.75);
    antiv0Cuts->SetEtaMax(0.75);
  } else if (suffix == "15") {
    v0Cuts->SetDaughterDCAMax(1.2);
    antiv0Cuts->SetDaughterDCAMax(1.2);
  } else if (suffix == "16") {
    v0Cuts->SetDaughterDCAMax(0.9);
    antiv0Cuts->SetDaughterDCAMax(0.9);
  } else if (suffix == "17") {
    v0Cuts->SetDaughterDCAtoPV(0.06);
    antiv0Cuts->SetDaughterDCAtoPV(0.06);
  } else if (suffix == "18") {
    v0Cuts->SetDaughterDCAtoPV(0.04);
    antiv0Cuts->SetDaughterDCAtoPV(0.04);
  } else if (suffix == "19") {
    v0Cuts->SetK0Rejection(0.48, 0.515);
    antiv0Cuts->SetK0Rejection(0.48, 0.515);
  } else if (suffix == "20") {
    v0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    antiv0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
  } else if (suffix == "21") {
    v0Cuts->SetV0RadiusMax(120.f);
    antiv0Cuts->SetV0RadiusMax(120.f);
    v0Cuts->SetV0DecayVertexMax(120.f);
    antiv0Cuts->SetV0DecayVertexMax(120.f);
  } else if (suffix == "22") {
    v0Cuts->SetV0RadiusMax(80.f);
    antiv0Cuts->SetV0RadiusMax(80.f);
  } else if (suffix == "23") {
    v0Cuts->SetV0RadiusMin(0.);
    antiv0Cuts->SetV0RadiusMin(0.);
  } else if (suffix == "24") {
    v0Cuts->SetV0RadiusMin(5.);
    antiv0Cuts->SetV0RadiusMin(5.);
  } else if (suffix == "25") {
    v0Cuts->SetV0RadiusMin(5.);
    antiv0Cuts->SetV0RadiusMin(5.);
  } else if (suffix == "26") {
    v0Cuts->SetV0DecayVertexMax(80.f);
    antiv0Cuts->SetV0DecayVertexMax(80.f);
  } else if (suffix == "27") {
    v0Cuts->SetV0DecayVertexMax(120.f);
    antiv0Cuts->SetV0DecayVertexMax(120.f);
  }
  if (suffix == "999") {
    v0Cuts->SetCheckCutsMC(true);
    antiv0Cuts->SetCheckCutsMC(true);
    v0Cuts->SetLightweight(false);
    antiv0Cuts->SetLightweight(false);
  }

  AliSigma0PhotonMotherCuts *sigmaCuts =
      AliSigma0PhotonMotherCuts::DefaultCuts();
  sigmaCuts->SetIsMC(isMC);
  sigmaCuts->SetPDG(3212, 3122, 22);
  sigmaCuts->SetLambdaCuts(v0Cuts);
  sigmaCuts->SetV0ReaderName(V0ReaderName.Data());
  if (suffix != "0" && suffix != "999") {
    sigmaCuts->SetLightweight(true);
  }

  AliSigma0PhotonMotherCuts *antiSigmaCuts =
      AliSigma0PhotonMotherCuts::DefaultCuts();
  antiSigmaCuts->SetIsMC(isMC);
  antiSigmaCuts->SetPDG(-3212, -3122, 22);
  antiSigmaCuts->SetLambdaCuts(antiv0Cuts);
  antiSigmaCuts->SetV0ReaderName(V0ReaderName.Data());
  if (suffix != "0" && suffix != "999") {
    antiSigmaCuts->SetLightweight(true);
  }
  if (suffix == "60") {
    sigmaCuts->SetPhotonMaxPt(1);
    antiSigmaCuts->SetPhotonMaxPt(1);
  } else if (suffix == "61") {
    sigmaCuts->SetPhotonMaxPt(999);
    antiSigmaCuts->SetPhotonMaxPt(999);
  }

  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3212);
  PDGParticles.push_back(3212);
  PDGParticles.push_back(3212);
  PDGParticles.push_back(3212);
  PDGParticles.push_back(3212);
  PDGParticles.push_back(3212);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(22);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(22);

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

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  for (int i = 0; i < 78; ++i) {
    NBins.push_back(750);
    kMin.push_back(0.);
    kMax.push_back(3.);
  }

  AliFemtoDreamCollConfig *config =
      new AliFemtoDreamCollConfig("Femto", "Femto");
  std::vector<int> MultBins;
  if (trigger == "kHighMultV0") {
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
    config->SetMultBins(MultBins);
  } else {
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
    MultBins.push_back(60);
    MultBins.push_back(80);
    config->SetMultBins(MultBins);
  }
  config->SetMultBinning(true);

  config->SetZBins(ZVtxBins);
  if (MomRes) {
    if (isMC) {
      config->SetMomentumResolution(true);
    } else {
      std::cout << "You are trying to request the Momentum Resolution without "
                   "MC Info; fix it wont work! \n";
    }
  }
  if (etaPhiPlotsAtTPCRadii) {
    if (isMC) {
      config->SetPhiEtaBinnign(true);
    } else {
      std::cout << "You are trying to request the Eta Phi Plots without MC "
                   "Info; fix it wont work! \n";
    }
  }
  config->SetdPhidEtaPlots(etaPhiPlotsAtTPCRadii);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(10);
  config->SetUseEventMixing(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kSPD);
  if (suffix != "0") {
    config->SetMinimalBookingME(true);
  }

  AliAnalysisTaskSigma0Femto *task =
      new AliAnalysisTaskSigma0Femto("AnalysisTaskSigma0Femto");
  if (trigger == "kINT7") {
    task->SetTrigger(AliVEvent::kINT7);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
  } else if (trigger == "kHighMultV0") {
    task->SetTrigger(AliVEvent::kHighMultV0);
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
  }
  task->SetV0ReaderName(V0ReaderName.Data());
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetV0Cuts(v0Cuts);
  task->SetAntiV0Cuts(antiv0Cuts);
  task->SetSigmaCuts(sigmaCuts);
  task->SetAntiSigmaCuts(antiSigmaCuts);
  task->SetCollectionConfig(config);

  if (suffix != "0" && suffix != "999") {
    task->SetLightweight(true);
  }
  if (suffix == "50") {
    task->SetPhotonLegPileUpCut(true);
  }

  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":Sigma0_Femto_";
  if (trigger == "kHighMultV0") containerName += "HighMultV0_";
  containerName += suffix;

  TString name = "histo_";
  if (trigger == "kHighMultV0") name += "HighMultV0_";
  name += suffix;
  AliAnalysisDataContainer *cOutputList = mgr->CreateContainer(
      name, TList::Class(), AliAnalysisManager::kOutputContainer,
      containerName.Data());

  name = "femto_";
  if (trigger == "kHighMultV0") name += "HighMultV0_";
  name += suffix;
  AliAnalysisDataContainer *cFemtoOutputList = mgr->CreateContainer(
      name, TList::Class(), AliAnalysisManager::kOutputContainer,
      containerName.Data());

  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, cOutputList);
  mgr->ConnectOutput(task, 2, cFemtoOutputList);

  return task;
}
