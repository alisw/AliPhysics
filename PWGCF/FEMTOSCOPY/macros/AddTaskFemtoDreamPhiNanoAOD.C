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
  evtCuts->SetSphericityCuts(0.7,1);

//    if (suffix == "1") {
//        evtCuts->SetSpherocityCuts(0.5,1);
//    }

//    if (suffix == "2") {
//        evtCuts->SetSpherocityCuts(0.6,1);
//    }

//    if (suffix == "3") {
//        evtCuts->SetSpherocityCuts(0.7,1);
//    }

//    if (suffix == "4") {
//        evtCuts->SetSpherocityCuts(0.8,1);
//    }

//    if (suffix == "5") {
//        evtCuts->SetSpherocityCuts(0.9,1);
//    }

  AliFemtoDreamTrackCuts *TrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetCutCharge(-1);

  if (suffix != "0") {
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
  }

  AliFemtoDreamTrackCuts *TrackPosKaonCuts = AliFemtoDreamTrackCuts::PrimKaonCuts(isMC);
  TrackPosKaonCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *TrackNegKaonCuts = AliFemtoDreamTrackCuts::PrimKaonCuts(isMC);
  TrackNegKaonCuts->SetCutCharge(-1);

  if (suffix != "0") {
    TrackPosKaonCuts->SetMinimalBooking(true);
    TrackNegKaonCuts->SetMinimalBooking(true);
  }

  if (suffix == "1") {
      TrackPosKaonCuts->SetPtRange(0.15,1.3);
      TrackNegKaonCuts->SetPtRange(0.15,1.3);
  }

  if (suffix == "2") {
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
  }

  if (suffix == "3") {
      TrackPosKaonCuts->SetPtRange(0.15,1.5);
      TrackNegKaonCuts->SetPtRange(0.15,1.5);
  }

  if (suffix == "4") {
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
  }

  if (suffix == "5") {
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPtExclusion(0.2,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.2,0.4);
  }

  if (suffix == "6") {
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.5);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.5);
  }

  if (suffix == "7") {
      TrackPosKaonCuts->SetPtRange(0.15,1.5);
      TrackNegKaonCuts->SetPtRange(0.15,1.5);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
  }

  if (suffix == "8") {
      TrackPosKaonCuts->SetPtRange(0.15,1.5);
      TrackNegKaonCuts->SetPtRange(0.15,1.5);
      TrackPosKaonCuts->SetPtExclusion(0.2,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.2,0.4);
  }

  if (suffix == "9") {
      TrackPosKaonCuts->SetPtRange(0.15,1.5);
      TrackNegKaonCuts->SetPtRange(0.15,1.5);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.5);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.5);
  }


  if (suffix == "10") {
      TrackPosKaonCuts->SetPtRange(0.15,1.3);
      TrackNegKaonCuts->SetPtRange(0.15,1.3);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
  }

  if (suffix == "11") {
      TrackPosKaonCuts->SetPtRange(0.15,1.3);
      TrackNegKaonCuts->SetPtRange(0.15,1.3);
      TrackPosKaonCuts->SetPtExclusion(0.2,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.2,0.4);
  }

  if (suffix == "12") {
      TrackPosKaonCuts->SetPtRange(0.15,1.3);
      TrackNegKaonCuts->SetPtRange(0.15,1.3);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.5);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.5);
  }

  if (suffix == "13") {
      TrackPosKaonCuts->SetPtRange(0.15,1.0);
      TrackNegKaonCuts->SetPtRange(0.15,1.0);
  }

  if (suffix == "14") {
      TrackPosKaonCuts->SetPtRange(0.15,1.2);
      TrackNegKaonCuts->SetPtRange(0.15,1.2);
  }

  if (suffix == "15") {
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, 3);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, 3);
  }

  if (suffix == "16") {
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, 3);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, 3);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
  }

  if (suffix == "17") {
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, 4);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, 4);
  }

  if (suffix == "18") {
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, 4);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, 4);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
  }

  if (suffix == "19") {
      TrackPosKaonCuts->SetFilterBit(128);
      TrackNegKaonCuts->SetFilterBit(128);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
  }

  if (suffix == "20") {
      TrackPosKaonCuts->SetFilterBit(128);
      TrackNegKaonCuts->SetFilterBit(128);
      TrackPosKaonCuts->SetPtRange(0.15,1.0);
      TrackNegKaonCuts->SetPtRange(0.15,1.0);
  }

  if (suffix == "21") {
      TrackPosKaonCuts->SetDCAVtxZ(2.0);
      TrackPosKaonCuts->SetDCAVtxXY(2.0);
      TrackNegKaonCuts->SetDCAVtxZ(2.0);
      TrackNegKaonCuts->SetDCAVtxXY(2.0);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
  }

  if (suffix == "22") {
      TrackPosKaonCuts->SetDCAVtxZ(1.0);
      TrackPosKaonCuts->SetDCAVtxXY(1.0);
      TrackNegKaonCuts->SetDCAVtxZ(1.0);
      TrackNegKaonCuts->SetDCAVtxXY(1.0);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
  }

  if (suffix == "23") {
      TrackPosKaonCuts->SetDCAVtxZ(1.0);
      TrackPosKaonCuts->SetDCAVtxXY(1.0);
      TrackNegKaonCuts->SetDCAVtxZ(1.0);
      TrackNegKaonCuts->SetDCAVtxXY(1.0);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
  }

  if (suffix == "24") {
      TrackPosKaonCuts->SetDCAVtxZ(0.5);
      TrackPosKaonCuts->SetDCAVtxXY(0.5);
      TrackNegKaonCuts->SetDCAVtxZ(0.5);
      TrackNegKaonCuts->SetDCAVtxXY(0.5);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
  }

  if (suffix == "25") {
      TrackPosKaonCuts->SetDCAVtxZ(0.5);
      TrackPosKaonCuts->SetDCAVtxXY(0.5);
      TrackNegKaonCuts->SetDCAVtxZ(0.5);
      TrackNegKaonCuts->SetDCAVtxXY(0.5);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
  }

  if (suffix == "26") {
      TrackPosKaonCuts->SetDCAVtxZ(0.8);
      TrackPosKaonCuts->SetDCAVtxXY(0.6);
      TrackNegKaonCuts->SetDCAVtxZ(0.8);
      TrackNegKaonCuts->SetDCAVtxXY(0.6);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
  }

  if (suffix == "27") {
      TrackPosKaonCuts->SetDCAVtxZ(0.8);
      TrackPosKaonCuts->SetDCAVtxXY(0.6);
      TrackNegKaonCuts->SetDCAVtxZ(0.8);
      TrackNegKaonCuts->SetDCAVtxXY(0.6);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
  }

  if (suffix == "28") {
      TrackPosKaonCuts->SetDCAVtxZ(0.3);
      TrackPosKaonCuts->SetDCAVtxXY(0.2);
      TrackNegKaonCuts->SetDCAVtxZ(0.3);
      TrackNegKaonCuts->SetDCAVtxXY(0.2);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
  }

  if (suffix == "29") {
      TrackPosKaonCuts->SetDCAVtxZ(0.3);
      TrackPosKaonCuts->SetDCAVtxXY(0.2);
      TrackNegKaonCuts->SetDCAVtxZ(0.3);
      TrackNegKaonCuts->SetDCAVtxXY(0.2);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
  }

  if (suffix == "30") {
      TrackPosKaonCuts->SetPtRange(0.15,1.0);
      TrackNegKaonCuts->SetPtRange(0.15,1.0);
      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, 3);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, 3);
  }

  if (suffix == "31") {
      TrackPosKaonCuts->SetPtRange(0.15,1.0);
      TrackNegKaonCuts->SetPtRange(0.15,1.0);
      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, 3);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, 3);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
  }

  if (suffix == "32") {
      TrackPosKaonCuts->SetPtRange(0.15,1.0);
      TrackNegKaonCuts->SetPtRange(0.15,1.0);
      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, 4);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, 4);
  }

  if (suffix == "33") {
      TrackPosKaonCuts->SetPtRange(0.15,1.0);
      TrackNegKaonCuts->SetPtRange(0.15,1.0);
      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, 4);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, 4);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
  }

  if (suffix == "34") {
      TrackPosKaonCuts->SetDCAVtxZ(1.0);
      TrackPosKaonCuts->SetDCAVtxXY(1.0);
      TrackNegKaonCuts->SetDCAVtxZ(1.0);
      TrackNegKaonCuts->SetDCAVtxXY(1.0);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetFilterBit(128);
      TrackNegKaonCuts->SetFilterBit(128);
  }

  if (suffix == "35") {
      TrackPosKaonCuts->SetDCAVtxZ(1.0);
      TrackPosKaonCuts->SetDCAVtxXY(1.0);
      TrackNegKaonCuts->SetDCAVtxZ(1.0);
      TrackNegKaonCuts->SetDCAVtxXY(1.0);
      TrackPosKaonCuts->SetPtRange(0.15,1.4);
      TrackNegKaonCuts->SetPtRange(0.15,1.4);
      TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
      TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
      TrackPosKaonCuts->SetFilterBit(128);
      TrackNegKaonCuts->SetFilterBit(128);
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
  double Phimass=TDatabasePDG::Instance()->GetParticle(333)->Mass();


  if (suffix != "0") {
    TrackCutsPhi->SetMinimalBooking(true);
  }


  // Now we define stuff we want for our Particle collection
  // Thanks, CINT - will not compile due to an illegal constructor
  // std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  // First we need to tell him about the particles we mix, from the
  // PDG code the mass is obtained.
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(333);

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

  // Number of bins
  std::vector<int> NBins;
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  std::vector<float> kMin;
  // minimum k* value
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  // maximum k* value
  std::vector<float> kMax;
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);

  AliFemtoDreamCollConfig *config =
      new AliFemtoDreamCollConfig("Femto", "Femto");
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetMultBinning(true);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(10);
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
      new AliAnalysisTaskNanoAODFemtoDreamPhi("AliAnalysisTaskNanoAODFemtoDreamPhi", isMC);
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
  coutputQA = mgr->CreateContainer(
      QAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}
