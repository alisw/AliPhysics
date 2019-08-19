AliAnalysisTaskSE *AddTaskFemtoDreamPhi(bool isMC = false,
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


  AliFemtoDreamTrackCuts *TrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetCutCharge(-1);

//  if (suffix != "0") {
//    TrackCuts->SetMinimalBooking(true);
//    AntiTrackCuts->SetMinimalBooking(true);
//  }

  AliFemtoDreamTrackCuts *TrackPosKaonCuts = AliFemtoDreamTrackCuts::PrimKaonCuts(isMC);
  TrackPosKaonCuts->SetCutCharge(1);
//  TrackPosKaonCuts->SetPID(AliPID::kKaon, 999, 5);
//  TrackPosKaonCuts->SetPlotTOFMass(true);

  AliFemtoDreamTrackCuts *TrackNegKaonCuts = AliFemtoDreamTrackCuts::PrimKaonCuts(isMC);
  TrackNegKaonCuts->SetCutCharge(-1);
//  TrackNegKaonCuts->SetPID(AliPID::kKaon, 999, 5);
//  TrackNegKaonCuts->SetPlotTOFMass(true);

//  if (suffix != "0") {
//    TrackPosKaonCuts->SetMinimalBooking(true);
//    TrackNegKaonCuts->SetMinimalBooking(true);
//  }

  AliFemtoDreamv0Cuts *TrackCutsPhi = new AliFemtoDreamv0Cuts();
  TrackCutsPhi->SetIsMonteCarlo(isMC);
  TrackCutsPhi->SetAxisInvMassPlots(400, 0.9, 1.2);
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

//  //cutwindow
//  if (suffix == "1") {
//    TrackCutsPhi->SetCutWindow(1.06,1.08);
//  }
//  //invmasscuts
//  if (suffix == "2") {
//    TrackCutsPhi->SetCutInvMass(0.006);
//  }
//  if (suffix == "3") {
//    TrackCutsPhi->SetCutInvMass(0.004);
//  }
//  //sphericitycuts
//  if (suffix == "4") {
//    evtCuts->SetSphericityCuts(0.6,1);
//  }
//  if (suffix == "5") {
//    evtCuts->SetSphericityCuts(0.65,1);
//  }
//  if (suffix == "6") {
//    evtCuts->SetSphericityCuts(0.7,1);
//  }
//  if (suffix == "7") {
//    evtCuts->SetSphericityCuts(0.75,1);
//  }
//  if (suffix == "8") {
//    evtCuts->SetSphericityCuts(0.8,1);
//  }
//  if (suffix == "9") {
//    evtCuts->SetSphericityCuts(0.85,1);
//  }
//  if (suffix == "10") {
//    evtCuts->SetSphericityCuts(0.9,1);
//  }
//  if (suffix == "11") {
//    evtCuts->SetSphericityCuts(0.95,1);
//  }
//  //sphericity and cutwindow
//  if (suffix == "12") {
//    TrackCutsPhi->SetCutWindow(1.08,1.13);
//    evtCuts->SetSphericityCuts(0.6,1);
//  }
//  if (suffix == "13") {
//    TrackCutsPhi->SetCutWindow(1.13,1.18);
//    evtCuts->SetSphericityCuts(0.6,1);
//  }
//  if (suffix == "14") {
//    TrackCutsPhi->SetCutWindow(1.08,1.18);
//    evtCuts->SetSphericityCuts(0.6,1);
//  }
//  if (suffix == "15") {
//    TrackCutsPhi->SetCutWindow(1.08,1.13);
//    evtCuts->SetSphericityCuts(0.7,1);
//  }
//  if (suffix == "16") {
//    TrackCutsPhi->SetCutWindow(1.13,1.18);
//    evtCuts->SetSphericityCuts(0.7,1);
//  }
//  if (suffix == "17") {
//    TrackCutsPhi->SetCutWindow(1.08,1.18);
//    evtCuts->SetSphericityCuts(0.7,1);
//  }
//  if (suffix == "18") {
//    TrackCutsPhi->SetCutWindow(1.08,1.13);
//    evtCuts->SetSphericityCuts(0.8,1);
//  }
//  if (suffix == "19") {
//    TrackCutsPhi->SetCutWindow(1.13,1.18);
//    evtCuts->SetSphericityCuts(0.8,1);
//  }
//  if (suffix == "20") {
//    TrackCutsPhi->SetCutWindow(1.08,1.18);
//    evtCuts->SetSphericityCuts(0.8,1);
//  }
//  if (suffix == "21") {
//    TrackCutsPhi->SetCutWindow(1.08,1.13);
//    evtCuts->SetSphericityCuts(0.9,1);
//  }
//  if (suffix == "22") {
//    TrackCutsPhi->SetCutWindow(1.13,1.18);
//    evtCuts->SetSphericityCuts(0.9,1);
//  }
//  if (suffix == "23") {
//    TrackCutsPhi->SetCutWindow(1.08,1.18);
//    evtCuts->SetSphericityCuts(0.9,1);
//  }

//  if (suffix != "0") {
//    TrackCutsPhi->SetMinimalBooking(true);
//  }


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

  std::vector<int> pairQA;
  pairQA.push_back(11);
  pairQA.push_back(0);
  pairQA.push_back(12);
  pairQA.push_back(11);
  pairQA.push_back(12);
  pairQA.push_back(0);

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
  config->SetExtendedQAPairs(pairQA);
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

    */
    if (isMC) {
      config->SetMomentumResolution(true);
    } else {
      std::cout << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
    }

    /*
    if (isMC) {
      config->SetPhiEtaBinnign(true);
    } else {
      std::cout << "You are trying to request the Eta Phi Plots without MC Info;
  fix it wont work! \n";
    }
  }
   */
  // now we create the task
  AliAnalysisTaskFemtoDreamPhi *task =
      new AliAnalysisTaskFemtoDreamPhi("AliAnalysisTaskFemtoDreamPhi", isMC);
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
  TString QAName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputQA = mgr->CreateContainer(
      QAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}
