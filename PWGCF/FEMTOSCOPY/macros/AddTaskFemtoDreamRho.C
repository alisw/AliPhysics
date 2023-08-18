AliAnalysisTaskSE* AddTaskFemtoDreamRho(bool isMC = false,
                                        TString CentEst = "kHM",
                                        float fdPhidEta = 0.01,
                                        bool MCtemplatefit = false,
                                        float fSpherDown = 0.7,
                                        bool doCleaning = false,
                                        bool rejectKaon = false,
                                        const char* cutVariation = "0") {
          
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
  evtCuts->SetSphericityCuts(0.0, 1);

  if (suffix == "1") {
    evtCuts->SetSphericityCuts(fSpherDown, 1.0, 0.5);
  }
                                        
  AliFemtoDreamTrackCuts *TrackPosProtonCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, true, true);
  TrackPosProtonCuts->SetCutCharge(1);
  //MC Template treatment
  if ( !MCtemplatefit ) {
    TrackPosProtonCuts->SetFilterBit(128); // Filterbit 5+6
    TrackPosProtonCuts->SetDCAVtxZ(0.2);
    TrackPosProtonCuts->SetDCAVtxXY(0.1);
  } else {
    TrackPosProtonCuts->SetFilterBit(128); // Filterbit 7 //Else the DCA is to tight
    TrackPosProtonCuts->SetDCAVtxZ(3.0);
    TrackPosProtonCuts->SetDCAVtxXY(3.0);
  }
  if ( isMC && MCtemplatefit ) {
    //TrackPosProtonCuts->SetPlotContrib(true);
    TrackPosProtonCuts->CheckParticleMothers(true);
    TrackPosProtonCuts->SetPlotDCADist(true);
    //TrackPosProtonCuts->SetOriginMultiplicityHists(true);
    TrackPosProtonCuts->SetFillQALater(false); //Be careful about this flag! When the MinimalBooking is set
  }
                                        
  AliFemtoDreamTrackCuts *TrackNegProtonCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, true, true);
  TrackNegProtonCuts->SetCutCharge(-1);
  //MC Template treatment
  if ( !MCtemplatefit ) {
    TrackNegProtonCuts->SetFilterBit(128); // Filterbit 5+6
    TrackNegProtonCuts->SetDCAVtxZ(0.2);
    TrackNegProtonCuts->SetDCAVtxXY(0.1);
  } else {
    TrackNegProtonCuts->SetFilterBit(128); // Filterbit 7 //Else the DCA is to tight
    TrackNegProtonCuts->SetDCAVtxZ(3.0);
    TrackNegProtonCuts->SetDCAVtxXY(3.0);
  }
  if ( isMC && MCtemplatefit ) {
    //TrackNegProtonCuts->SetPlotContrib(true);
    TrackNegProtonCuts->CheckParticleMothers(true);
    TrackNegProtonCuts->SetPlotDCADist(true);
    //TrackNegProtonCuts->SetOriginMultiplicityHists(true);
    TrackNegProtonCuts->SetFillQALater(false); //Be careful about this flag! When the MinimalBooking is set
  }

  if (suffix != "0") {
      TrackPosProtonCuts->SetMinimalBooking(false);
      TrackNegProtonCuts->SetMinimalBooking(false);
  }

  AliFemtoDreamTrackCuts *TrackPosPionCuts =
      AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, true, true);
  TrackPosPionCuts->SetCutCharge(1);
  //MC Template treatment
  if ( !MCtemplatefit ) {
    TrackPosPionCuts->SetFilterBit(96); // Filterbit 5+6
    TrackPosPionCuts->SetDCAVtxZ(0.3);
    TrackPosPionCuts->SetDCAVtxXY(0.3);
  } else {
    TrackPosPionCuts->SetFilterBit(128); // Filterbit 7 //Else the DCA is to tight
    TrackPosPionCuts->SetDCAVtxZ(3.0);
    TrackPosPionCuts->SetDCAVtxXY(3.0);
  }
  if ( isMC && MCtemplatefit ) {
    //TrackPosPionCuts->SetPlotContrib(true);
    TrackPosPionCuts->CheckParticleMothers(true);
    TrackPosPionCuts->SetPlotDCADist(true);
    //TrackPosPionCuts->SetOriginMultiplicityHists(true);
    TrackPosPionCuts->SetFillQALater(false); //Be careful about this flag! When the MinimalBooking is set
  }
                                        
  AliFemtoDreamTrackCuts *TrackNegPionCuts =
      AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, true, true);
  TrackNegPionCuts->SetCutCharge(-1);
  //MC Template treatment
  if ( !MCtemplatefit ) {
    TrackNegPionCuts->SetFilterBit(96); // Filterbit 5+6
    TrackNegPionCuts->SetDCAVtxZ(0.3);
    TrackNegPionCuts->SetDCAVtxXY(0.3);
  } else {
    TrackNegPionCuts->SetFilterBit(128); // Filterbit 7 //Else the DCA is to tight
    TrackNegPionCuts->SetDCAVtxZ(3.0);
    TrackNegPionCuts->SetDCAVtxXY(3.0);
  }
  //MC Template treatment
  if ( isMC && MCtemplatefit ) {
    //TrackNegPionCuts->SetPlotContrib(true);
    TrackNegPionCuts->CheckParticleMothers(true);
    TrackNegPionCuts->SetPlotDCADist(true);
    //TrackNegPionCuts->SetOriginMultiplicityHists(true);
    TrackNegPionCuts->SetFillQALater(false); //Be careful about this flag! When the MinimalBooking is set
  }

  if (suffix != "0") {
    TrackPosPionCuts->SetMinimalBooking(false);
    TrackNegPionCuts->SetMinimalBooking(false);
  }
                                        
  AliFemtoDreamv0Cuts *TrackCutsRho = new AliFemtoDreamv0Cuts();
  TrackCutsRho->SetIsMonteCarlo(isMC);
  TrackCutsRho->SetAxisInvMassPlots(5000, 0, 5);
  TrackCutsRho->SetCutInvMass(0.150); //Should be determined from fit to the rho peak
  AliFemtoDreamTrackCuts *dummyCutsPos = new AliFemtoDreamTrackCuts();
  dummyCutsPos->SetIsMonteCarlo(isMC);
  AliFemtoDreamTrackCuts *dummyCutsNeg = new AliFemtoDreamTrackCuts();
  dummyCutsNeg->SetIsMonteCarlo(isMC);
  TrackCutsRho->SetPosDaugterTrackCuts(dummyCutsPos); //These are needed but will not be used
  TrackCutsRho->SetNegDaugterTrackCuts(dummyCutsNeg);
  TrackCutsRho->SetPDGCodePosDaug(211);
  TrackCutsRho->SetPDGCodeNegDaug(211);
  TrackCutsRho->SetPDGCodev0(113);
  //Reject Kaons present in the sample:
  if (rejectKaon)TrackCutsRho->SetKaonRejection(0.4,0.6); //inv mass down and up
  // now we create the task
  AliAnalysisTaskFemtoDreamRho *task =
      new AliAnalysisTaskFemtoDreamRho("AliAnalysisTaskFemtoDreamRho", isMC, doCleaning);
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
  
  // Now we define stuff we want for our Particle collection
  // Thanks, CINT - will not compile due to an illegal constructor
  // std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  // First we need to tell him about the particles we mix, from the
  // PDG code the mass is obtained.
  std::vector<int> PDGParticles;
  PDGParticles.push_back(211);  // 0 //pion+
  PDGParticles.push_back(211);  // 1 //pion-
  PDGParticles.push_back(113);  // 2 //rho
  PDGParticles.push_back(2212); // 3 //proton+
  PDGParticles.push_back(2212); // 4 //proton-
                                        
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

  //The next part is for the result histograms. The order of hist. is the following:
  //                Particle1     Particle2
  //Particle 1       Hist 1         Hist2
  //
  //Particle 2                      Hist3
  //The same way the values for binning, minimum and maximum k* range have to be set!
  // Number of bins
  std::vector<int> NBins;
  //  NBins.push_back(750);

  // minimum k* value
  std::vector<float> kMin;
  //  kMin.push_back(0.);

  // maximum k* value
  std::vector<float> kMax;
  //  kMax.push_back(3.);

  std::vector<int> pairQA;
  //  pairQA.push_back(11);

  for (int i = 0; i < ( 5*(5+1)/2 ); i++) {
    NBins.push_back(750);
    kMin.push_back(0.);
    kMax.push_back(3.);
    pairQA.push_back(0);
  }
  
  pairQA[0] = 11;   // pp
  pairQA[1] = 11;   // pap
  pairQA[2] = 11;   // prho
  pairQA[3] = 11;   // apap
  pairQA[4] = 11;   // apphi
  pairQA[5] = 11;   // phiphi
  pairQA[6] = 11;   // pp
  pairQA[7] = 11;   // pap
  pairQA[8] = 11;   // prho
  pairQA[9] = 11;   // apap
  pairQA[10] = 11;  // apphi
  pairQA[11] = 11;  // phiphi
  pairQA[12] = 11;  // pp
  pairQA[13] = 11;  // pap
  pairQA[14] = 11;  // prho

  if (isMC) {
    pairQA[0] = 0;  // pp
    pairQA[1] = 0;  // pap
    pairQA[2] = 0;  // prho
    pairQA[3] = 0;  // apap
    pairQA[4] = 0;  // apphi
    pairQA[5] = 0;  // phiphi
    pairQA[6] = 0;  // pp
    pairQA[7] = 0;  // pap
    pairQA[8] = 0;  // prho
    pairQA[9] = 0;  // apap
    pairQA[10] = 0; // apphi
    pairQA[11] = 0; // phiphi
    pairQA[12] = 0; // pp
    pairQA[13] = 0; // pap
    pairQA[14] = 0; // prho
  }

  //pair rejection
  std::vector<bool> closeRejection;
  closeRejection.push_back(true);   // pi+ pi+
  closeRejection.push_back(false);  // pi+ pi- 
  closeRejection.push_back(true);   // pi- pi-
  closeRejection.push_back(false);  // rho pi+
  closeRejection.push_back(false);  // rho pi-
  closeRejection.push_back(false);  // rho rho
  closeRejection.push_back(true);   // p   pi+
  closeRejection.push_back(false);  // p   pi-
  closeRejection.push_back(false);  // p   rho
  closeRejection.push_back(true);   // p   p
  closeRejection.push_back(false);  // p-  pi+
  closeRejection.push_back(true);   // p-  pi-
  closeRejection.push_back(false);  // p-  rho
  closeRejection.push_back(false);  // p-  p
  closeRejection.push_back(true);   // p-  p-
  
  if (suffix == "99") {
    //Deactivate the ClosePairRejection
    fdPhidEta=0.;
    closeRejection.clear();
    closeRejection.push_back(false);  // pi+ pi+
    closeRejection.push_back(false);  // pi+ pi- 
    closeRejection.push_back(false);  // pi- pi-
    closeRejection.push_back(false);  // rho pi+
    closeRejection.push_back(false);  // rho pi-
    closeRejection.push_back(false);  // rho rho
    closeRejection.push_back(false);  // p   pi+
    closeRejection.push_back(false);  // p   pi-
    closeRejection.push_back(false);  // p   rho
    closeRejection.push_back(false);  // p   p
    closeRejection.push_back(false);  // p-  pi+
    closeRejection.push_back(false);  // p-  pi-
    closeRejection.push_back(false);  // p-  rho
    closeRejection.push_back(false);  // p-  p
    closeRejection.push_back(false);  // p-  p-
  }


  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto");
  config->SetPtQA(true);
  config->SetMassQA(true);
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetMultBinning(true);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetDeltaEtaMax(fdPhidEta); // https://alice-notes.web.cern.ch/system/files/notes/analysis/616/2018-08-10-NotepPb.pdf
  config->SetDeltaPhiMax(fdPhidEta);
  config->SetMixingDepth(10); // AN
  config->SetkTBinning(true);
  config->SetmTBinning(true);
  config->SetExtendedQAPairs(pairQA);
  config->SetUseEventMixing(true); // Check this flag
  config->SetMinimalBookingME(false);
  config->SetdPhidEtaPlots(true);
  config->SetdPhidEtaPlotsSmallK(true);
                                        
  if (isMC) {
    config->SetMomentumResolution(true);
    config->SetPhiEtaBinnign(true);
  } else {
    std::cout << "You are trying to request the Momentum Resolution without MC "
                 "Info; fix it wont work! \n";
  }

  // Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackPosProtonCuts);
  task->SetAntiProtonCuts(TrackNegProtonCuts);
  task->SetRhoCuts(TrackCutsRho);
  task->SetPosPionCuts(TrackPosPionCuts);
  task->SetNegPionCuts(TrackNegPionCuts);
  task->SetCollectionConfig(config);
  task->SetDoCleaning(doCleaning);
  task->SetIsMC(isMC);
  
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
  coutputQA = mgr->CreateContainer(QAName.Data(), TList::Class(),
                                   AliAnalysisManager::kOutputContainer,
                                   Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}
