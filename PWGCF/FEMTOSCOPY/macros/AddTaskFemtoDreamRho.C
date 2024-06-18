AliAnalysisTaskSE *AddTaskFemtoDreamRho(bool isMC = false,
                                        bool doMcTruth = false,
                                        bool doAncestors = false,
                                        TString CentEst = "kHM",
                                        float fdPhidEta = 0.01,
                                        bool MCtemplatefit = false,
                                        float fSpherDown = 0.7,
                                        bool doCleaning = false,
                                        bool rejectKaon = false,
                                        bool doSpherSelect = false,
                                        bool doClosePairRejection = false,
                                        bool doProjections = false,
                                        float rhoPtThreshold = 0.,
                                        float rhoPtThresholdupper = 4.,
                                        float rhoCandInvMassLow = 0.075,
                                        float rhoCandInvMassHigh = 0.075,
                                        bool isSameCharge = false,
                                        bool isMCTrueRhoCombBkrg = false,
                                        bool isMCcheckedCombs = false,
                                        bool useNegativePairs = false,
                                        const char *cutVariation = "0")
{
  TString suffix = TString::Format("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler())
  {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);
  evtCuts->SetCutMinContrib(2);
  evtCuts->SetZVtxPosition(-10., 10.);

  printf("This task uses StandardCutsRun2!\n");

  evtCuts->SetSphericityCuts(0.0, 1);

  float fdPhi = 0;
  float fdEta = 0;

  if (doClosePairRejection)
  {
    fdPhi = fdPhidEta;
    fdEta = fdPhidEta;
  }
  else
  {
    fdPhi = 0;
    fdEta = 0;
  }

  if (doSpherSelect)
  {
    evtCuts->SetSphericityCuts(fSpherDown, 1.0, 0.5);
  }

  AliFemtoDreamTrackCuts *TrackPosProtonCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, true, true);
  TrackPosProtonCuts->SetCutCharge(1);
  // MC Template treatment
  if (!MCtemplatefit)
  {
    TrackPosProtonCuts->SetFilterBit(128); // Filterbit 5+6
    TrackPosProtonCuts->SetDCAVtxZ(0.2);
    TrackPosProtonCuts->SetDCAVtxXY(0.1);
  }
  else
  {
    TrackPosProtonCuts->SetFilterBit(128); // Filterbit 7 //Else the DCA is to tight
    TrackPosProtonCuts->SetDCAVtxZ(3.0);
    TrackPosProtonCuts->SetDCAVtxXY(3.0);
  }
  if (isMC && MCtemplatefit)
  {
    // TrackPosProtonCuts->SetPlotContrib(true);
    TrackPosProtonCuts->CheckParticleMothers(true);
    TrackPosProtonCuts->SetPlotDCADist(true);
    // TrackPosProtonCuts->SetOriginMultiplicityHists(true);
    TrackPosProtonCuts->SetFillQALater(false); // Be careful about this flag! When the MinimalBooking is set
  }

  AliFemtoDreamTrackCuts *TrackNegProtonCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, true, true);
  TrackNegProtonCuts->SetCutCharge(-1);
  // MC Template treatment
  if (!MCtemplatefit)
  {
    TrackNegProtonCuts->SetFilterBit(128); // Filterbit 5+6
    TrackNegProtonCuts->SetDCAVtxZ(0.2);
    TrackNegProtonCuts->SetDCAVtxXY(0.1);
  }
  else
  {
    TrackNegProtonCuts->SetFilterBit(128); // Filterbit 7 //Else the DCA is to tight
    TrackNegProtonCuts->SetDCAVtxZ(3.0);
    TrackNegProtonCuts->SetDCAVtxXY(3.0);
  }
  if (isMC && MCtemplatefit)
  {
    // TrackNegProtonCuts->SetPlotContrib(true);
    TrackNegProtonCuts->CheckParticleMothers(true);
    TrackNegProtonCuts->SetPlotDCADist(true);
    // TrackNegProtonCuts->SetOriginMultiplicityHists(true);
    TrackNegProtonCuts->SetFillQALater(false); // Be careful about this flag! When the MinimalBooking is set
  }

  if (suffix != "0")
  {
    TrackPosProtonCuts->SetMinimalBooking(false);
    TrackNegProtonCuts->SetMinimalBooking(false);
  }

  AliFemtoDreamTrackCuts *TrackPosPionCuts =
      AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, true, true);
  TrackPosPionCuts->SetCutCharge(1);
  // TrackPosPionCuts->SetPtRange(0.14, 10);
  //  MC Template treatment
  if (!MCtemplatefit)
  {
    TrackPosPionCuts->SetFilterBit(96); // Filterbit 5+6
    TrackPosPionCuts->SetDCAVtxZ(0.3);
    TrackPosPionCuts->SetDCAVtxXY(0.3);
  }
  else
  {
    TrackPosPionCuts->SetFilterBit(128); // Filterbit 7 //Else the DCA is to tight
    TrackPosPionCuts->SetDCAVtxZ(3.0);
    TrackPosPionCuts->SetDCAVtxXY(3.0);
  }
  if (isMC && MCtemplatefit)
  {
    // TrackPosPionCuts->SetPlotContrib(true);
    TrackPosPionCuts->CheckParticleMothers(true);
    TrackPosPionCuts->SetPlotDCADist(true);
    // TrackPosPionCuts->SetOriginMultiplicityHists(true);
    TrackPosPionCuts->SetFillQALater(false); // Be careful about this flag! When the MinimalBooking is set
  }

  AliFemtoDreamTrackCuts *TrackNegPionCuts =
      AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, true, true);
  TrackNegPionCuts->SetCutCharge(-1);
  // TrackNegPionCuts->SetPtRange(0.14, 10);
  //  MC Template treatment
  if (!MCtemplatefit)
  {
    TrackNegPionCuts->SetFilterBit(96); // Filterbit 5+6
    TrackNegPionCuts->SetDCAVtxZ(0.3);
    TrackNegPionCuts->SetDCAVtxXY(0.3);
  }
  else
  {
    TrackNegPionCuts->SetFilterBit(128); // Filterbit 7 //Else the DCA is to tight
    TrackNegPionCuts->SetDCAVtxZ(3.0);
    TrackNegPionCuts->SetDCAVtxXY(3.0);
  }
  // MC Template treatment
  if (isMC && MCtemplatefit)
  {
    // TrackNegPionCuts->SetPlotContrib(true);
    TrackNegPionCuts->CheckParticleMothers(true);
    TrackNegPionCuts->SetPlotDCADist(true);
    // TrackNegPionCuts->SetOriginMultiplicityHists(true);
    TrackNegPionCuts->SetFillQALater(false); // Be careful about this flag! When the MinimalBooking is set
  }

  // Only used if the doProjection flag is set
  AliFemtoDreamTrackCuts *TrackPosPionMinvCuts =
      AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, true, true);
  TrackPosPionMinvCuts->SetCutCharge(1);
  // MC Template treatment
  if (!MCtemplatefit)
  {
    TrackPosPionMinvCuts->SetFilterBit(96); // Filterbit 5+6
    TrackPosPionMinvCuts->SetDCAVtxZ(0.3);
    TrackPosPionMinvCuts->SetDCAVtxXY(0.3);
  }
  else
  {
    TrackPosPionMinvCuts->SetFilterBit(128); // Filterbit 7 //Else the DCA is to tight
    TrackPosPionMinvCuts->SetDCAVtxZ(3.0);
    TrackPosPionMinvCuts->SetDCAVtxXY(3.0);
  }
  if (isMC && MCtemplatefit)
  {
    // TrackPosPionCuts->SetPlotContrib(true);
    TrackPosPionMinvCuts->CheckParticleMothers(true);
    TrackPosPionMinvCuts->SetPlotDCADist(true);
    // TrackPosPionCuts->SetOriginMultiplicityHists(true);
    TrackPosPionMinvCuts->SetFillQALater(false); // Be careful about this flag! When the MinimalBooking is set
  }

  AliFemtoDreamTrackCuts *TrackNegPionMinvCuts =
      AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, true, true);
  TrackNegPionMinvCuts->SetCutCharge(-1);
  // MC Template treatment
  if (!MCtemplatefit)
  {
    TrackNegPionMinvCuts->SetFilterBit(96); // Filterbit 5+6
    TrackNegPionMinvCuts->SetDCAVtxZ(0.3);
    TrackNegPionMinvCuts->SetDCAVtxXY(0.3);
  }
  else
  {
    TrackNegPionMinvCuts->SetFilterBit(128); // Filterbit 7 //Else the DCA is to tight
    TrackNegPionMinvCuts->SetDCAVtxZ(3.0);
    TrackNegPionMinvCuts->SetDCAVtxXY(3.0);
  }
  // MC Template treatment
  if (isMC && MCtemplatefit)
  {
    // TrackNegPionCuts->SetPlotContrib(true);
    TrackNegPionMinvCuts->CheckParticleMothers(true);
    TrackNegPionMinvCuts->SetPlotDCADist(true);
    // TrackNegPionCuts->SetOriginMultiplicityHists(true);
    TrackNegPionMinvCuts->SetFillQALater(false); // Be careful about this flag! When the MinimalBooking is set
  }

  if (suffix != "0")
  {
    TrackPosPionCuts->SetMinimalBooking(false);
    TrackNegPionCuts->SetMinimalBooking(false);
  }

  AliFemtoDreamv0Cuts *TrackCutsRho = new AliFemtoDreamv0Cuts();
  TrackCutsRho->SetIsMonteCarlo(isMC);
  TrackCutsRho->SetAxisInvMassPlots(5000, 0, 5);
  TrackCutsRho->SetCutInvMass(0.150 / 2); // Should be determined from fit to the rho peak
  AliFemtoDreamTrackCuts *dummyCutsPos = new AliFemtoDreamTrackCuts();
  dummyCutsPos->SetIsMonteCarlo(isMC);
  AliFemtoDreamTrackCuts *dummyCutsNeg = new AliFemtoDreamTrackCuts();
  dummyCutsNeg->SetIsMonteCarlo(isMC);
  TrackCutsRho->SetPosDaugterTrackCuts(dummyCutsPos); // These are needed but will not be used
  TrackCutsRho->SetNegDaugterTrackCuts(dummyCutsNeg);
  TrackCutsRho->SetPDGCodePosDaug(211);
  TrackCutsRho->SetPDGCodeNegDaug(211);
  TrackCutsRho->SetPDGCodev0(113);
  TrackCutsRho->SetPtRange(rhoPtThreshold, rhoPtThresholdupper);
  TrackCutsRho->SetFineInvMassPtBins(true);

  if (suffix == "0")
  {
    TrackCutsRho->SetCutInvMass(0.150 / 2);
  }
  else if (suffix == "1")
  {
    TrackCutsRho->SetCutWindow(0.850, 0.900);
  }
  else if (suffix == "2")
  {
    TrackCutsRho->SetCutWindow(0.820, 0.900);
  }
  else if (suffix == "3")
  {
    TrackCutsRho->SetCutWindow(0.650, 0.730);
  }
  else if (suffix == "4")
  {
    TrackCutsRho->SetCutWindow(0.650, 0.700);
  }
  else if (suffix == "999")
  {
    TrackCutsRho->SetCutWindow(0., 5.0);
  }
  else
  {
    TrackCutsRho->SetCutWindow(rhoCandInvMassLow, rhoCandInvMassHigh);
  }

  // This needs further implementation
  // TrackCutsRho->SetEtaRange(-0.8, 0.8); //For now we also cut on the reconstructed rho, as to not be "too" forward (these are anyway probably not usable for femto due to the large Deta with the protons)
  //  Reject Kaons present in the sample:
  if (rejectKaon)
  {
    TrackCutsRho->SetKaonRejection(0.4, 0.6); // inv mass down and up
  }

  AliFemtoDreamTrackCuts *TrackCutsRhoMCTrue = new AliFemtoDreamTrackCuts();

  // now we create the task
  AliAnalysisTaskFemtoDreamRho *task =
      new AliAnalysisTaskFemtoDreamRho("AliAnalysisTaskFemtoDreamRho", isMC, doMcTruth, doCleaning, doAncestors, doProjections, rhoPtThreshold, isSameCharge, useNegativePairs, isMCTrueRhoCombBkrg, isMCcheckedCombs);
  // THIS IS VERY IMPORTANT ELSE YOU DONT PROCESS ANY EVENTS
  // kINT7 == Minimum bias
  // kHighMultV0 high multiplicity triggered by the V0 detector
  if (CentEst == "kInt7")
  {
    task->SetTrigger(AliVEvent::kINT7);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  }
  else if (CentEst == "kHM")
  {
    task->SetTrigger(AliVEvent::kHighMultV0);
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
  }
  else
  {
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
  if (doProjections)
  {
    PDGParticles.push_back(211); // 5 //pion+
    PDGParticles.push_back(211); // 6 //pion-
  }
  // if (doMCtrueRho)
  //{
  //   PDGParticles.push_back(113); // 5 //rhoMCTrue
  // }

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

  // The next part is for the result histograms. The order of hist. is the following:
  //                 Particle1     Particle2
  // Particle 1       Hist 1         Hist2
  //
  // Particle 2                      Hist3

  //          pi+       pi-       rho       p       Ap
  // pi+      pi+pi+1   pi+pi-2   pi+rho4   pi+p7   pi+Ap11
  //
  // pi-                pi+pi+3   pi-rho5   pi-p8   pi-Ap12
  //
  // rho                          rhorho6   rhop9   rhoAp13
  //
  // p                                      pp10    pAp14
  //
  // Ap                                             ApAp15

  //          pi+       pi-       rho       p       Ap
  // pi+      pi+pi+1   pi+pi-2   pi+rho3   pi+p4   pi+Ap5
  //
  // pi-                pi-pi-6   pi-rho7   pi-p8   pi-Ap9
  //
  // rho                          rhorho10  rhop11  rhoAp12
  //
  // p                                      pp13    pAp14
  //
  // Ap                                             ApAp15

  // Here with the DoProjection
  //           pi+       pi-       rho       p       Ap       pi+Minv           pi-Minv
  //  pi+      pi+pi+1   pi+pi-2   pi+rho3   pi+p4   pi+Ap5   pi+pi+Minv6       pi+pi-Minv7
  //
  //  pi-                pi-pi-8   pi-rho9   pi-p10  pi-Ap11  pi-pi+Minv12      pi-pi-Minv13
  //
  //  rho                          rhorho14  rhop15  rhoAp16  rhopi+Minv17      rhopi-Minv18
  //
  //  p                                      pp19    pAp20    ppi+Minv21        ppi-Minv22
  //
  //  Ap                                             ApAp23   Appi+Minv24       Appi-Minv25
  //
  //  pi+Minv                                                 pi+Minvpi+Minv26  pi+Minvpi-Minv27
  //
  //  pi-Minv                                                                   pi-Minvpi-Minv28

  // The same way the values for binning, minimum and maximum k* range have to be set!
  //  Number of bins
  std::vector<int> NBins;
  //  NBins.push_back(750);

  // minimum k* value
  std::vector<float> kMin;
  //  kMin.push_back(0.);

  // maximum k* value
  std::vector<float> kMax;

  std::vector<int> pairQA;
  //  pairQA.push_back(11);

  bool doMCtrueRho = false;
  float numberOfPart = (doProjections) ? 7. : 5.;

  for (int i = 0; i < (numberOfPart * (numberOfPart + 1) / 2); i++) // change from 5 to 6 is fIsMCTrue is set
  {
    std::cout << "check int i: " << i << std::endl;
    NBins.push_back(750);
    kMin.push_back(0.);
    kMax.push_back(3.);
    pairQA.push_back(0);
  }

  pairQA[0] = 11;  /// pi+ pi+     //pi+pi+
  pairQA[1] = 11;  /// pi+ pi-     //pi+pi-
  pairQA[2] = 12;  /// pi- pi-     //pi+rho //this needs the CPR check how the framwork handels these cases
  pairQA[3] = 11;  /// rho pi+     //pi+p
  pairQA[4] = 11;  /// rho pi-     //pi+Ap
  pairQA[5] = 11;  /// rho rho     //pi-pi-
  pairQA[6] = 12;  /// p   pi+     //pi-rho
  pairQA[7] = 11;  /// p   pi-     //pi-p
  pairQA[8] = 11;  /// p   rho     //pi-Ap
  pairQA[9] = 22;  /// p   p       //rhorho
  pairQA[10] = 21; /// p-  pi+     //rhop  //May not need CPR
  pairQA[11] = 21; /// p-  pi-     //rhoAp
  pairQA[12] = 11; /// p-  rho     //pp
  pairQA[13] = 11; /// p-  p       //pAp
  pairQA[14] = 11; /// p-  p-      //ApAp

  if (doProjections)
  {
    // override each entry
    pairQA[0] = 11; /// pi+ pi+     //pi+pi+
    pairQA[1] = 11; /// pi+ pi-     //pi+pi-
    pairQA[2] = 12; /// pi- pi-     //pi+rho //this needs the CPR check how the framwork handels these cases
    pairQA[3] = 11; /// rho pi+     //pi+p
    pairQA[4] = 11; /// rho pi-     //pi+Ap
    pairQA[5] = 11; /// rho rho     //pi+pi+Minv
    pairQA[6] = 11; /// p   pi+     //pi+pi-Minv

    pairQA[7] = 11;  /// rho rho     //pi-pi-
    pairQA[8] = 12;  /// p   pi+     //pi-rho
    pairQA[9] = 11;  /// p   pi-     //pi-p
    pairQA[10] = 11; /// p   rho     //pi-Ap
    pairQA[11] = 11; /// p-  pi-     //pi-pi+Minv
    pairQA[12] = 11; /// p-  rho     //pi-pi-Minv

    pairQA[13] = 22; /// p   p       //rhorho
    pairQA[14] = 21; /// p-  pi+     //rhop  //May not need CPR
    pairQA[15] = 21; /// p-  pi-     //rhoAp
    pairQA[16] = 21; /// p-  rho     //rhopi+Minv
    pairQA[17] = 21; /// p-  p       //rhopi-Minv

    pairQA[18] = 11; /// p-  rho     //pp
    pairQA[19] = 11; /// p-  p       //pAp
    pairQA[20] = 11; /// p-  rho     //ppi+Minv
    pairQA[21] = 11; /// p-  p       //ppi-Minv

    pairQA[22] = 11; /// p-  p-      //ApAp
    pairQA[23] = 11; /// p-  rho     //Appi+Minv
    pairQA[24] = 11; /// p-  p       //Appi-Minv

    pairQA[25] = 11; /// p-  rho     //pi+Minvpi+Minv
    pairQA[26] = 11; /// p-  p       //pi+Minvpi-Minv

    pairQA[27] = 11; /// p-  p       //pi-Minvpi-Minv
  }

  // for (int i = 0; i < (numberOfPart * (numberOfPart + 1) / 2); i++) // change from 5 to 6 is fIsMCTrue is set
  //{
  //   std::cout << "check2 int i: " << i << std::endl;
  //   pairQA[i] = 11; // phiphi
  //  }

  /*if (isMC) // override the values in case the MC True is used we don't need any of that (we may need this in order to perform the corss-check) //check if this breaks
  {
    pairQA[0] = 00;  // pi+pi+
    pairQA[1] = 00;  // pi+pi-
    pairQA[2] = 00;  // pi+rho
    pairQA[3] = 00;  // pi+p
    pairQA[4] = 00;  // pi+Ap
    pairQA[5] = 00;  // pi-pi-
    pairQA[6] = 00;  // pi-rho
    pairQA[7] = 00;  // pi-p
    pairQA[8] = 00;  // pi-Ap
    pairQA[9] = 00;  // rhorho
    pairQA[10] = 00; // rhop
    pairQA[11] = 00; // rhoAp
    pairQA[12] = 00; // pp
    pairQA[13] = 00; // pAp
    pairQA[14] = 00; // ApAp
  }*/

  std::vector<bool> closeRejection; // for now we don't use it, verify with MC if that is needed

  // pair rejection
  closeRejection.push_back(false); // pi+ pi+
  closeRejection.push_back(false); // pi+ pi-
  closeRejection.push_back(false); // pi+rho
  closeRejection.push_back(false); // pi+p
  closeRejection.push_back(false); // pi+Ap
  closeRejection.push_back(false); // pi-pi-
  closeRejection.push_back(false); // pi-rho
  closeRejection.push_back(false); // pi-p
  closeRejection.push_back(false); // pi-Ap
  closeRejection.push_back(false); // rhorho
  closeRejection.push_back(false); // rhop
  closeRejection.push_back(false); // rhoAp
  closeRejection.push_back(false); // pp
  closeRejection.push_back(false); // pAp
  closeRejection.push_back(false); // ApAp

  if (doProjections)
  {
    closeRejection.push_back(false); // pi+ pi+
    closeRejection.push_back(false); // pi+ pi-
    closeRejection.push_back(false); // pi+rho
    closeRejection.push_back(false); // pi+p
    closeRejection.push_back(false); // pi+Ap
    closeRejection.push_back(false); // pi-pi-
    closeRejection.push_back(false); // pi-rho
    closeRejection.push_back(false); // pi-p
    closeRejection.push_back(false); // pi-Ap
    closeRejection.push_back(false); // rhorho
    closeRejection.push_back(false); // rhop
    closeRejection.push_back(false); // rhoAp
    closeRejection.push_back(false); // pp
  }

  // for (int i = 0; i < (numberOfPart * (numberOfPart + 1) / 2); i++) // change from 5 to 6 is fIsMCTrue is set
  //{
  //   if (i > 14)
  //   {
  //     std::cout << "check2 int i: " << i << std::endl;
  //    closeRejection.push_back(true);
  //   }
  //}

  if (suffix == "99")
  {
    // Deactivate the ClosePairRejection
    fdPhidEta = 0.;
    closeRejection.clear();
    closeRejection.push_back(false); // pi+ pi+
    closeRejection.push_back(false); // pi+ pi-
    closeRejection.push_back(false); // pi- pi-
    closeRejection.push_back(false); // rho pi+
    closeRejection.push_back(false); // rho pi-
    closeRejection.push_back(false); // rho rho
    closeRejection.push_back(false); // p   pi+
    closeRejection.push_back(false); // p   pi-
    closeRejection.push_back(false); // p   rho
    closeRejection.push_back(false); // p   p
    closeRejection.push_back(false); // p-  pi+
    closeRejection.push_back(false); // p-  pi-
    closeRejection.push_back(false); // p-  rho
    closeRejection.push_back(false); // p-  p
    closeRejection.push_back(false); // p-  p-
    if (false)
    {
      closeRejection.push_back(false); // pi+ rhoMC
      closeRejection.push_back(false); // pi+ rhoMC
      closeRejection.push_back(false); // prho
      closeRejection.push_back(false); // apap
      closeRejection.push_back(false); // apphi
      closeRejection.push_back(false); // phiphi
    }
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
  config->SetDeltaEtaMax(fdEta);
  config->SetDeltaPhiMax(fdPhi);
  config->SetMixingDepth(10); // AN
  config->SetkTBinning(true);
  config->SetmTBinning(true);
  config->SetExtendedQAPairs(pairQA);
  config->SetUseEventMixing(true);
  config->SetMinimalBookingME(false);
  config->SetdPhidEtaPlots(true);
  config->SetdPhidEtaPlotsSmallK(true);
  config->SetMinvKtandRelativeKBinning(true); // Check this in case of MC

  std::cout << " Check the addTask config" << std::endl;
  std::cout << " ZVtxBins.size() " << ZVtxBins.size() << std::endl;
  std::cout << " MultBins.size() " << MultBins.size() << std::endl;
  std::cout << " PDGParticles.size() " << PDGParticles.size() << std::endl;
  std::cout << " NBins.size() " << NBins.size() << std::endl;
  std::cout << " kMin.size() " << kMin.size() << std::endl;
  std::cout << " kMax.size() " << kMax.size() << std::endl;
  std::cout << " closeRejection.size() " << closeRejection.size() << std::endl;
  std::cout << " pairQA.size() " << pairQA.size() << std::endl;
  std::cout << " config->SetDeltaEtaMax(fdEta); " << fdEta << std::endl;
  std::cout << " config->SetDeltaPhiMax(fdPhi); " << fdPhi << std::endl;

  if (isMC)
  {
    config->SetMomentumResolution(true);
    config->SetPhiEtaBinnign(true);
    if (doAncestors)
    {
      config->SetAncestors(doAncestors);
    }
  }
  else
  {
    std::cout << "You are trying to request the Momentum Resolution without MC "
                 "Info; fix it wont work! \n";
  }

  if (doMcTruth) // turn the troublesome histos i.e. detadphi* off for now, can be fixed by setting the proper values at each TPC rad manually.
  {
    config->SetPhiEtaBinnign(false);
    config->SetdPhidEtaPlots(false);
    config->SetdPhidEtaPlotsSmallK(false);
    // config->SetMinimalBookingME(true);
  }

  // Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackPosProtonCuts);
  task->SetAntiProtonCuts(TrackNegProtonCuts);
  task->SetRhoCuts(TrackCutsRho);
  // if (TrackCutsRhoMCTrue)
  // {
  //   task->SetRhoMCTrueCuts(TrackCutsRhoMCTrue);
  //  }
  task->SetPosPionCuts(TrackPosPionCuts);
  task->SetNegPionCuts(TrackNegPionCuts);
  if (doProjections)
  {
    task->SetPosPionMinvCuts(TrackPosPionMinvCuts);
    task->SetNegPionMinvCuts(TrackNegPionMinvCuts);
  }
  task->SetCollectionConfig(config);
  task->SetDoCleaning(doCleaning);
  task->SetIsMC(isMC);
  task->SetDoMcTruth(doMcTruth);
  task->SetDoAncestors(doAncestors);

  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);

  AliAnalysisDataContainer *coutputQA;
  TString addon = "";
  if (CentEst == "kInt7")
  {
    addon += "MB";
  }
  else if (CentEst == "kHM")
  {
    addon += "HM";
  }
  TString QAName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputQA = mgr->CreateContainer(QAName.Data(), TList::Class(),
                                   AliAnalysisManager::kOutputContainer,
                                   Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}
