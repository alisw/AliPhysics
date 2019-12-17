#include "TROOT.h"
#include "TSystem.h"
AliAnalysisTaskSE* AddTaskFemtoDreamDeuteron(
  bool isMC = false,//1
  TString CentEst = "kInt7",//2
  bool DCAPlots = false,//3
  bool CombSigma = false,//4
  bool ContributionSplitting = false//5
) {

  //Framework specific blabla
  // the manager is static, so get the existing manager via the static method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  // just to see if all went well, check if the input event handler has been
  // connected
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  //Now we need to setup the Event cuts, we use the Don't worry event cuts
  //from the ALICE DPG which is the offical thing to use.
  AliFemtoDreamEventCuts *evtCuts =
    AliFemtoDreamEventCuts::StandardCutsRun2();
  //This sets the method we want to use to clean up events with negative or too
  //low multiplicity. Usually you use the matching multiplicity estiamtor in your
  //event collection
  evtCuts->CleanUpMult(false, false, false, true);

  //Track Cuts are defined here for deuterons
  AliFemtoDreamTrackCuts *TrackCutsDeuteronDCA = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsDeuteronDCA->SetCutCharge(1);
  AliFemtoDreamTrackCuts *TrackCutsDeuteronMass = new AliFemtoDreamTrackCuts();
  TrackCutsDeuteronMass->SetPlotDCADist(DCAPlots);
  TrackCutsDeuteronMass->SetPlotCombSigma(CombSigma);
  TrackCutsDeuteronMass->SetPlotContrib(ContributionSplitting);
  TrackCutsDeuteronMass->SetIsMonteCarlo(isMC);
  TrackCutsDeuteronMass->SetCutCharge(1);
  TrackCutsDeuteronMass->SetFilterBit(128);
  //You want to optimize this to only select regions with reasonable high purity
  TrackCutsDeuteronMass->SetPtRange(0.4, 4.0);
  //You want to ensure good tracking quality
  TrackCutsDeuteronMass->SetEtaRange(-0.8, 0.8);
  TrackCutsDeuteronMass->SetNClsTPC(80);
  TrackCutsDeuteronMass->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
  //You want to select primaries!
  TrackCutsDeuteronMass->SetDCAVtxZ(0.2);
  TrackCutsDeuteronMass->SetDCAVtxXY(0.1);
  TrackCutsDeuteronMass->SetCutSharedCls(true);
  TrackCutsDeuteronMass->SetCutTPCCrossedRows(true, 70, 0.83);
  //Here you define the PID
  TrackCutsDeuteronMass->SetPID(AliPID::kDeuteron, 999.);
  //We are looking for pions rejecting them would be obstructive.
  TrackCutsDeuteronMass->SetRejLowPtPionsTOF(false);
  //this checks if the sigma of the wanted hypothesis is the smallest, and if
  //another particle has a smaller sigma, the track is rejected.
  TrackCutsDeuteronMass->SetCutSmallestSig(false);

  //The same things for anti deuterons
  AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronDCA = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsAntiDeuteronDCA->SetCutCharge(-1);
  AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronMass = new AliFemtoDreamTrackCuts();
  TrackCutsAntiDeuteronMass->SetPlotDCADist(DCAPlots);
  TrackCutsAntiDeuteronMass->SetPlotCombSigma(CombSigma);
  TrackCutsAntiDeuteronMass->SetPlotContrib(ContributionSplitting);
  TrackCutsAntiDeuteronMass->SetIsMonteCarlo(isMC);
  TrackCutsAntiDeuteronMass->SetCutCharge(-1);
  TrackCutsAntiDeuteronMass->SetFilterBit(128);
  //You want to optimize this to only select regions with reasonable high purity
  TrackCutsAntiDeuteronMass->SetPtRange(0.4, 4.0);
  //You want to ensure good tracking quality
  TrackCutsAntiDeuteronMass->SetEtaRange(-0.8, 0.8);
  TrackCutsAntiDeuteronMass->SetNClsTPC(80);
  TrackCutsAntiDeuteronMass->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
  //You want to select primaries!
  TrackCutsAntiDeuteronMass->SetDCAVtxZ(0.2);
  TrackCutsAntiDeuteronMass->SetDCAVtxXY(0.1);
  TrackCutsAntiDeuteronMass->SetCutSharedCls(true);
  TrackCutsAntiDeuteronMass->SetCutTPCCrossedRows(true, 70, 0.83);
  //Here you define the PID
  TrackCutsAntiDeuteronMass->SetPID(AliPID::kDeuteron, 999.);
  //We are looking for pions rejecting them would be obstructive.
  TrackCutsAntiDeuteronMass->SetRejLowPtPionsTOF(false);
  //this checks if the sigma of the wanted hypothesis is the smallest, and if
  //another particle has a smaller sigma, the track is rejected.
  TrackCutsAntiDeuteronMass->SetCutSmallestSig(false);

  //Track Cuts are defined here for deuterons
  AliFemtoDreamTrackCuts *TrackCutsProtonDCA = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsProtonDCA->SetCutCharge(1);
  AliFemtoDreamTrackCuts *TrackCutsProtonMass = new AliFemtoDreamTrackCuts();
  TrackCutsProtonMass->SetPlotDCADist(DCAPlots);
  TrackCutsProtonMass->SetPlotCombSigma(CombSigma);
  TrackCutsProtonMass->SetPlotContrib(ContributionSplitting);
  TrackCutsProtonMass->SetIsMonteCarlo(isMC);
  TrackCutsProtonMass->SetCutCharge(1);
  TrackCutsProtonMass->SetFilterBit(128);
  //You want to optimize this to only select regions with reasonable high purity
  TrackCutsProtonMass->SetPtRange(0.5, 4.05);
  //You want to ensure good tracking quality
  TrackCutsProtonMass->SetEtaRange(-0.8, 0.8);
  TrackCutsProtonMass->SetNClsTPC(80);
  TrackCutsProtonMass->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
  //You want to select primaries!
  TrackCutsProtonMass->SetDCAVtxZ(0.2);
  TrackCutsProtonMass->SetDCAVtxXY(0.1);
  TrackCutsProtonMass->SetCutSharedCls(true);
  TrackCutsProtonMass->SetCutTPCCrossedRows(true, 70, 0.83);
  //Here you define the PID
  TrackCutsProtonMass->SetPID(AliPID::kProton, 999.);
  TrackCutsProtonMass->SetRejLowPtPionsTOF(false);
  TrackCutsProtonMass->SetCutSmallestSig(false);

  //The same things for anti deuterons
  AliFemtoDreamTrackCuts *TrackCutsAntiProtonDCA = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsAntiProtonDCA->SetCutCharge(-1);
  AliFemtoDreamTrackCuts *TrackCutsAntiProtonMass = new AliFemtoDreamTrackCuts();
  TrackCutsAntiProtonMass->SetPlotDCADist(DCAPlots);
  TrackCutsAntiProtonMass->SetPlotCombSigma(CombSigma);
  TrackCutsAntiProtonMass->SetPlotContrib(ContributionSplitting);
  TrackCutsAntiProtonMass->SetIsMonteCarlo(isMC);
  TrackCutsAntiProtonMass->SetCutCharge(-1);
  TrackCutsAntiProtonMass->SetFilterBit(128);
  //You want to optimize this to only select regions with reasonable high purity
  TrackCutsAntiProtonMass->SetPtRange(0.5, 4.05);
  //You want to ensure good tracking quality
  TrackCutsAntiProtonMass->SetEtaRange(-0.8, 0.8);
  TrackCutsAntiProtonMass->SetNClsTPC(80);
  TrackCutsAntiProtonMass->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
  //You want to select primaries!
  TrackCutsAntiProtonMass->SetDCAVtxZ(0.2);
  TrackCutsAntiProtonMass->SetDCAVtxXY(0.1);
  TrackCutsAntiProtonMass->SetCutSharedCls(true);
  TrackCutsAntiProtonMass->SetCutTPCCrossedRows(true, 70, 0.83);
  //Here you define the PID
  TrackCutsAntiProtonMass->SetPID(AliPID::kProton, 999.);
  TrackCutsAntiProtonMass->SetRejLowPtPionsTOF(false);
  TrackCutsAntiProtonMass->SetCutSmallestSig(false);


  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(1000010020);
  PDGParticles.push_back(1000010020);

  std::vector<bool> closeRejection;
  //pairs:
  // pp             0
  // p bar p        1
  // p d            2
  // p bar d        3
  // bar p bar p    4
  // bar p d        5
  // bar p bar d    6
  // d d            7
  // d bar d        8
  // bar d bar d    9

  const int nPairs = 10;
  for (int i = 0; i < nPairs; ++i) {

    closeRejection.push_back(false);
  }

  closeRejection[0] = true;  // pp
  closeRejection[2] = true;  // pd
  closeRejection[4] = true;  // barp barp
  closeRejection[6] = true;  // barp bard
  closeRejection[7] = true;  // dd
  closeRejection[9] = true;  // bard bard



  //We need to set the ZVtx bins
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
  //The Multiplicity bins are set here
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

  //The next part is for the result histograms. The order of hist. is the following:
  //                Particle1     Particle2
  //Particle 1       Hist 1         Hist2
  //
  //Particle 2                      Hist3
  //The same way the values for binning, minimum and maximum k* range have to be set!
  //Number of bins
  std::vector<int> NBins;
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  std::vector<float> kMin;
  //minimum k* value
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  //maximum k* value
  std::vector<float> kMax;
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);

  //To put all this into the task we add it to our collection config object in
  //the following way:
  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto");
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  //Do you want to have an explicit binning of the correlation function for each multiplicity
  //bin set above?
  config->SetMultBinning(true);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetDeltaEtaMax(0.012); // and here you set the actual values
  config->SetDeltaPhiMax(0.012); // and here you set the actual values
  config->SetMixingDepth(10);
  /*
  //This is just to show off what would be possible in case you are interested, don't be confused by this at the beginning
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
      std::cout << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
    }
    if (isMC) {
      config->SetPhiEtaBinnign(true);
    } else {
      std::cout << "You are trying to request the Eta Phi Plots without MC Info; fix it wont work! \n";
    }
  }
   */
  //now we create the task
  AliAnalysisTaskFemtoDreamDeuteron *task =
    new AliAnalysisTaskFemtoDreamDeuteron("FemtoDreamDefault", isMC);
  //THIS IS VERY IMPORTANT ELSE YOU DONT PROCESS ANY EVENTS
  //kINT7 == Minimum bias
  //kHighMultV0 high multiplicity triggered by the V0 detector
  if (CentEst == "kInt7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (CentEst == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
  } else {
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
  }

  //Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetTrackCutsDeuteronDCA(TrackCutsDeuteronDCA);
  task->SetTrackCutsDeuteronMass(TrackCutsDeuteronMass);
  task->SetTrackCutsAntiDeuteronDCA(TrackCutsAntiDeuteronDCA);
  task->SetTrackCutsAntiDeuteronMass(TrackCutsAntiDeuteronMass);
  task->SetTrackCutsProtonDCA(TrackCutsProtonDCA);
  task->SetTrackCutsProtonMass(TrackCutsProtonMass);
  task->SetTrackCutsAntiProtonDCA(TrackCutsAntiProtonDCA);
  task->SetTrackCutsAntiProtonMass(TrackCutsAntiProtonMass);
  task->SetCollectionConfig(config);
  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);
  TString addon = "";
  if (CentEst == "kInt7") {
    addon += "MB";
  } else if (CentEst == "kHM") {
    addon += "HM";
  }

  AliAnalysisDataContainer *coutputQA;
  TString QAName = Form("%sQA", addon.Data());
  coutputQA = mgr->CreateContainer(
                QAName.Data(), TList::Class(),
                AliAnalysisManager::kOutputContainer,
                Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}

