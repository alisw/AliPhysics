#include "TROOT.h"
#include "TSystem.h"
AliAnalysisTaskSE* AddTaskFemtoDreamDeuteron(
    bool isMC=false,
    TString CentEst="kInt7")
{
  //Framework specific blabla
  // the manager is static, so get the existing manager via the static method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
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
  AliFemtoDreamEventCuts *evtCuts=
      AliFemtoDreamEventCuts::StandardCutsRun2();
  //This sets the method we want to use to clean up events with negative or too
  //low multiplicity. Usually you use the matching multiplicity estiamtor in your
  //event collection
  evtCuts->CleanUpMult(false,false,false,true);

  //Track Cuts are defined here
  AliFemtoDreamTrackCuts *TrackCuts1=new AliFemtoDreamTrackCuts();
  //This is used for the DCA distribution to estimate the fractions of
  //primaries, secondaries etc. At the moment the splitting of secondary contributions
  //is only done for protons!
  TrackCuts1->SetPlotDCADist(false);
  //A combined Sigma plot of the Sigma_TPC vs. Sigma_TOF in a TH3F, eats memory like
  //a student on D-Day.
  TrackCuts1->SetPlotCombSigma(false);

  //Ill only comment on the non self speaking ones!
  TrackCuts1->SetIsMonteCarlo(isMC);
  TrackCuts1->SetCutCharge(1);
  TrackCuts1->SetFilterBit(128);
  //You want to optimize this to only select regions with reasonable high purity
  TrackCuts1->SetPtRange(0.4, 3.0);
  //You want to ensure good tracking quality
  TrackCuts1->SetEtaRange(-0.8, 0.8);
  TrackCuts1->SetNClsTPC(80);
  TrackCuts1->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
  //You want to select primaries!
  TrackCuts1->SetDCAVtxZ(0.2);
  TrackCuts1->SetDCAVtxXY(0.1);
  TrackCuts1->SetCutSharedCls(true);
  TrackCuts1->SetCutTPCCrossedRows(true,70,0.83);
  //Here you define the PID
  TrackCuts1->SetPID(AliPID::kDeuteron, 0.8);
  //We are looking for pions rejecting them would be obstructive.
  TrackCuts1->SetRejLowPtPionsTOF(true);
  //this checks if the sigma of the wanted hypothesis is the smallest, and if
  //another particle has a smaller sigma, the track is rejected.
  TrackCuts1->SetCutSmallestSig(true);

  //The same things for negative pions
  AliFemtoDreamTrackCuts *TrackCuts2=new AliFemtoDreamTrackCuts();
  TrackCuts2->SetIsMonteCarlo(isMC);
  TrackCuts2->SetCutCharge(-1);
  TrackCuts2->SetFilterBit(128);
  TrackCuts2->SetPtRange(0.4, 3.0);
  TrackCuts2->SetEtaRange(-0.8, 0.8);
  TrackCuts2->SetNClsTPC(80);
  TrackCuts2->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
  TrackCuts2->SetDCAVtxZ(0.2);
  TrackCuts2->SetDCAVtxXY(0.1);
  TrackCuts2->SetCutSharedCls(true);
  TrackCuts2->SetCutTPCCrossedRows(true,70,0.83);
  TrackCuts2->SetPID(AliPID::kDeuteron, 0.8);
  TrackCuts2->SetRejLowPtPionsTOF(true);
  TrackCuts2->SetCutSmallestSig(true);

  AliFemtoDreamTrackCuts *TrackCuts3=new AliFemtoDreamTrackCuts();
  TrackCuts3->SetIsMonteCarlo(isMC);
  TrackCuts3->SetCutCharge(1);
  TrackCuts3->SetFilterBit(128);
  TrackCuts3->SetPtRange(0.5, 4.05);
  TrackCuts3->SetEtaRange(-0.8, 0.8);
  TrackCuts3->SetNClsTPC(80);
  TrackCuts3->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
  TrackCuts3->SetDCAVtxZ(0.2);
  TrackCuts3->SetDCAVtxXY(0.1);
  TrackCuts3->SetCutSharedCls(true);
  TrackCuts3->SetCutTPCCrossedRows(true,70,0.83);
  TrackCuts3->SetPID(AliPID::kProton, 0.75);
  TrackCuts3->SetRejLowPtPionsTOF(true);
  TrackCuts3->SetCutSmallestSig(true);

  AliFemtoDreamTrackCuts *TrackCuts4=new AliFemtoDreamTrackCuts();
  TrackCuts4->SetIsMonteCarlo(isMC);
  TrackCuts4->SetCutCharge(1);
  TrackCuts4->SetFilterBit(128);
  TrackCuts4->SetPtRange(0.5, 4.05);
  TrackCuts4->SetEtaRange(-0.8, 0.8);
  TrackCuts4->SetNClsTPC(80);
  TrackCuts4->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
  TrackCuts4->SetDCAVtxZ(0.2);
  TrackCuts4->SetDCAVtxXY(0.1);
  TrackCuts4->SetCutSharedCls(true);
  TrackCuts4->SetCutTPCCrossedRows(true,70,0.83);
  TrackCuts4->SetPID(AliPID::kProton, 0.75);
  TrackCuts4->SetRejLowPtPionsTOF(true);
  TrackCuts4->SetCutSmallestSig(true);

  //Now we define stuff we want for our Particle collection
  //Thanks, CINT - will not compile due to an illegal constructor
  //std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  //First we need to tell him about the particles we mix, from the
  //PDG code the mass is obtained.
  std::vector<int> PDGParticles;
  PDGParticles.push_back(211);
  PDGParticles.push_back(211);

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
  std::vector<float> kMin;
  //minimum k* value
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  //maximum k* value
  std::vector<float> kMax;
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);

  //To put all this into the task we add it to our collection config object in
  //the following way:
  AliFemtoDreamCollConfig *config=new AliFemtoDreamCollConfig("Femto","Femto");
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  //Do you want to have an explicit binning of the correlation function for each multiplicity
  //bin set above?
  config->SetMultBinning(true);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  //Here we set the mixing depth.
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
  AliAnalysisTaskFemtoDreamDeuteron *task=
      new AliAnalysisTaskFemtoDreamDeuteron("FemtoDreamDefault",isMC);
  //THIS IS VERY IMPORTANT ELSE YOU DONT PROCESS ANY EVENTS
  //kINT7 == Minimum bias
  //kHighMultV0 high multiplicity triggered by the V0 detector
  if(CentEst == "kInt7"){
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (CentEst == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
  }else{
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
  }

  //Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetTrackCutsPart1(TrackCuts1);
  task->SetTrackCutsPart2(TrackCuts2);
  task->SetTrackCutsPart3(TrackCuts3);
  task->SetTrackCutsPart4(TrackCuts4);
  task->SetCollectionConfig(config);

  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);

  AliAnalysisDataContainer *coutputQA;
  TString QAName = Form("MyTask");
  coutputQA = mgr->CreateContainer(
      QAName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}

