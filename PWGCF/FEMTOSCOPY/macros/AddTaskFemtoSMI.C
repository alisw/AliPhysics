#include "TROOT.h"
#include "TSystem.h"
// additional inclusions following: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/macros/AddTaskFemtoLoton.C
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskFemtoSMI.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"

AliAnalysisTaskSE* AddTaskFemtoSMI(
    bool isMC = false, // for now ignore implementing MC options
    TString CentEst = "kInt7",
    bool ControlLambdaCuts = false, // this and following not in older version
    int    taskCounter = 0,
    double InvMassCut = 0.004,
    double CPACut = 0.99,
    double KaonRejectionLowerCut = 0.48,
    double KaonRejectionUpperCut = 0.515,
    double DaugToPrimCut = 0.05,
    double DaugTov0Cut = 1.5
    )
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

  // covering the protons
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);
  // leaving out anti-particles for now
  
  //Lambda Cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  
  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);  //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);  //Lambda
  //leaving out anti-lambda def's for now  

  if(ControlLambdaCuts) {
    //Adding additional settings to access from train setup
    v0Cuts->SetCutInvMass(InvMassCut);  // default argument: 0.004
    v0Cuts->SetCutCPA(CPACut);  // default: 0.99
    v0Cuts->SetKaonRejection(KaonRejectionLowerCut,KaonRejectionUpperCut);  // default: 0.48, 0.515
    v0Cuts->SetCutDCADaugToPrimVtx(DaugToPrimCut);  // default: 0.05
    v0Cuts->SetCutDCADaugTov0Vtx(DaugTov0Cut);  // default: 1.5
  }

  //Now we define stuff we want for our Particle collection
  //Thanks, CINT - will not compile due to an illegal constructor
  //std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  //First we need to tell him about the particles we mix, from the
  //PDG code the mass is obtained.
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3122);

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
  AliFemtoDreamCollConfig *config=new AliFemtoDreamCollConfig("Femto","Femto"); //IC: third option "false"
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
  //Added config:
  

  //config->SetClosePairRejection(false); // needed??
  config->SetDeltaEtaMax(0.017); // ??
  config->SetDeltaPhiMax(0.017); // ??
  // config->SetExtendedQAPairs(pairQA); 
  // config->SetmTBins(mTBins);
  // config->SetDomTMultBinning(true);
  // config->SetUseEventMixing(true);
  // config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  // need something from the "Systematic"-section in FemtoLoton.C??
  //End of added config

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
    config->SetZBins(ZVtxBins); // already included above

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
  AliAnalysisTaskFemtoSMI *task=
      new AliAnalysisTaskFemtoSMI("FemtoDreamPLDefault",isMC);
  //THIS IS VERY IMPORTANT ELSE YOU DONT PROCESS ANY EVENTS
  //kINT7 == Minimum bias
  //kHighMultV0 high multiplicity triggered by the V0 detector
  if(CentEst == "kInt7"){
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    task->SetTrigger(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (CentEst == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    task->SetTrigger(AliVEvent::kHighMultV0);
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
  // task->SetTrackCutsPart1(TrackCuts1);
  // task->SetTrackCutsPart2(TrackCuts2);
  task->SetCollectionConfig(config);
  // new additions for task:
  task->SetProtonCuts(TrackCuts);
  task->Setv0Cuts(v0Cuts);
  // needed ?? : task->SetCorrelationConfig(config); // We'll see, I guess
 

  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);

  AliAnalysisDataContainer *coutputQA;
  TString QAName;
  if(taskCounter > 0) 
	QAName = Form("MyTask%d",taskCounter);
  else
	QAName = Form("MyTask");

  coutputQA = mgr->CreateContainer(
      QAName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}

