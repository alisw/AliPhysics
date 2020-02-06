AliAnalysisTaskSE* AddTaskFemtoDreamSigPi(
    bool isMC=true,
    TString CentEst="kHM", const char *cutVar = "0") {

  TString suffix = TString::Format("%s", cutVar);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  AliFemtoDreamEventCuts *evtCuts=
      AliFemtoDreamEventCuts::StandardCutsRun2();
  //This sets the method we want to use to clean up events with negative or too
  //low multiplicity. Usually you use the matching multiplicity estiamtor in your
  //event collection
  evtCuts->CleanUpMult(false,false,false,true);
  evtCuts->SetZVtxPosition(-10., 10.);

  //Track Cuts are defined here, no cuts other than MC since we only want to know how many are particles are produced and generate the SE distributions.
  //positive pions
  AliFemtoDreamTrackCuts *fTrackCutsPosPion=new AliFemtoDreamTrackCuts();
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  //The same things for negative pions
  AliFemtoDreamTrackCuts *fTrackCutsNegPion=new AliFemtoDreamTrackCuts();
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  //The same things for pion0
  AliFemtoDreamTrackCuts *fTrackCutsPion0=new AliFemtoDreamTrackCuts();
  fTrackCutsPion0->SetIsMonteCarlo(isMC);
  //The same things for sigma+
  AliFemtoDreamTrackCuts *fTrackCutsSigmaPlus=new AliFemtoDreamTrackCuts();
  fTrackCutsSigmaPlus->SetIsMonteCarlo(isMC);
  //The same things for sigma-
  AliFemtoDreamTrackCuts *fTrackCutsSigmaMinus=new AliFemtoDreamTrackCuts();
  fTrackCutsSigmaMinus->SetIsMonteCarlo(isMC);
  //The same things for sigma0
  AliFemtoDreamTrackCuts *fTrackCutsSigma0=new AliFemtoDreamTrackCuts();
  fTrackCutsSigma0->SetIsMonteCarlo(isMC);
 

  //Now we define stuff we want for our Particle collection
  //Thanks, CINT - will not compile due to an illegal constructor
  //std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  //First we need to tell him about the particles we mix, from the
  //PDG code the mass is obtained.
  //Order must mymik the order as the particles are added to the PairCleaner
  std::vector<int> PDGParticles;
  PDGParticles.push_back(211);  // pi+
  PDGParticles.push_back(-211); // pi-
  PDGParticles.push_back(111);  // pi0
  PDGParticles.push_back(3222); // Sig+
  PDGParticles.push_back(3112); // Sig-
  PDGParticles.push_back(3212); // Sig0
 

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
  kMax.push_back(3.);
  //pair rejection
  std::vector<bool> closeRejection;
  closeRejection.push_back(false);
  closeRejection.push_back(false);  
  closeRejection.push_back(false);
  closeRejection.push_back(false);
  closeRejection.push_back(false);  
  closeRejection.push_back(false);
  closeRejection.push_back(false);
  closeRejection.push_back(false);  
  closeRejection.push_back(false);
  closeRejection.push_back(false);
  closeRejection.push_back(false);  
  closeRejection.push_back(false);
  closeRejection.push_back(false);
  closeRejection.push_back(false);  
  closeRejection.push_back(false);
  closeRejection.push_back(false);
  closeRejection.push_back(false);  
  closeRejection.push_back(false);
  closeRejection.push_back(false);
  closeRejection.push_back(false);  
  closeRejection.push_back(false);
  //QA plots for tracks
  std::vector<int> pairQA;
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);
  pairQA.push_back(11);

  //To put all this into the task we add it to our collection config object in
  //the following way:
  AliFemtoDreamCollConfig *config=new AliFemtoDreamCollConfig("Femto","Femto");
  config->SetZBins(ZVtxBins);
  //Do you want to have an explicit binning of the correlation function for each multiplicity
  //bin set above?
  config->SetMultBinning(false);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetExtendedQAPairs(pairQA);
  //Here we set the mixing depth.
  config->SetMixingDepth(10); // AN
  config->SetMinimalBookingME(false);
  
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
  
  //now we create the task
  AliAnalysisTaskFemtoDreamSigPi *task=
      new AliAnalysisTaskFemtoDreamSigPi("FemtoDreamDefaultSigPi",isMC);
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
  task->SetTrackCutsPosPion(fTrackCutsPosPion);
  task->SetTrackCutsNegPion(fTrackCutsNegPion);
  task->SetTrackCutsPion0(fTrackCutsPion0);
  task->SetTrackCutsSigmaPlus(fTrackCutsSigmaPlus);
  task->SetTrackCutsSigmaMinus(fTrackCutsSigmaMinus);
  task->SetTrackCutsSigma0(fTrackCutsSigma0);
  task->SetCollectionConfig(config);
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
  coutputQA = mgr->CreateContainer(
    QAName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}
