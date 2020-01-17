/*
 * AddTaskLeuteron.C
 *
 *  Created on:	10 December 2019
 *	Author:	Michael Jung
 */

AliAnalysisTaskSE *AddTaskLeuteron(
  bool fullBlastQA = false,
  bool isMC = false,
  bool isNanoAOD = true,
  bool BruteForceDebugging = false,
  TString CentEst = "kHM"){
// kINT7


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if(!mgr){							  // check if the analysis manager is there
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }

  if(!mgr->GetInputEventHandler()){				  // check if the input event handler is there
    printf("Input event handler not connected!\n");
    return nullptr;
  }

  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();

  if(!evtCuts){
    printf("StandardCutsRun2 not found\n");			  // check if the event cuts are there
    return nullptr;
  }

  evtCuts->CleanUpMult(false,false,false,true);
    // CleanUpMult(1,2,3,4)
    // 1. (boolian) use multiplicity measured by the SPD
    // 2. (boolian) use multiplicity measured by the v0A
    // 3. (boolian) use multiplicity measured by the v0C
    // 4. (boolian) use the reference multiplicity

  AliFemtoDreamTrackCuts *TrackCuts1 = new AliFemtoDreamTrackCuts();

  if(!TrackCuts1){						  // check if TrackCuts1 is there
    printf("TrackCuts1 not found\n");
    return nullptr;
  }

  TrackCuts1->SetPlotDCADist(false);			// plot DCA_xy vs. pT
  TrackCuts1->SetPlotCombSigma(true);			// plot combined sigma: nSigmaTOF vs. nSigmaTPC vs. momentum
  TrackCuts1->SetIsMonteCarlo(isMC);
  TrackCuts1->SetCutCharge(1);				// set electrical charge of particle 1
  TrackCuts1->SetFilterBit(128);			// 128 is TPC only
  TrackCuts1->SetPtRange(0.8,2.5);			// set range for the transverse momentum (GeV/c)
  TrackCuts1->SetEtaRange(-0.8,0.8);			// set range of the pseudo-rapidity
  TrackCuts1->SetNClsTPC(70);				// set lower limit of clusters per track in the TPC
  TrackCuts1->SetDCAReCalculation(true);		// recalculate the DCA by PropagateToVertex or use information stored in AOD
  TrackCuts1->SetDCAVtxZ(1.0);				// DCA from track to z-coordiante of primary vertex (cm)
  TrackCuts1->SetDCAVtxXY(1.0);				// DCA from track to x-y-plane of primary vertex (cm)
  TrackCuts1->SetCutSharedCls(true);			// cut out tracks which have shared clusters in the TPC
  TrackCuts1->SetCutTPCCrossedRows(true,70,0.83);
    // SetCutTPCCrossedRows(1,2,3)
    // 1. argument (boolian) use it or not
    // 2. agrument (integer) number of crossed rows
    // 3. argument (float) fraction of crossed rows over findable clusters

  TrackCuts1->SetPID(AliPID::kDeuteron,1.5);		// maximum momentum of the particle at its entrance point to the TPC (not pt) measured only(!) in the TPC (GeV/c)
							// above threshold use TPC + TOF 
  TrackCuts1->SetRejLowPtPionsTOF(true);		// reject pions with low transverse momentum measured in the TOF
  TrackCuts1->SetCutSmallestSig(true);			// reject tracks which have a lower sigma for other particles 
  TrackCuts1->SetMinimalBooking(false);			// set minimal booking (get only pt spectrum)


  AliFemtoDreamTrackCuts *TrackCuts2 = new AliFemtoDreamTrackCuts();
  
  if(!TrackCuts2){					// check if TrackCuts2 is there
    printf("TrackCuts2 not found\n");
    return nullptr;
  }

  TrackCuts2->SetPlotDCADist(false);
  TrackCuts2->SetPlotCombSigma(false);
  TrackCuts2->SetIsMonteCarlo(isMC);
  TrackCuts2->SetCutCharge(-1);
  TrackCuts2->SetFilterBit(128);
  TrackCuts2->SetPtRange(0.8,2.5);
  TrackCuts2->SetEtaRange(-0.8,0.8);
  TrackCuts2->SetNClsTPC(70);			
  TrackCuts2->SetDCAReCalculation(true);
  TrackCuts2->SetDCAVtxZ(1.0);
  TrackCuts2->SetDCAVtxXY(1.0);
  TrackCuts2->SetCutSharedCls(true);
  TrackCuts2->SetCutTPCCrossedRows(true,70,0.83);
  TrackCuts2->SetPID(AliPID::kDeuteron,1.5);
  TrackCuts2->SetRejLowPtPionsTOF(true);
  TrackCuts2->SetCutSmallestSig(true);
  TrackCuts2->SetMinimalBooking(false);


  AliFemtoDreamv0Cuts *LambdaCuts3 = AliFemtoDreamv0Cuts::LambdaCuts(isMC,true,false);
    // LambdaCuts(1,2,3)
    // 1. (boolian) run over data or Monte Carlo data
    // 2. (boolian) plot CPA vs. pT (CPA: cosine of pointing angle) 
    // 3. (boolian) study cuts in detail (works only for Monte Carlo data)
  
  if(!LambdaCuts3){				    // check if LambdaCuts3 is there
    printf("LambdaCuts3 not found\n");
    return nullptr;
  }

  AliFemtoDreamTrackCuts *TrackCuts3a = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC,true,false);
    // DecayProtonCuts(1,2,3)
    // 1. (boolian) run over data or Monte Carlo data
    // 2. (boolian) enable pile-up rejection
    // 3. (boolian) study cuts in detail (works only for Monte Carlo data)
  
  if(!TrackCuts3a){				    // check if TrackCuts3a for the decay products is there
    printf("TrackCuts3a not found\n");
    return nullptr;
  }

  AliFemtoDreamTrackCuts *TrackCuts3b = AliFemtoDreamTrackCuts::DecayPionCuts(isMC,true,false);
    // DecayPionCuts(1,2,3)
    // 1. (boolian) run over data or Monte Carlo data
    // 2. (boolian) enable pile-up rejection
    // 3. (boolian) study cuts in detail (works only for Monte Carlo data)

  if(!TrackCuts3b){				    // check if TrackCuts3b for the decay product is there
    printf("TrackCuts3b not found\n");
    return nullptr;
  }

  TrackCuts3a->SetCutCharge(1);	    // Proton
  TrackCuts3b->SetCutCharge(-1);    // negative Pion
  
  LambdaCuts3->SetPosDaugterTrackCuts(TrackCuts3a); // it is "Daugter" and not "Daughter", check AliFemtoDreamv0Cuts.h
  LambdaCuts3->SetNegDaugterTrackCuts(TrackCuts3b);
  LambdaCuts3->SetPDGCodePosDaug(2212);	  // Proton
  LambdaCuts3->SetPDGCodeNegDaug(-211);	  // negative Pion
  LambdaCuts3->SetPDGCodev0(3122);	  // Lambda


  AliFemtoDreamv0Cuts *LambdaCuts4 = AliFemtoDreamv0Cuts::LambdaCuts(isMC,true,false);
  
  if(!LambdaCuts4){					    // check if LambdaCuts4 is there
    printf("LambdaCuts4 not found\n");
    return nullptr;
  }

  AliFemtoDreamTrackCuts *TrackCuts4a = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC,true,false);

  if(!TrackCuts4a){					    // check if TrackCuts4a for the decay product is there
    printf("TrackCuts4a not found\n");
    return nullptr;
  }

  AliFemtoDreamTrackCuts *TrackCuts4b = AliFemtoDreamTrackCuts::DecayPionCuts(isMC,true,false);

  if(!TrackCuts4b){					    // check if TrackCuts4b for the decay product is there
    printf("TrackCuts4b not found\n");
    return nullptr;
  }

  TrackCuts4a->SetCutCharge(-1);    // Antiproton
  TrackCuts4b->SetCutCharge(1);	    // positive Pion

  LambdaCuts4->SetNegDaugterTrackCuts(TrackCuts4a);
  LambdaCuts4->SetPosDaugterTrackCuts(TrackCuts4b);
  LambdaCuts4->SetPDGCodePosDaug(211);	  // positive Pion
  LambdaCuts4->SetPDGCodeNegDaug(-2212);  // Antiproton
  LambdaCuts4->SetPDGCodev0(-3122);	  // Antilambda

  if(!fullBlastQA){
    evtCuts->SetMinimalBooking(true);
    TrackCuts1->SetMinimalBooking(true);
    TrackCuts2->SetMinimalBooking(true);
    LambdaCuts3->SetMinimalBooking(true);
    LambdaCuts4->SetMinimalBooking(true);
  }

  std::vector<int> NBins;		// number of bins
  std::vector<float> kMin;		// minimum value of k*
  std::vector<float> kMax;		// maximum value of k*
  std::vector<int> PDGParticles;	// PDG codes of the particles
  std::vector<float> ZVtxBins;		// set the number of bins for the z-component of the primary vertex
  std::vector<int> MultBins;

  const int nPairs = 10;		// number of correlation between particles
    // 1: Deuteron - Deuteron
    // 2: Deuteron - Lambda
    // 3: Lambda - Lambda

  // set Nbins, kMin and kMax
  for(int i = 1;i < nPairs + 1;++i){
    NBins.push_back(750);
    kMin.push_back(0.0);
    kMax.push_back(3.0);
  }

  // set PDG codes, in the same order as PairCleaner->StoreParticle() 
  PDGParticles.push_back(1000010020);	// Deuteron
  PDGParticles.push_back(1000010020);	// Antideuteron
  PDGParticles.push_back(3312);		// Lambda
  PDGParticles.push_back(3312);		// Antilambda

  // set ZVtxBins
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

  // set MultBins
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

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto","Femto");
  
  if(!config){
    printf("config not found\n");
    return nullptr;
  }

  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetPDGCodes(PDGParticles);
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetMultBinning(true);		  // enable explicit binning of the correlation function for each multiplicity
  config->SetMixingDepth(10);		  // the number of saved events for the event mixing
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);	  // reference multiplicity estimator


  if(BruteForceDebugging){
    printf("NBins, kMin, kMax, PDGParticles, ZVtxBins and MultBins added to config\n");
  }

  AliAnalysisTaskLeuteronNanoAOD *taskNanoAOD;
  AliAnalysisTaskLeuteronAOD *taskAOD;

  if(isNanoAOD){

    taskNanoAOD = new AliAnalysisTaskLeuteronNanoAOD("FemtoLeuteronNanoAOD",isMC);

    if(!taskNanoAOD){				  // check if the NanoAOD task is there
      printf("taskNanoAOD not found\n");
      return nullptr;
    }

    if(CentEst == "kInt7"){				  // kInt7: minimum bias
      taskNanoAOD->SelectCollisionCandidates(AliVEvent::kINT7);
      std::cout << "Added kINT7 Trigger" << std::endl;
    } else if(CentEst == "kHM"){				  // kHighMultV0: high multiplicity triggered by the V0 detector
      taskNanoAOD->SelectCollisionCandidates(AliVEvent::kHighMultV0);
      std::cout << "Added kHighMultV0 Trigger" << std::endl;
    } else{
      std::cout << "=====================================================================" << std::endl;
      std::cout << "=====================================================================" << std::endl;
      std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
      std::cout << "=====================================================================" << std::endl;
      std::cout << "=====================================================================" << std::endl;
    }

    taskNanoAOD->SetEventCuts(evtCuts);
    taskNanoAOD->SetTrackCutsPart1(TrackCuts1);
    taskNanoAOD->SetTrackCutsPart2(TrackCuts2);
    taskNanoAOD->Setv0CutsPart3(LambdaCuts3);
    taskNanoAOD->Setv0CutsPart4(LambdaCuts4);
    taskNanoAOD->SetCollectionConfig(config);
    mgr->AddTask(taskNanoAOD);

    if(BruteForceDebugging){
      printf("evtCuts, TrackCuts1, TrackCuts2, LambdaCuts3, LambdaCuts4 and config added to the NanoAOD task\n");
    }
  
  } else{

    taskAOD = new AliAnalysisTaskLeuteronAOD("FemtoLeuteronAOD",isMC);

    if(!taskAOD){				  // check if the AOD task is there
      printf("taskAOD not found\n");
      return nullptr;
    }

    if(CentEst == "kInt7"){				  // kInt7: minimum bias
      taskAOD->SelectCollisionCandidates(AliVEvent::kINT7);
      std::cout << "Added kINT7 Trigger" << std::endl;
    } else if(CentEst == "kHM"){				  // kHighMultV0: high multiplicity triggered by the V0 detector
      taskAOD->SelectCollisionCandidates(AliVEvent::kHighMultV0);
      std::cout << "Added kHighMultV0 Trigger" << std::endl;
    } else{
      std::cout << "=====================================================================" << std::endl;
      std::cout << "=====================================================================" << std::endl;
      std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
      std::cout << "=====================================================================" << std::endl;
      std::cout << "=====================================================================" << std::endl;
    }

    taskAOD->SetEventCuts(evtCuts);
    taskAOD->SetTrackCutsPart1(TrackCuts1);
    taskAOD->SetTrackCutsPart2(TrackCuts2);
    taskAOD->Setv0CutsPart3(LambdaCuts3);
    taskAOD->Setv0CutsPart4(LambdaCuts4);
    taskAOD->SetCollectionConfig(config);
    mgr->AddTask(taskAOD);

    if(BruteForceDebugging){
      printf("evtCuts, TrackCuts1, TrackCuts2, LambdaCuts3, LambdaCuts4 and config added to the AOD task\n");
    }




  }


    

  if(BruteForceDebugging){
    printf("task added to the analysis manager\n");
  }

  TString file = AliAnalysisManager::GetCommonFileName();
                                                                                                                         
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  if(!cinput){
    printf("cinput not found\n");
    return nullptr;
  }

  AliAnalysisDataContainer *coutputQA;
  TString QAName = Form("MyTask");
  coutputQA = mgr->CreateContainer(
    QAName.Data(), TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), QAName.Data()));

  if(BruteForceDebugging){
    printf("coutputQA container created\n");
  }

  if(isNanoAOD){
    mgr->ConnectInput(taskNanoAOD, 0, cinput);
    mgr->ConnectOutput(taskNanoAOD, 1, coutputQA);
    return taskNanoAOD;
  } else{
    mgr->ConnectInput(taskAOD, 0, cinput);
    mgr->ConnectOutput(taskAOD, 1, coutputQA);

    if(BruteForceDebugging){
      printf("taskAOD connected to Input and Output\n");
    }
    return taskAOD;

    if(BruteForceDebugging){
      printf("taskAOD returned\n");
    }

  }
 
  if(BruteForceDebugging){
    printf("End of AddTask\n");
  }

}
