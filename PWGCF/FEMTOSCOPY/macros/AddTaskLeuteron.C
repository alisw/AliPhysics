/*
 * AddTaskLeuteron.C
 *
 *  Created on:	10 December 2019
 *	Author:	Michael Jung
 */

AliAnalysisTaskSE *AddTaskLeuteron(
  bool fullBlastQA = false,
  bool isMC = false,
  bool isHighMultV0 = true,
  bool isNanoAOD = true,
  bool BruteForceDebugging = false){

  // isHighMultV0:
  // (false)  kINT7:	    minimum bias trigger
  // (true)   kHighMultV0:  high multiplicity trigger

  if(BruteForceDebugging){
    printf("x-x-> AddTaskLeuteron: Begin of the AddTask\n");
  }

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
    // 1. (boolean) use multiplicity measured by the SPD
    // 2. (boolean) use multiplicity measured by the v0A
    // 3. (boolean) use multiplicity measured by the v0C
    // 4. (boolean) use the reference multiplicity
  
  AliFemtoDreamTrackCuts *TrackCuts1 = new AliFemtoDreamTrackCuts();

  if(!TrackCuts1){				        // check if TrackCuts1 is there
    printf("TrackCuts1 not found\n");
    return nullptr;
  }

  TrackCuts1->SetPlotDCADist(false);			// plot DCA_xy vs. pT
  TrackCuts1->SetPlotCombSigma(false);			// plot combined sigma: nSigmaTOF vs. nSigmaTPC vs. momentum
  TrackCuts1->SetIsMonteCarlo(isMC);
  TrackCuts1->SetCutCharge(1);				// set electrical charge of particle 1
  TrackCuts1->SetFilterBit(128);			// 128 is TPC only
  TrackCuts1->SetPtRange(0.5,4.05);			// set range for the transverse momentum (GeV/c)
  TrackCuts1->SetEtaRange(-0.8,0.8);			// set range of the pseudo-rapidity
  TrackCuts1->SetNClsTPC(80);				// set lower limit of clusters per track in the TPC
  TrackCuts1->SetDCAReCalculation(true);		// recalculate the DCA by PropagateToVertex or use information stored in AOD
  TrackCuts1->SetDCAVtxZ(2.0);				// DCA from track to z-coordiante of primary vertex (cm)
  TrackCuts1->SetDCAVtxXY(1.0);				// DCA from track to x-y-plane of primary vertex (cm)
  TrackCuts1->SetCutSharedCls(true);			// cut out tracks which have shared clusters in the TPC
  TrackCuts1->SetCutTPCCrossedRows(true,70,0.83);
    // SetCutTPCCrossedRows(1,2,3)
    // 1. argument (boolean) use it or not
    // 2. agrument (integer) lower limit for the number of crossed rows
    // 3. argument (float) lower limit for the fraction of crossed rows over findable clusters

  TrackCuts1->SetPID(AliPID::kProton,0.75);		// maximum momentum of the particle at its entrance point to the TPC (not pt) measured only(!) in the TPC (GeV/c)
							// above threshold use TPC + TOF 
  TrackCuts1->SetRejLowPtPionsTOF(true);		// reject pions with low transverse momentum measured in the TOF
  TrackCuts1->SetCutSmallestSig(true);			// reject tracks which have a lower sigma for other particles 
  TrackCuts1->SetMinimalBooking(false);			// set minimal booking (get only pt spectrum)

  if(BruteForceDebugging){
    printf("x-x-> AddTaskLeuteron: Cuts for the Proton (TrackCuts1) set\n");
  }

  AliFemtoDreamTrackCuts *TrackCuts2 = new AliFemtoDreamTrackCuts();

  if(!TrackCuts2){
    printf("TrackCuts2 not found\n");
    return nullptr;
  }

  TrackCuts2->SetPlotDCADist(false);
  TrackCuts2->SetPlotCombSigma(false);
  TrackCuts2->SetIsMonteCarlo(isMC);
  TrackCuts2->SetCutCharge(-1);
  TrackCuts2->SetFilterBit(128);
  TrackCuts2->SetPtRange(0.5,4.05);
  TrackCuts2->SetEtaRange(-0.8,0.8);
  TrackCuts2->SetNClsTPC(80);
  TrackCuts2->SetDCAReCalculation(true);
  TrackCuts2->SetDCAVtxZ(0.2);
  TrackCuts2->SetDCAVtxXY(0.1);
  TrackCuts2->SetCutSharedCls(true);
  TrackCuts2->SetCutTPCCrossedRows(true,70,0.83);
  TrackCuts2->SetPID(AliPID::kProton,0.75);
  TrackCuts2->SetRejLowPtPionsTOF(true);
  TrackCuts2->SetCutSmallestSig(true);
  TrackCuts2->SetMinimalBooking(false);

  if(BruteForceDebugging){
    printf("x-x-> AddTaskLeuteron: Cuts for the Antiproton (TrackCuts2) set\n");
  }

  AliFemtoDreamTrackCuts *TrackCuts3 = new AliFemtoDreamTrackCuts();

  if(!TrackCuts3){
    printf("TrackCuts3 not found\n");
    return nullptr;
  }

  TrackCuts3->SetPlotDCADist(false);
  TrackCuts3->SetPlotCombSigma(false);
  TrackCuts3->SetIsMonteCarlo(isMC);
  TrackCuts3->SetCutCharge(1);
  TrackCuts3->SetFilterBit(128);
  TrackCuts3->SetPtRange(0.8,2.5);
  TrackCuts3->SetEtaRange(-0.8,0.8);
  TrackCuts3->SetNClsTPC(70);
  TrackCuts3->SetDCAReCalculation(true);
  TrackCuts3->SetDCAVtxZ(1.0);
  TrackCuts3->SetDCAVtxXY(1.0);
  TrackCuts3->SetCutSharedCls(true);
  TrackCuts3->SetCutTPCCrossedRows(true,70,0.83);
  TrackCuts3->SetPID(AliPID::kDeuteron,1.5);
  TrackCuts3->SetRejLowPtPionsTOF(true);
  TrackCuts3->SetCutSmallestSig(true);
  TrackCuts3->SetMinimalBooking(false);

  if(BruteForceDebugging){
    printf("x-x-> AddTaskLeuteron: Cuts for the Deuteron (TrackCuts3) set\n");
  }


  AliFemtoDreamTrackCuts *TrackCuts4 = new AliFemtoDreamTrackCuts();
  
  if(!TrackCuts4){					// check if TrackCuts4 is there
    printf("TrackCuts4 not found\n");
    return nullptr;
  }

  TrackCuts4->SetPlotDCADist(false);
  TrackCuts4->SetPlotCombSigma(false);
  TrackCuts4->SetIsMonteCarlo(isMC);
  TrackCuts4->SetCutCharge(-1);
  TrackCuts4->SetFilterBit(128);
  TrackCuts4->SetPtRange(0.8,2.5);
  TrackCuts4->SetEtaRange(-0.8,0.8);
  TrackCuts4->SetNClsTPC(70);			
  TrackCuts4->SetDCAReCalculation(true);
  TrackCuts4->SetDCAVtxZ(1.0);
  TrackCuts4->SetDCAVtxXY(1.0);
  TrackCuts4->SetCutSharedCls(true);
  TrackCuts4->SetCutTPCCrossedRows(true,70,0.83);
  TrackCuts4->SetPID(AliPID::kDeuteron,1.5);
  TrackCuts4->SetRejLowPtPionsTOF(true);
  TrackCuts4->SetCutSmallestSig(true);
  TrackCuts4->SetMinimalBooking(false);

  if(BruteForceDebugging){
    printf("x-x-> AddTaskLeuteron: Cuts for the Antideuteron (TrackCuts4) set\n");
  }

  AliFemtoDreamv0Cuts *LambdaCuts5 = AliFemtoDreamv0Cuts::LambdaCuts(isMC,true,false);
    // LambdaCuts(1,2,3)
    // 1. (boolean) run over data or Monte Carlo data
    // 2. (boolean) plot CPA vs. pT (CPA: cosine of pointing angle) 
    // 3. (boolean) study cuts in detail (works only for Monte Carlo data)
  
  if(!LambdaCuts5){				    // check if LambdaCuts5 is there
    printf("LambdaCuts5 not found\n");
    return nullptr;
  }

  AliFemtoDreamTrackCuts *TrackCuts5a = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC,true,false);
    // DecayProtonCuts(1,2,3)
    // 1. (boolean) run over data or Monte Carlo data
    // 2. (boolean) enable pile-up rejection
    // 3. (boolean) study cuts in detail (works only for Monte Carlo data)
  
  if(!TrackCuts5a){				    // check if TrackCuts5a for the decay products is there
    printf("TrackCuts5a not found\n");
    return nullptr;
  }

  AliFemtoDreamTrackCuts *TrackCuts5b = AliFemtoDreamTrackCuts::DecayPionCuts(isMC,true,false);
    // DecayPionCuts(1,2,3)
    // 1. (boolean) run over data or Monte Carlo data
    // 2. (boolean) enable pile-up rejection
    // 3. (boolean) study cuts in detail (works only for Monte Carlo data)

  if(!TrackCuts5b){				    // check if TrackCuts5b for the decay product is there
    printf("TrackCuts5b not found\n");
    return nullptr;
  }

  TrackCuts5a->SetCutCharge(1);	    // Proton
  TrackCuts5b->SetCutCharge(-1);    // negative Pion
  
  LambdaCuts5->SetPosDaugterTrackCuts(TrackCuts5a); // it is "Daugter" and not "Daughter", check AliFemtoDreamv0Cuts.h
  LambdaCuts5->SetNegDaugterTrackCuts(TrackCuts5b);
  LambdaCuts5->SetPDGCodePosDaug(2212);	  // Proton
  LambdaCuts5->SetPDGCodeNegDaug(-211);	  // negative Pion
  LambdaCuts5->SetPDGCodev0(3122);	  // Lambda

  if(BruteForceDebugging){
    printf("x-x-> AddTaskLeuteron: Cuts for the Lambda (LambdaCuts5) set\n");
  }

  AliFemtoDreamv0Cuts *LambdaCuts6 = AliFemtoDreamv0Cuts::LambdaCuts(isMC,true,false);
  
  if(!LambdaCuts6){					    // check if LambdaCuts6 is there
    printf("LambdaCuts6 not found\n");
    return nullptr;
  }

  AliFemtoDreamTrackCuts *TrackCuts6a = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC,true,false);

  if(!TrackCuts6a){					    // check if TrackCuts6a for the decay product is there
    printf("TrackCuts6a not found\n");
    return nullptr;
  }

  AliFemtoDreamTrackCuts *TrackCuts6b = AliFemtoDreamTrackCuts::DecayPionCuts(isMC,true,false);

  if(!TrackCuts6b){					    // check if TrackCuts6b for the decay product is there
    printf("TrackCuts6b not found\n");
    return nullptr;
  }

  TrackCuts6a->SetCutCharge(-1);    // Antiproton
  TrackCuts6b->SetCutCharge(1);	    // positive Pion

  LambdaCuts6->SetNegDaugterTrackCuts(TrackCuts6a);
  LambdaCuts6->SetPosDaugterTrackCuts(TrackCuts6b);
  LambdaCuts6->SetPDGCodePosDaug(211);	  // positive Pion
  LambdaCuts6->SetPDGCodeNegDaug(-2212);  // Antiproton
  LambdaCuts6->SetPDGCodev0(-3122);	  // Antilambda
  
  if(BruteForceDebugging){
    printf("x-x-> AddTaskLeuteron: Cuts for the Antilambda (LambdaCuts6) set\n");
  }

  if(!fullBlastQA){
    evtCuts->SetMinimalBooking(true);
    TrackCuts1->SetMinimalBooking(true);
    TrackCuts2->SetMinimalBooking(true);
    TrackCuts3->SetMinimalBooking(true);
    TrackCuts4->SetMinimalBooking(true);
    LambdaCuts5->SetMinimalBooking(true);
    LambdaCuts6->SetMinimalBooking(true);
  }

  std::vector<int> NBins;		// number of bins
  std::vector<float> kMin;		// minimum value of k*
  std::vector<float> kMax;		// maximum value of k*
  std::vector<int> PDGParticles;	// PDG codes of the particles
  std::vector<float> ZVtxBins;		// set the number of bins for the z-component of the primary vertex
  std::vector<int> MultBins;

  const int nPairs = 21;		// number of correlation between particles

  //		    | proton | antiproton | deuteron | antideuteron | lambda | antilambda |
  //  -------------------------------------------------------------------------------------
  //  proton	    |	x    |	   x	  |	x    |	    x	    |	x    |	    x	  |  6
  //  antiproton    |	     |	   x	  |	x    |	    x	    |	x    |	    x	  |  5
  //  deuteron	    |	     |		  |	x    |	    x	    |	x    |	    x	  |  4
  //  antideuteron  |	     |		  |	     |	    x	    |	x    |	    x	  |  3
  //  lambda	    |	     |		  |	     |		    |	x    |	    x	  |  2
  //  antilambda    |	     |		  |	     |		    |	     |	    x	  |  1
  //											    --
  //											    21


  // set Nbins, kMin and kMax
  for(int i = 1;i < nPairs + 1;++i){
    NBins.push_back(750);
    kMin.push_back(0.0);
    kMax.push_back(3.0);
  }

  // set PDG codes, in the same order as PairCleaner->StoreParticle() 
  PDGParticles.push_back(2212);		// Proton
  PDGParticles.push_back(2212);		// Antiproton
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
  config->SetMultBinning(true);					  // enable explicit binning of the correlation function for each multiplicity
  config->SetMixingDepth(10);					  // the number of saved events for the event mixing
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);	  // reference multiplicity estimator

  if(BruteForceDebugging){
    printf("x-x-> AddTaskLeuteron: Values handed over to the config\n");
  }

  AliAnalysisTaskLeuteronNanoAOD *taskNanoAOD;
  AliAnalysisTaskLeuteronAOD *taskAOD;


  if(isNanoAOD){

    taskNanoAOD = new AliAnalysisTaskLeuteronNanoAOD("FemtoLeuteronNanoAOD",isMC,isHighMultV0,BruteForceDebugging);

    if(!taskNanoAOD){				  // check if the NanoAOD task is there
      printf("taskNanoAOD not found\n");
      return nullptr;
    }

    if(isHighMultV0){

      taskNanoAOD->SelectCollisionCandidates(AliVEvent::kHighMultV0);
      std::cout << "Added kHighMultV0 Trigger" << std::endl;

    } else{

      taskNanoAOD->SelectCollisionCandidates(AliVEvent::kINT7);
      std::cout << "Added kINT7 Trigger" << std::endl;

    }

    taskNanoAOD->SetEventCuts(evtCuts);
    taskNanoAOD->SetTrackCutsPart1(TrackCuts1);
    taskNanoAOD->SetTrackCutsPart2(TrackCuts2);
    taskNanoAOD->SetTrackCutsPart3(TrackCuts3);
    taskNanoAOD->SetTrackCutsPart4(TrackCuts4);
    taskNanoAOD->Setv0CutsPart5(LambdaCuts5);
    taskNanoAOD->Setv0CutsPart6(LambdaCuts6);
    taskNanoAOD->SetCollectionConfig(config);
    mgr->AddTask(taskNanoAOD);

  
  } else{

    taskAOD = new AliAnalysisTaskLeuteronAOD("FemtoLeuteronAOD",isMC,isHighMultV0,BruteForceDebugging);

    if(!taskAOD){				  // check if the AOD task is there
      printf("taskAOD not found\n");
      return nullptr;
    }

    if(isHighMultV0){

      taskAOD->SelectCollisionCandidates(AliVEvent::kHighMultV0);
      std::cout << "Added kHighMultV0 Trigger" << std::endl;

    } else{

      taskAOD->SelectCollisionCandidates(AliVEvent::kINT7);
      std::cout << "Added kINT7 Trigger" << std::endl;

    }

    taskAOD->SetEventCuts(evtCuts);
    taskAOD->SetTrackCutsPart1(TrackCuts1);
    taskAOD->SetTrackCutsPart2(TrackCuts2);
    taskAOD->SetTrackCutsPart3(TrackCuts3);
    taskAOD->SetTrackCutsPart4(TrackCuts4);
    taskAOD->Setv0CutsPart5(LambdaCuts5);
    taskAOD->Setv0CutsPart6(LambdaCuts6);
    taskAOD->SetCollectionConfig(config);
      mgr->AddTask(taskAOD);

  }



                                                                                                                         
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  if(!cinput){
    printf("cinput not found\n");
    return nullptr;
  }


  TString addon = "";
  TString suffix = "";
  TString file = AliAnalysisManager::GetCommonFileName();

  TString coutputEventCutsName = Form("%sEventCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEventCuts = mgr->CreateContainer(
    coutputEventCutsName.Data(), 
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), coutputEventCutsName.Data())
  );

  TString coutputProtonCutsName = Form("%sProtonCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputProtonCuts = mgr->CreateContainer(
    coutputProtonCutsName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), coutputProtonCutsName.Data())
  );

  TString coutputAntiprotonCutsName = Form("%sAntiprotonCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiprotonCuts = mgr->CreateContainer(
    coutputAntiprotonCutsName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), coutputAntiprotonCutsName.Data())
  );

  TString coutputDeuteronCutsName = Form("%sDeuteronCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputDeuteronCuts = mgr->CreateContainer(
    coutputDeuteronCutsName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), coutputDeuteronCutsName.Data())
  );
  
  TString coutputAntideuteronCutsName = Form("%sAntideuteronCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntideuteronCuts = mgr->CreateContainer(
    coutputAntideuteronCutsName.Data(), 
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), coutputAntideuteronCutsName.Data())
  );
  
  TString coutputLambdaCutsName = Form("%sLambdaCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputLambdaCuts = mgr->CreateContainer(
    coutputLambdaCutsName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), coutputLambdaCutsName.Data())
  );
  
  TString coutputAntilambdaCutsName = Form("%sAntilambdaCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntilambdaCuts = mgr->CreateContainer(
    coutputAntilambdaCutsName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), coutputAntilambdaCutsName.Data())
  );
  
  TString coutputPairCleanerName = Form("%sPairCleaner%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputPairCleaner = mgr->CreateContainer(
    coutputPairCleanerName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), coutputPairCleanerName.Data())
    );

  TString coutputResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputResults = mgr->CreateContainer(
    coutputResultsName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), coutputResultsName.Data())
  );

  TString coutputResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputResultsQA = mgr->CreateContainer(
    coutputResultsQAName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), coutputResultsQAName.Data())
  );



  if(isNanoAOD){
    mgr->ConnectInput(taskNanoAOD,  0, cinput);
    mgr->ConnectOutput(taskNanoAOD, 1, coutputEventCuts);
    mgr->ConnectOutput(taskNanoAOD, 2, coutputProtonCuts);
    mgr->ConnectOutput(taskNanoAOD, 3, coutputAntiprotonCuts);
    mgr->ConnectOutput(taskNanoAOD, 4, coutputDeuteronCuts);
    mgr->ConnectOutput(taskNanoAOD, 5, coutputAntideuteronCuts);
    mgr->ConnectOutput(taskNanoAOD, 6, coutputLambdaCuts);
    mgr->ConnectOutput(taskNanoAOD, 7, coutputAntilambdaCuts);
    mgr->ConnectOutput(taskNanoAOD, 8, coutputPairCleaner);
    mgr->ConnectOutput(taskNanoAOD, 9, coutputResults);
    mgr->ConnectOutput(taskNanoAOD,10, coutputResultsQA);
    
    if(BruteForceDebugging){
        printf("x-x-> AddTaskLeuteron: Only the taskNanoAOD left to return\n");
    }

    return taskNanoAOD;
    
  } else{
      mgr->ConnectInput(taskAOD,  0, cinput);
      mgr->ConnectOutput(taskAOD, 1, coutputEventCuts);
      mgr->ConnectOutput(taskAOD, 2, coutputProtonCuts);
      mgr->ConnectOutput(taskAOD, 3, coutputAntiprotonCuts);
      mgr->ConnectOutput(taskAOD, 4, coutputDeuteronCuts);
      mgr->ConnectOutput(taskAOD, 5, coutputAntideuteronCuts);
      mgr->ConnectOutput(taskAOD, 6, coutputLambdaCuts);
      mgr->ConnectOutput(taskAOD, 7, coutputAntilambdaCuts);
      mgr->ConnectOutput(taskAOD, 8, coutputPairCleaner);
      mgr->ConnectOutput(taskAOD, 9, coutputResults);
      mgr->ConnectOutput(taskAOD,10, coutputResultsQA);

      if(BruteForceDebugging){
        printf("x-x-> AddTaskLeuteron: Only the taskAOD left to return\n");
      }

      return taskAOD;    
    }




}
