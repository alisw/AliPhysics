AliAnalysisTask *AddTask_oezdemir_ElectronEfficiency(Bool_t getFromAlien=kFALSE,TString cFileName = "Configpp2010ElectronEfficiency.C"){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_oezdemir_ElectronEfficiency", "No analysis manager found.");
    return 0;
  }

//Get the current train configuration
  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  TString configBasePath("$TRAIN_ROOT/oezdemir_LOWMASS/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";

  if (getFromAlien &&
      (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/m/mozdemir/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data())))
     ) {
        configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;

  //load dielectron configuration file
  gROOT->LoadMacro(configFilePath.Data());

  std::cout << "computing binning..." << std::endl;
  Double_t EtaBins[nBinsEta+1];
  for(Int_t i=0;i<=nBinsEta;i++) { EtaBins[i] = EtaMin + i*(EtaMax-EtaMin)/nBinsEta; }
  Double_t PhiBins[nBinsPhi+1];
  for(Int_t i=0;i<=nBinsPhi;i++) { PhiBins[i] = PhiMin + i*(PhiMax-PhiMin)/nBinsPhi; }

  const Int_t nBinsPt =  ( sizeof(PtBins) / sizeof(PtBins[0]) )-1;

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  std::cout << "hasMC = " << hasMC << std::endl;
 
  // Electron efficiency task
  AliAnalysisTaskElectronEfficiency *task = new AliAnalysisTaskElectronEfficiency("oezdemir_ElectronEfficiency");
  std::cout << "task created: " << task->GetName() << std::endl;
  
  Int_t triggerNames=(AliVEvent::kMB);
  Bool_t deactivateTree=kFALSE;
  Bool_t forcePhysSelAndTrigMask=kFALSE;
  //event related
  // Note: event cuts are identical for all analysis 'cutInstance's that run together!
  if (!hasMC) task->UsePhysicsSelection();
  if (!hasMC) task->SetTriggerMask(triggerNames);
  if (forcePhysSelAndTrigMask) { // should be kFALSE by default! (may be needed for new MC productions according to Mahmut)
    task->UsePhysicsSelection();
    task->SetTriggerMask(triggerNames);
  }
  task->SetEventFilter(SetupEventCuts()); //returns eventCuts from Config.
  task->SetCentralityRange(CentMin, CentMax);
  task->SetNminEleInEventForRej(NminEleInEventForRej);
  //track related
  task->SetCheckV0daughterElectron(checkV0dauEle);
  task->SetEtaRangeGEN(EtaMinGEN, EtaMaxGEN);
  task->SetPtRangeGEN(PtMinGEN, PtMaxGEN);
  //MC related
  task->SetCutInjectedSignal(CutInjectedSignals);
  //output related
  task->SetBins(nBinsPt,PtBins,nBinsEta,EtaBins,nBinsPhi,PhiBins);
  task->SetRunBins(sRuns);
  if (deactivateTree) task->SetWriteTree(kFALSE);
  else                task->SetWriteTree(writeTree);
  task->SetSupportedCutInstance(supportedCutInstance);
  task->CreateHistoGen();
  
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliAnalysisFilter *trackCuts = SetupTrackCutsAndSettings(i); // main function in config file
    if (!trackCuts) { std::cout << "WARNING: no TrackCuts given - skipping this Cutset ('"<<arrNames->At(i)->GetName()<<"')!" << std::endl; continue; }
    if (isPrefilterCutset) {
      Int_t success = SetupPrefilterPairCuts(i);
      if (!success) { std::cout << "WARNING: no/bad Prefilter PairCuts given - skipping this Cutset ('"<<arrNames->At(i)->GetName()<<"')!" << std::endl; continue; }
    }
    //
    // fill std vectors with all information which is individual per track setting:
    task->AttachTrackCuts(trackCuts);
    task->AttachDoPrefilterEff(isPrefilterCutset);
    task->AttachRejCutMee(rejCutMee);
    task->AttachRejCutTheta(rejCutTheta);
    task->AttachRejCutPhiV(rejCutPhiV);
    
    task->CreateHistograms(names,i);
  }
  
  mgr->AddTask(task);
  
  //
  // Create containers for input/output
  //
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("oezdemir_ElectronEfficiency", TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,"LMEE.root");
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("oezdemir_supportHistos", TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,"LMEE.root");
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("oezdemir_EffTree", TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer,"LMEE.root");
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("oezdemir_stats", TH1D::Class(),
                                                            AliAnalysisManager::kOutputContainer,"LMEE.root");                                                          

  //connect input/output
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);

  return task;

}//AddTask
 
