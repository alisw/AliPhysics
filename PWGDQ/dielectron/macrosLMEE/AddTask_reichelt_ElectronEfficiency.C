AliAnalysisTask *AddTask_reichelt_ElectronEfficiency(Char_t* outputFileName="LMEEoutput.root", 
                                                     Bool_t getFromAlien=kFALSE,
                                                     Bool_t deactivateTree=kFALSE // enabling this has priority over 'writeTree' in config file! (enable for LEGO trains)
                                                     )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_reichelt_ElectronEfficiency", "No analysis manager found.");
    return 0;
  }
  
  TString configBasePath("$TRAIN_ROOT/reichelt_lowmass/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  
  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/p/preichel/PWGDQ/dielectron/macrosLMEE/Config_reichelt_ElectronEfficiency.C ."))
      && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/p/preichel/PWGDQ/dielectron/macrosLMEE/LMEECutLib_reichelt.C ."))
      // Task files (.cxx & .h) not possible because they need to be compiled...
      ) {
    configBasePath=Form("%s/",gSystem->pwd());
  }
  
  TString configFile("Config_reichelt_ElectronEfficiency.C");
  TString configFilePath(configBasePath+configFile);
  if (gSystem->Exec(Form("ls %s", configFilePath.Data()))==0) {
    std::cout << "loading config: " << configFilePath.Data() << std::endl;
    gROOT->LoadMacro(configFilePath.Data());
  } else {
    std::cout << "config not found: " << configFilePath.Data() << std::endl;
    return 0; // if return is not called, the job will fail instead of running wihout this task... (good for local tests, bad for train)
  }
  TString configLMEECutLib("LMEECutLib_reichelt.C");
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);
  if (gSystem->Exec(Form("ls %s", configLMEECutLibPath.Data()))==0) {
    std::cout << "loading config: " << configLMEECutLibPath.Data() << std::endl;
    gROOT->LoadMacro(configLMEECutLibPath.Data());
  } else std::cout << "config not found: " << configLMEECutLibPath.Data() << std::endl;
  
  std::cout << "computing binning..." << std::endl;
  Double_t EtaBins[nBinsEta+1];
  for(Int_t i=0;i<=nBinsEta;i++) { EtaBins[i] = EtaMin + i*(EtaMax-EtaMin)/nBinsEta; }
  Double_t PhiBins[nBinsPhi+1];
  for(Int_t i=0;i<=nBinsPhi;i++) { PhiBins[i] = PhiMin + i*(PhiMax-PhiMin)/nBinsPhi; }

  const Int_t nBinsPt =  ( sizeof(PtBins) / sizeof(PtBins[0]) )-1;
  
  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);
  std::cout << "hasMC = " << hasMC << std::endl;

  // Electron efficiency task
  AliAnalysisTaskElectronEfficiency *task = new AliAnalysisTaskElectronEfficiency("reichelt_ElectronEfficiency");
  std::cout << "task created: " << task->GetName() << std::endl;
  //event related
  task->SetMC(hasMC);
  task->SetRequireVertex(reqVertex);
  task->SetMaxVertexZ(vertexZcut);
  task->SetCentralityRange(CentMin, CentMax);
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
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("reichelt_ElectronEfficiency", TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("reichelt_supportHistos", TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("reichelt_EffTree", TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);

  //connect input/output
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);

  return task;

}//AddTask
