AliAnalysisTask *AddTask_reichelt_ElectronEfficiency(TString configFile="Config_reichelt_ElectronEfficiency.C",
                                                     Double_t centMin=0, Double_t centMax=50,
                                                     Bool_t getFromAlien=kFALSE,
                                                     Bool_t deactivateTree=kFALSE, // enabling this has priority over 'writeTree' in config file! (enable for LEGO trains)
                                                     Char_t* outputFileName="reichelt_ElectronEfficiency.root",
                                                     Bool_t forcePhysSelAndTrigMask=kFALSE, // possibility to activate UsePhysicsSelection and SetTriggerMask for MC (may be needed for new MC productions according to Mahmut) as well as for AOD data.
                                                     Int_t triggerNames=(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral),
                                                     Int_t collCands=AliVEvent::kAny
                                                     )
{
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_reichelt_ElectronEfficiency", "No analysis manager found.");
    return 0;
  }

  // environment for testing/running at GSI:
  TString configBasePath("$TRAIN_ROOT/reichelt_lowmass/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  // typical Aliroot environment:
  if (trainRoot.IsNull()) configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";

  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/p/preichel/PWGDQ/dielectron/macrosLMEE/%s .",configFile.Data())))
      && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/p/preichel/PWGDQ/dielectron/macrosLMEE/LMEECutLib_reichelt.C ."))
      ) {
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+configFile);
  if (gSystem->Exec(Form("ls %s", configFilePath.Data()))==0) {
    std::cout << "loading config: " << configFilePath.Data() << std::endl;
    gROOT->LoadMacro(configFilePath.Data());
  } else {
    std::cout << "config not found: " << configFilePath.Data() << std::endl;
    return 0x0; // if return is not called, the job will fail instead of running without this task... (good for local tests, bad for train)
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

  Bool_t bESDANA=kFALSE; //Autodetect via InputHandler
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTaskElectronEfficiency","running on AODs.");
  }
  else if (mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()){
    ::Info("AddTaskElectronEfficiency","switching on ESD specific code, make sure ESD cuts are used.");
    bESDANA=kTRUE;
  }
  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);
  std::cout << "hasMC = " << hasMC << std::endl;

  // Electron efficiency task
  AliAnalysisTaskElectronEfficiency *task = new AliAnalysisTaskElectronEfficiency("reichelt_ElectronEfficiency");
  std::cout << "task created: " << task->GetName() << std::endl;

  //event related
  // Note: event cuts are identical for all analysis 'cutInstance's that run together!
  if (!hasMC) task->UsePhysicsSelection();
  if (!hasMC) task->SetTriggerMask(triggerNames);
  if (forcePhysSelAndTrigMask) { // should be kFALSE by default! (may be needed for new MC productions according to Mahmut)
    task->UsePhysicsSelection();
    task->SetTriggerMask(triggerNames);
  }
  task->SelectCollisionCandidates(collCands); // didnt check its meaning and effect...!
  //---
  task->SetEventFilter(SetupEventCuts(bESDANA)); //returns eventCuts from Config.
  task->SetCentralityRange(centMin, centMax);
  task->SetNminEleInEventForRej(NminEleInEventForRej);
  //track related
  task->SetEtaRangeGEN(EtaMinGEN, EtaMaxGEN);
  task->SetPtRangeGEN(PtMinGEN, PtMaxGEN);
  task->SetCalcEfficiencyPoslabel(CalcEfficiencyPoslabel);

  // resolution calculation
  task->SetCalcResolution(CalcResolution);
  if(CalcResolution || CalcEfficiencyRec) task->SetResolutionCuts(SetupTrackCutsAndSettings(99));
  if(CalcResolution) {
    task->SetDeltaMomBinning(NbinsDeltaMom,DeltaMomMin,DeltaMomMax);
    task->SetRelMomBinning(NbinsRelMom,RelMomMin,RelMomMax);
    task->SetDeltaEtaBinning(NbinsDeltaEta,DeltaEtaMin,DeltaEtaMax);
    task->SetDeltaThetaBinning(NbinsDeltaTheta,DeltaThetaMin,DeltaThetaMax);
    task->SetDeltaPhiBinning(NbinsDeltaPhi,DeltaPhiMin,DeltaPhiMax);
  }
  // resolution usage
  task->SetCalcEfficiencyGen(CalcEfficiencyGen); // will do both if Gen and Rec are set to kTRUE.
  task->SetCalcEfficiencyRec(CalcEfficiencyRec);
  if (resolutionfile.Contains("CENTRALITY")) resolutionfile.ReplaceAll("CENTRALITY",Form("%02.0f%02.0f",centMin,centMax));
  if(CalcEfficiencyRec && !resolutionfile.IsNull() &&
     (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/p/preichel/PWGDQ/dielectron/supportFiles/%s .",resolutionfile.Data()))) ){
    std::cout << "using resolution file: " << resolutionfile.Data() << std::endl;
    TFile *fRes = TFile::Open(Form("%s/%s",gSystem->pwd(),resolutionfile.Data()),"READ");
    if(bUseRelPResolution) task->SetResolutionP ((TObjArray*) fRes->Get("RelPResArr"),  kTRUE);
    else                   task->SetResolutionP ((TObjArray*) fRes->Get("DeltaPResArr"),kFALSE);
    if(bUseEtaResolution)  task->SetResolutionEta  ( (TObjArray*) fRes->Get("EtaResArr"));
    else                   task->SetResolutionTheta( (TObjArray*) fRes->Get("ThetaResArr"));
    task->SetResolutionPhi( (TObjArray*) fRes->Get("PhiEleResArr"), (TObjArray*) fRes->Get("PhiPosResArr"));
  }

  // pair efficiency
  task->SetDoPairing(doPairing);
  if(doPairing){
    task->SetKineTrackCuts(SetupTrackCutsAndSettings(100));
    //task->SetPairCuts(SetupTrackCutsAndSettings(101));
    SetupTrackCutsAndSettings(101); // this fills the pair cuts into rejCutMee,rejCutTheta,rejCutPhiV
    task->SetPairCutMee(rejCutMee);
    task->SetPairCutTheta(rejCutTheta);
    task->SetPairCutPhiV(rejCutPhiV);
    Double_t MeeBins[nBinsMee+1];
    for(Int_t i=0;i<=nBinsMee;i++) { MeeBins[i] = MeeMin + i*(MeeMax-MeeMin)/nBinsMee; }
    Double_t PteeBins[nBinsPtee+1];
    for(Int_t i=0;i<=nBinsPtee;i++) { PteeBins[i] = PteeMin + i*(PteeMax-PteeMin)/nBinsPtee; }
    task->SetBins(nBinsPt,PtBins,nBinsEta,EtaBins,nBinsPhi,PhiBins,nBinsMee,MeeBins,nBinsPtee,PteeBins);
  } else {
    task->SetBins(nBinsPt,PtBins,nBinsEta,EtaBins,nBinsPhi,PhiBins);
  }
  //output related
  task->SetRunBins(sRuns);
  if (deactivateTree) task->SetWriteTree(kFALSE);
  else                task->SetWriteTree(writeTree);
  task->SetSupportedCutInstance(supportedCutInstance);
  task->CreateHistoGen();

  // Post PID correction
  SetupITSSigmaEleCorrection(task);
  SetupTPCSigmaEleCorrection(task);

  // Monte Carlo Signals
  SetupMCSignals(task);

  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliAnalysisFilter *trackCuts = SetupTrackCutsAndSettings(i, bESDANA); // main function in config file
    if (!trackCuts) { std::cout << "WARNING: no TrackCuts given - skipping this Cutset ('"<<arrNames->At(i)->GetName()<<"')!" << std::endl; continue; }
    if (isPrefilterCutset) {
      Int_t success = SetupPrefilterPairCuts(i);
      if (!success) { std::cout << "WARNING: no/bad Prefilter PairCuts given - skipping this Cutset ('"<<arrNames->At(i)->GetName()<<"')!" << std::endl; continue; }
    }
    // fill std vectors with all information which is individual per track setting:
    task->AttachTrackCuts(trackCuts);
    task->AttachDoPrefilterEff(isPrefilterCutset);
    task->AttachExtraTrackCuts(anaFilterExtra);
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
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("reichelt_stats", TH1D::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);

  //connect input/output
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);

  return task;

}//AddTask
