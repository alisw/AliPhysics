AliAnalysisTask *AddTask_jjung_ElectronEfficiency(Bool_t getFromAlien=kFALSE,
                                                     TString cFileName = "Config_jjung_ElectronEfficiency.C",
                                                     Char_t* outputFileName="LMEE.root",
                                                     Bool_t deactivateTree=kFALSE, // enabling this has priority over 'writeTree'! (needed for LEGO trains)
                                                     TString resolutionfile = ""
                                                     )
{

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jjung_ElectronEfficiency", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  //TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  TString configBasePath= "/data4/jung/localLegotrainEfficiency/";
  if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/j/jjung/%s .",cFileName.Data()))) ){
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);
  std::cout << "Configpath:  " << configFilePath << std::endl;
  if (gSystem->Exec(Form("ls %s", configFilePath.Data()))==0) {
    std::cout << "loading config: " << configFilePath.Data() << std::endl;
    gROOT->LoadMacro(configFilePath.Data());
  } else {
    std::cout << "config not found: " << configFilePath.Data() << std::endl;
    return 0; // if return is not called, the job will fail instead of running wihout this task... (good for local tests, bad for train)
  }

  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);
  std::cout << "hasMC = " << hasMC << std::endl;


  std::cout << "computing binning..." << std::endl;
  Double_t EtaBins[nBinsEta+1];
  for(Int_t i=0;i<=nBinsEta;i++) { EtaBins[i] = EtaMin + i*(EtaMax-EtaMin)/nBinsEta; }
  Double_t PhiBins[nBinsPhi+1];
  for(Int_t i=0;i<=nBinsPhi;i++) { PhiBins[i] = PhiMin + i*(PhiMax-PhiMin)/nBinsPhi; }

  const Int_t nBinsPt = ( sizeof(PtBins) / sizeof(PtBins[0]) )-1;


  // Electron efficiency task
  AliAnalysisTaskElectronEfficiency *task = new AliAnalysisTaskElectronEfficiency("jjung_ElectronEfficiency");
  std::cout << "task created: " << task->GetName() << std::endl;

  if(CalcEfficiencyRec && !resolutionfile.IsNull() &&
     (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/j/jjung/supportFiles/%s .",resolutionfile.Data()))) ){
    TFile *fRes = TFile::Open(Form("%s/%s",gSystem->pwd(),resolutionfile.Data()),"READ");
    if(bUseRelPResolution) task->SetResolutionP ((TObjArray*) fRes->Get("RelPResArr"),  kTRUE);
    else                   task->SetResolutionP ((TObjArray*) fRes->Get("DeltaPResArr"),kFALSE);
    if(bUseEtaResolution)  task->SetResolutionEta  ((TObjArray*) fRes->Get("EtaResArr"));
    else                   task->SetResolutionTheta((TObjArray*) fRes->Get("ThetaResArr"));
    task->SetResolutionPhi( (TObjArray*) fRes->Get("PhiEleResArr"), (TObjArray*) fRes->Get("PhiPosResArr"));
  }

  if(CalcEfficiencyRec) task->SetResolutionCuts(SetupTrackCutsAndSettings(-1));
  task->SetCalcEfficiencyRec(CalcEfficiencyRec);
  task->SetCalcEfficiencyPoslabel(CalcEfficiencyPoslabel);
  SetupMCSignals(task);
  //event related
  task->SetEventFilter(SetupEventCuts()); //returns eventCuts from Config. //cutlib->GetEventCuts(LMEECutLib::kPbPb2011_TPCTOF_Semi1)
  task->SetCentralityRange(CentMin, CentMax);
  task->SetNminEleInEventForRej(NminEleInEventForRej);
  //track related
  task->SetEtaRangeGEN(EtaMinGEN, EtaMaxGEN);
  task->SetPtRangeGEN(PtMinGEN, PtMaxGEN);
  //MC related
  //task->SetCutInjectedSignal(CutInjectedSignals);
  // resolution calculation
  task->SetCalcResolution(CalcResolution);
  if(CalcResolution) task->SetResolutionCuts(SetupTrackCutsAndSettings(-1));
  task->SetDeltaMomBinning(NbinsDeltaMom,DeltaMomMin,DeltaMomMax);
  task->SetRelMomBinning(NbinsRelMom,RelMomMin,RelMomMax);
  task->SetDeltaEtaBinning(NbinsDeltaEta,DeltaEtaMin,DeltaEtaMax);
  task->SetDeltaThetaBinning(NbinsDeltaTheta,DeltaThetaMin,DeltaThetaMax);
  task->SetDeltaPhiBinning(NbinsDeltaPhi,DeltaPhiMin,DeltaPhiMax);
  task->SetDeltaAngleBinning(NbinsDeltaAngle,DeltaAngleMin,DeltaAngleMax);
  // pair efficiency
  if(doPairing){
    task->SetKineTrackCuts(SetupTrackCutsAndSettings(900));
    //task->SetPairCuts(SetupTrackCutsAndSettings(101));
    SetupTrackCutsAndSettings(901);
    task->SetPairCutMee(rejCutMee);
    task->SetPairCutTheta(rejCutTheta);
    task->SetPairCutPhiV(rejCutPhiV);
    task->SetBins(nBinsPt,PtBins,nBinsEta,EtaBins,nBinsPhi,PhiBins,nBinsMee,MeeBins,nBinsPtee,PteeBins);
    task->SetDoPairing(kTRUE);
  }
  else
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
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("jjung_ElectronEfficiency", TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("jjung_supportHistos", TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("jjung_EffTree", TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("jjung_stats", TH1D::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);

  //connect input/output
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);

  return task;

}//AddTask
