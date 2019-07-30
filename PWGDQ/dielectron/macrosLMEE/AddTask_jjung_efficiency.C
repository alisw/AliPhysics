
AliAnalysisTaskElectronEfficiencyV2* AddTask_jjung_efficiency(TString name = "name", Bool_t isAOD, Bool_t getFromAlien = kFALSE, TString configFile="Config_jjung_lowmass.C", Bool_t DoCentralityCorrection = kFALSE, Int_t wagonnr = 0) {

  std::cout << "########################################\nADDTASK of ANALYSIS started\n########################################" << std::endl;

  // #########################################################
  // #########################################################
  // Configuring Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName = "AnalysisResults.root"; // create a subfolder in the file

  // #########################################################
  // #########################################################
  // Loading individual config file either local or from Alien

  // TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  TString configBasePath= "/data4/jung/localLegotrainNewEfficiency/";
  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/j/jjung/%s .",configFile.Data())))
      ) {
    configBasePath=Form("%s/",gSystem->pwd());
  }
  TString configFilePath(configBasePath+configFile);

  // Loading config and cutlib
  Bool_t err=kFALSE;
  err |= gROOT->LoadMacro(configFilePath.Data());
  if (err) { Error("AddTask_jjung_ElectronEfficiency_v2","Config(s) could not be loaded!"); return 0x0; }

  // Download resolution file (configured in your config.C)
  // if (GetResolutionFromAlien == kTRUE)
  //   std::cout << "Trying to download resolution file" << std::endl;
  //   gSystem->Exec(Form("alien_cp alien://%s .",resoFilenameFromAlien.c_str()));
  //   std::cout << "Load resolution file from AliEn" << std::endl;
  // }
  //
  // // Download centrality file (configured in your config.C)
  // if (GetCentralityFromAlien == kTRUE && !gSystem->Exec(Form("alien_cp alien://%s .",CentralityFilenameFromAlien.c_str()))){
  //   std::cout << "Load centrality file from AliEn" << std::endl;
  // }

  // #########################################################
  // #########################################################
  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(name.Data());


  // #########################################################
  // #########################################################
  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE);
 
  task->SetTriggerMask(triggerNames); 
  task->SetEventFilter(SetupEventCuts(wagonnr)); //returns eventCuts from Config.
  task->SetCentrality(centMin, centMax);
  

  // #########################################################
  // #########################################################
  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  // Do not set here your analysis pt-cuts
  task->SetMinPtGen(minGenPt);
  task->SetMaxPtGen(maxGenPt);
  task->SetMinEtaGen(minGenEta);
  task->SetMaxEtaGen(maxGenEta);

  // #########################################################
  // #########################################################
  // 4D single efficiency from pairs
  task->SetWriteLegsFromPair(WriteLegsFromPair);
  task->SetPtMinLegsFromPair(ptMinLegsFromPair);
  task->SetPtMaxLegsFromPair(ptMaxLegsFromPair);
  task->SetPtNBinsLegsFromPair(ptNBinsLegsFromPair);
  task->SetEtaMinLegsFromPair(etaMinLegsFromPair);
  task->SetEtaMaxLegsFromPair(etaMaxLegsFromPair);
  task->SetEtaNBinsLegsFromPair(etaNBinsLegsFromPair);
  task->SetPhiMinLegsFromPair(phiMinLegsFromPair);
  task->SetPhiMaxLegsFromPair(phiMaxLegsFromPair);
  task->SetPhiNBinsLegsFromPair(phiNBinsLegsFromPair);
  task->SetOpAngleMinLegsFromPair(opAngleMinLegsFromPair);
  task->SetOpAngleMaxLegsFromPair(opAngleMaxLegsFromPair);
  task->SetOpAngleNBinsLegsFromPair(opAngleNBinsLegsFromPair);


  // #########################################################
  // #########################################################
  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  task->SetKinematicCuts(minPtCut, maxPtCut, minEtaCut, maxEtaCut);

  // #########################################################
  // #########################################################
  // Set Binning
  if (usePtVector == true) {
    std::vector<double> ptBinsVec;
    for (unsigned int i = 0; i < nBinsPt+1; ++i){
      ptBinsVec.push_back(PtBins[i]);
    }
    task->SetPtBins(ptBinsVec);
  }
  else task->SetPtBinsLinear   (minPtBin,  maxPtBin, stepsPtBin);
  task->SetEtaBinsLinear  (minEtaBin, maxEtaBin, stepsEtaBin);
  task->SetPhiBinsLinear  (minPhiBin, maxPhiBin, stepsPhiBin);
  task->SetThetaBinsLinear(minThetaBin, maxThetaBin, stepsThetaBin);
  task->SetMassBinsLinear (minMassBin, maxMassBin, stepsMassBin);
  task->SetPairPtBinsLinear(minPairPtBin, maxPairPtBin, stepsPairPtBin);

  // #########################################################
  // #########################################################
  // Resolution File, If resoFilename = "" no correction is applied
  task->SetResolutionFile(resoFilename);
  task->SetResolutionFileFromAlien(resoFilenameFromAlien);
  task->SetResolutionDeltaPtBinsLinear   (DeltaMomMin, DeltaMomMax, NbinsDeltaMom);
  task->SetResolutionRelPtBinsLinear   (RelMomMin, RelMomMax, NbinsRelMom);
  task->SetResolutionEtaBinsLinear  (DeltaEtaMin, DeltaEtaMax, NbinsDeltaEta);
  task->SetResolutionPhiBinsLinear  (DeltaPhiMin, DeltaPhiMax, NbinsDeltaPhi);
  task->SetResolutionThetaBinsLinear(DeltaThetaMin, DeltaThetaMax, NbinsDeltaTheta);


  // #########################################################
  // #########################################################
  // Set centrality correction. If resoFilename = "" no correction is applied
  task->SetCentralityFile(centralityFilename);

  // #########################################################
  // #########################################################
  // Pairing related config
  task->SetDoPairing(DoPairing);
  task->SetULSandLS(DoULSLS);

  // #########################################################
  // #########################################################
  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  AddSingleLegMCSignal(task);
  AddPairMCSignal(task);


  // #########################################################
  // #########################################################
  // Adding cutsettings
  //TObjArray*  arrNames=names.Tokenize(";");
  const Int_t nDie=arrNames->GetEntriesFast();

  for (int iCut = 0; iCut < nDie; ++iCut){
    //TString cutDefinition(arrNames->At(iCut)->GetName());
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(iCut);
    task->AddTrackCuts(filter);
  }


  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("efficiency", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;
}
