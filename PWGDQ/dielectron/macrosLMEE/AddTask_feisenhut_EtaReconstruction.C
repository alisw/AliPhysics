AliAnalysisTaskEtaReconstruction* AddTask_feisenhut_EtaReconstruction(TString name = "name",
                                                                Bool_t isAOD,
                                                                Bool_t getFromAlien = kFALSE,
                                                                TString configFile="Config_feisenhut_EtaReconstruction.C",
                                                                Bool_t DoCentralityCorrection = kFALSE,
                                                                Bool_t cutlibPreloaded = kFALSE,
                                                                Int_t wagonnr = 0,
                                                                Int_t centrality = 4) {


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
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/f/feisenhu/PWGDQ/dielectron/macrosLMEE/%s .",configFile.Data())))
      && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/f/feisenhu/PWGDQ/dielectron/macrosLMEE/LMEECutLib_feisenhut.C ."))
      ) {
    configBasePath=Form("%s/",gSystem->pwd());
  }
  TString configFilePath(configBasePath+configFile);
  TString configLMEECutLib("LMEECutLib_feisenhut.C");
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);

  // Loading config and cutlib
  // err |= gROOT->LoadMacro(configLMEECutLibPath.Data());
  // err |= gROOT->LoadMacro(configFilePath.Data());
  // if (err) { Error("AddTask_caklein_ElectronEfficiency_v2","Config(s) could not be loaded!"); return 0x0; }

  Bool_t err=kFALSE;
  if (!cutlibPreloaded) { // should not be needed but seems to be...
    std::cout << "Cutlib was not preloaded" << std::endl;
    err |= gROOT->LoadMacro(configLMEECutLibPath.Data());
    err |= gROOT->LoadMacro(configFilePath.Data());
  }
  else{
    std::cout << "Cutlib was preloaded in a previous task" << std::endl;
  }


  // #########################################################
  // #########################################################
  // Creating an instance of the task
  AliAnalysisTaskEtaReconstruction* task = new AliAnalysisTaskEtaReconstruction(Form("%s%d",name.Data(), wagonnr));

  // #########################################################
  // #########################################################
  // Possibility to set generator. If nothing set all generators are taken into account
  // task->SetGeneratorName(generatorName);
  task->SetGeneratorMCSignalName(generatorNameForMCSignal);
  task->SetGeneratorULSSignalName(generatorNameForULSSignal);

  // #########################################################
  // #########################################################
  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE);
  task->SetTriggerMask(triggerNames);
  task->SetEventFilter(SetupEventCuts(isAOD)); //returns eventCuts from Config.

  double centMin = 0.;
  double centMax = 60.;
  GetCentrality(centrality, centMin, centMax);
  std::cout << "CentMin = " << centMin << "  CentMax = " << centMax << std::endl;
  task->SetCentrality(centMin, centMax);

  // #########################################################
  // #########################################################
  // Set debug variable for debugging code
  task->SetDebug(debug);

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
  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  task->SetKinematicCuts(minPtCut, maxPtCut, minEtaCut, maxEtaCut);

  // #########################################################
  // #########################################################
  // Set mass cuts for primary and secondary pairs
  task->SetMassCut(DoMassCut);
  task->SetPhotonMass(photonMass);
  task->SetUpperMassCutPrimaries(upperMassCutPrimaries);
  task->SetLowerMassCutPrimaries(lowerMassCutPrimaries);
  task->SetUpperPreFilterMass(upperPreFilterMass);
  task->SetLowerPreFilterMass(lowerPreFilterMass);
  task->SetMassCutSecondaries(massCutSecondaries);



  // #########################################################
  // #########################################################
  // Set Binning
  if (usePtVector == true) {
    std::vector<double> ptBinsVec;
    for (unsigned int i = 0; i < nBinsPt+1; ++i){
      ptBinsVec.push_back(ptBins[i]);
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
  // Set MCSignal and Cutsetting to fill the support histograms
  task->SetSupportHistoMCSignalAndCutsetting(nMCSignal, nCutsetting);

  // #########################################################
  // #########################################################
  // Set Cocktail weighting
  task->SetDoCocktailWeighting(DoCocktailWeighting);
  task->SetCocktailWeighting(CocktailFilename);
  task->SetCocktailWeightingFromAlien(CocktailFilenameFromAlien);

  // #########################################################
  // #########################################################
  // Pairing related config
  task->SetDoPairing(DoPairing);
  task->SetDoFourPairing(DoFourPairing);
  // task->SetULSandLS(DoULSLS);
  task->SetUsePreFilter(UsePreFilter);

  // #########################################################
  // #########################################################
  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  AddSinglePrimaryLegMCSignal(task);
  AddSingleSecondaryLegMCSignal(task);
  AddPrimaryPairMCSignal(task);
  AddSecondaryPairMCSignal(task);
  AddFourPairMCSignal(task);
  std::vector<bool> PrimaryDielectronsPairNotFromSameMother = AddSinglePrimaryLegMCSignal(task);
  std::vector<bool> SecondaryDielectronsPairNotFromSameMother = AddSingleSecondaryLegMCSignal(task);
  task->AddMCSignalsWherePrimaryDielectronPairNotFromSameMother(PrimaryDielectronsPairNotFromSameMother);
  task->AddMCSignalsWhereSecondaryDielectronPairNotFromSameMother(SecondaryDielectronsPairNotFromSameMother);


  // #########################################################
  // #########################################################
  // Set mean and width correction for ITS, TPC and TOF


  // #########################################################
  // #                set track cutsettings                  #
  // #########################################################
  // Adding standard primary electron track cutsettings
  TObjArray*  arrNames_prim=names_Prim_Track_standard_Cuts.Tokenize(";");

  const Int_t nDie=arrNames_prim->GetEntriesFast();
  for (int iCut = 0; iCut < nDie; ++iCut){
    TString cutDefinition(arrNames_prim->At(iCut)->GetName());
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(cutDefinition, isAOD);
    task->AddTrackCuts_primary_standard(filter);
    DoAdditionalWork(task);
  }

  // Adding PreFilter primary electron track cutsettings
  TObjArray*  arrNames_prim=names_Prim_Track_PreFilter_Cuts.Tokenize(";");
  const Int_t nDie=arrNames_prim->GetEntriesFast();

  for (int iCut = 0; iCut < nDie; ++iCut){
    TString cutDefinition(arrNames_prim->At(iCut)->GetName());
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(cutDefinition, isAOD);
    task->AddTrackCuts_primary_PreFilter(filter);
    DoAdditionalWork(task);
  }

  // #########################################################
  // #                 set pair cutsettings                  #
  // #########################################################
  // Adding primary pair cutsettings
  TObjArray*  arrNames_prim=names_Prim_Pair_Cuts.Tokenize(";");
  const Int_t nDie=arrNames_prim->GetEntriesFast();

  for (int iCut = 0; iCut < nDie; ++iCut){
    TString cutDefinition(arrNames_prim->At(iCut)->GetName());
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(cutDefinition, isAOD);
    task->AddPairCuts_primary(filter);
    DoAdditionalWork(task);
  }

  // Adding standard secondary pair cutsettings
  TObjArray*  arrNames_sec=names_Sec_Pair_standard_Cuts.Tokenize(";");
  const Int_t nDie=arrNames_sec->GetEntriesFast();

  for (int iCut = 0; iCut < nDie; ++iCut){
    TString cutDefinition(arrNames_sec->At(iCut)->GetName());
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(cutDefinition, isAOD);
    task->AddPairCuts_secondary_standard(filter);
    DoAdditionalWork(task);
  }

  // Adding PreFilter secondary pair cutsettings
  TObjArray*  arrNames_sec_PreFilter=names_Sec_Pair_PreFilter_Cuts.Tokenize(";");
  const Int_t nDie=arrNames_sec_PreFilter->GetEntriesFast();

  for (int iCut = 0; iCut < nDie; ++iCut){
    TString cutDefinition(arrNames_sec_PreFilter->At(iCut)->GetName());
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(cutDefinition, isAOD);
    task->AddPairCuts_secondary_PreFilter(filter);
    DoAdditionalWork(task);
  }



  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("efficiency%d", wagonnr), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;
}
