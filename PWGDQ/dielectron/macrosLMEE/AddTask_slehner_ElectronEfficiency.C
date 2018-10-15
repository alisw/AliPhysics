//Names should contain a comma seperated list of cut settings
//Current options: all, electrons, kCutSet1, TTreeCuts, V0_TPCcorr, V0_ITScorr
AliAnalysisTaskElectronEfficiencyV2* AddTask_slehner_ElectronEfficiency(
                                                                Int_t trackCut=0,
                                                                Int_t PIDCut=0,
                                                                Int_t evCut=0,
                                                                Double_t centMin=0.,
                                                                Double_t centMax=100.,
                                                                Bool_t PIDCorr=kFALSE,
                                                                Bool_t useAODFilterCuts=kFALSE,
                                                                TString TMVAweight = "TMVAClassification_BDTG.weights_094.xml" 
        ) {

  std::cout << "########################################\nADDTASK of ANALYSIS started\n########################################" << std::endl;
  TString configFile="Config_slehner_Efficiency.C";
  Bool_t getFromAlien=kFALSE;
  // #########################################################
  // #########################################################
  // Configuring Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  TString fileName        = AliAnalysisManager::GetCommonFileName();
  fileName = "AnalysisResults.root"; // create a subfolder in the file

  // #########################################################
  // #########################################################
  // Loading individual config file either local or from Alien

  // TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
//  //Load updated macros from private ALIEN path
  
  if (getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/s/slehner/PWGDQ/dielectron/macrosLMEE/%s .",configFile.Data())))
      && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/s/slehner/PWGDQ/dielectron/macrosLMEE/LMEECutLib_slehner.C ."))
      ) {
    configBasePath=Form("%s/",gSystem->pwd());
  }
  TString configFilePath(configBasePath+configFile);
  TString configLMEECutLib("LMEECutLib_slehner.C");
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);

  Bool_t err = kFALSE;
  err |= gROOT->LoadMacro(configLMEECutLibPath.Data());
  err |= gROOT->LoadMacro(configFilePath.Data());
  // #########################################################
  // #########################################################
  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(TString::Format("MaxTrCuts_MaxPIDCuts").Data());

  // #########################################################
  // #########################################################
  // Possibility to set generator. If nothing set all generators are taken into account
  // task->SetGeneratorName(generatorName);
  task->SetGeneratorMCSignalName(generatorNameForMCSignal);
  task->SetGeneratorULSSignalName(generatorNameForULSSignal);

  // #########################################################
  // #########################################################
	// Cut lib
  LMEECutLib* cutlib = new LMEECutLib();
  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE);
  task->SetTriggerMask(triggerNames);
  task->SetEventFilter(cutlib->GetEventCuts(centMin, centMax)); // All cut sets have same event cuts
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  
//  maybe redundant since already set in eventcuts above
  std::cout << "CentMin = " << centMin << "  CentMax = " << centMax << std::endl;
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
  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  task->SetKinematicCuts(minGenPt, maxGenPt, minGenEta, maxGenEta);

  // #########################################################
  // #########################################################
  // Set Binning
//  if(usePtVector == kTRUE){
//    std::vector<Double_t> ptBinsVec;
//    for (UInt_t i = 0; i < nBinsPt+1; ++i){
//      ptBinsVec.push_back(ptBins[i]);
//    }
//    task->SetPtBins(ptBinsVec);
//  }
//  else
  task->SetPtBinsLinear   (minPtBin,  maxPtBin, stepsPtBin);
  task->SetEtaBinsLinear  (minEtaBin, maxEtaBin, stepsEtaBin);
  task->SetPhiBinsLinear  (minPhiBin, maxPhiBin, stepsPhiBin);
  task->SetThetaBinsLinear(minThetaBin, maxThetaBin, stepsThetaBin);
  task->SetMassBinsLinear (0, 5, 500);
  task->SetPairPtBinsLinear(minPairPtBin, maxPairPtBin, stepsPairPtBin);
  
//  double mbinsarr[] = { 0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.18, 0.22, 0.30, 0.38, 0.46, 0.62, 0.7, 0.86, 1.1, 1.70, 2.30, 2.70, 2.90, 3.00, 3.10, 3.30, 4.00, 5.00};
  double mbinsarr[] = { 0.00, 0.02 ,0.04 ,0.08 ,0.14 ,0.22 ,0.38 ,0.54 ,1.1 ,1.7 ,2.5 ,2.9 ,3.0 ,3.1 ,3.3 ,3.5 ,4.0 ,5.07, 0.86, 1.1, 1.70, 2.30, 2.70, 2.90, 3.00, 3.10, 3.30, 4.00, 5.00}; //Carsten's current
  vector<double>mbins;
  for(int i=0; i< sizeof(mbinsarr) / sizeof(mbinsarr[0]); i++){ mbins.push_back(mbinsarr[i]); }
  task->SetMassBins(mbins);
  
  double ptbinsarr[]= {0.0,0.4,0.6,1,2.5,8};
  vector<double>ptbins;
  for(int i=0; i< sizeof(ptbinsarr) / sizeof(ptbinsarr[0]); i++){ ptbins.push_back(ptbinsarr[i]);} 
 
  task->SetPairPtBins(ptbins);

  // #########################################################
  // Resolution File, If resoFilename = "" no correction is applied
  task->SetResolutionFile(resoFilename);
  task->SetResolutionFileFromAlien(resoFilenameFromAlien);
  task->SetSmearGenerated(SetGeneratedSmearingHistos);
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
  task->SetULSandLS(DoULSLS);
  task->SetDeactivateLS(DeactivateLS);

  // #########################################################
  // #########################################################
  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  AddSingleLegMCSignal(task);
  AddPairMCSignal(task);
  std::vector<Bool_t> DielectronsPairNotFromSameMother = AddSingleLegMCSignal(task);
  task->AddMCSignalsWhereDielectronPairNotFromSameMother(DielectronsPairNotFromSameMother);

  // #########################################################
  // Adding multiple cutsettings
  for(Int_t MVACut = 0; MVACut <= 10; ++MVACut){
    std::cout << "CutTr: "<<trackCut<<" CutPID: "<<PIDCut<<" MVA Cut: "<<MVACut*0.2<<" added"<< std::endl;
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(trackCut, PIDCut, MVACut, useAODFilterCuts,TMVAweight);
    task->AddTrackCuts(filter);
    }
  
 
// Adding multiple cutsettings
//  std::cout << "CutTr: "<<trackCut<<" CutPID: "<<PIDCut<<" being added"<< std::endl;
//  AliAnalysisFilter* filter = SetupTrackCutsAndSettings(trackCut, PIDCut, useAODFilterCuts);
//  task->AddTrackCuts(filter);
    
  if(PIDCorr) setPIDCorrections(task);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("efficiency%d", 0), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;
}
