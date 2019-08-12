// Names should contain a comma seperated list of cut settings
// Current options: all, electrons, kCutSet1, TTreeCuts, V0_TPCcorr, V0_ITScorr
AliAnalysisTaskElectronEfficiencyV2* AddTask_acapon_Efficiency(TString names                = "kCutSet1",
                                                               Int_t whichGen               = 0, // 0=all sources, 1=Jpsi, 2=HS
                                                               Int_t wagonnr                = 0,
                                                               Int_t centrality             = 0,
                                                               Bool_t SDDstatus             = kTRUE,
                                                               Bool_t applyPIDcorr          = kTRUE,
                                                               std::string resoFilename     = "", // Leave blank to not use resolution files
                                                               Bool_t DoPairing             = kTRUE,
                                                               Bool_t DoULSLS               = kTRUE,
                                                               Bool_t DeactivateLS          = kTRUE,
                                                               Bool_t DoCocktailWeighting   = kFALSE,
                                                               Bool_t GetCocktailFromAlien  = kFALSE,
                                                               std::string CocktailFilename = "",
                                                               Bool_t useRun1binning        = kFALSE,
                                                               Bool_t cutlibPreloaded       = kFALSE,
                                                               Bool_t getFromAlien          = kFALSE)
{

  std::cout << "########################################\nADDTASK of ANALYSIS started\n########################################" << std::endl;

  TObjArray *arrNames = names.Tokenize(";");
  Int_t nDie          = arrNames->GetEntries();
  Printf("Number of implemented cuts: %i", nDie);

  std::string resoFilenameFromAlien = "/alice/cern.ch/user/a/acapon/ResolutionFiles/";
  resoFilenameFromAlien.append(resoFilename);

  std::string CocktailFilenameFromAlien = "/alice/cern.ch/user/a/acapon/Cocktails/";
  CocktailFilenameFromAlien.append(CocktailFilename);

  // #########################################################
  // Configuring Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  TString fileName        = AliAnalysisManager::GetCommonFileName();
  fileName = "AnalysisResults.root"; // create a subfolder in the file

  // #########################################################
  // Loading individual config file either local or from Alien
  TString configBasePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");
  TString configFile("Config_acapon_Efficiency.C");
  TString configLMEECutLib("LMEECutLib_acapon.C");
  // Load updated macros from private ALIEN path
  if(getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/a/acapon/PWGDQ/dielectron/macrosLMEE/%s .",configFile.Data())))
      && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PWGDQ/dielectron/macrosLMEE/LMEECutLib_acapon.C ."))
      ){
    configBasePath=Form("%s/",gSystem->pwd());
  }
  TString configFilePath(configBasePath+configFile);
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);

  // Load dielectron configuration files
  if(!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data())){
    gROOT->LoadMacro(configLMEECutLibPath.Data());
  }
  if(!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data())){
    gROOT->LoadMacro(configFilePath.Data());
  }

  // #########################################################
  // #########################################################
  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(Form("%s%d",names.Data(), wagonnr));

  // #########################################################
  // #########################################################
  // Possibility to set generator. If nothing set all generators are taken into account
  if(whichGen == 0){
    std::cout << "No generator specified. Looking at all sources" << std::endl;
    task->SetGeneratorMCSignalName("");
    task->SetGeneratorULSSignalName("");
  }
  else if(whichGen == 1){
    std::cout << "Generator names specified -> Jpsi2ee_1 and B2Jpsi2ee_1" << std::endl;
    task->SetGeneratorMCSignalName("Jpsi2ee_1;B2Jpsi2ee_1");
    task->SetGeneratorULSSignalName("Jpsi2ee_1;B2Jpsi2ee_1");
  }
  else if(whichGen == 2){
    std::cout << "Generator names specified -> Pythia CC_1, Pythia BB_1 and Pythia B_1" << std::endl;
    task->SetGeneratorMCSignalName("Pythia CC_1;Pythia BB_1;Pythia B_1");
    task->SetGeneratorULSSignalName("Pythia CC_1;Pythia BB_1;Pythia B_1");
  }

  // #########################################################
  // #########################################################
  // Cut lib
  LMEECutLib* cutlib = new LMEECutLib(SDDstatus);
  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE);
  task->SetTriggerMask(triggerNames);
  task->SetEventFilter(cutlib->GetEventCuts()); // All cut sets have same event cuts

  Double_t centMin = -99.;
  Double_t centMax = -90.;
  GetCentrality(centrality, centMin, centMax);
  std::cout << "CentMin = " << centMin << "  CentMax = " << centMax << std::endl;
  task->SetCentrality(centMin, centMax);
  if(centrality == 8){
    task->SetRun1Analysis(kTRUE);
  }

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
  // Set kinematic cuts for pairing
  task->SetKinematicCuts(minPtCut, maxPtCut, minEtaCut, maxEtaCut);

  // #########################################################
  // #########################################################
  // Set Binning
  if(usePtVector == kTRUE){
    std::vector<Double_t> ptBinsVec;
    for(UInt_t i = 0; i < nBinsPt+1; ++i){
      ptBinsVec.push_back(ptBins[i]);
    }
    task->SetPtBins(ptBinsVec);
  }
  else{
    task->SetPtBinsLinear(minPtBin,  maxPtBin, stepsPtBin);
  }
  task->SetEtaBinsLinear  (minEtaBin, maxEtaBin, stepsEtaBin);
  task->SetPhiBinsLinear  (minPhiBin, maxPhiBin, stepsPhiBin);
  task->SetThetaBinsLinear(minThetaBin, maxThetaBin, stepsThetaBin);
  /* task->SetMassBinsLinear (minMassBin, maxMassBin, stepsMassBin); */
  /* task->SetPairPtBinsLinear(minPairPtBin, maxPairPtBin, stepsPairPtBin); */

  // Use non linear binning for mass and pair pt
  Double_t massBinsArr[] = {0.00, 0.02, 0.04, 0.10, 0.14, 0.18, 0.24, 0.28, 0.34, 0.38,
                            0.44, 0.50, 0.60, 0.70, 0.76, 0.80, 0.86, 0.90, 0.96, 1.00,
                            1.04, 1.10, 1.40, 1.70, 2.00, 2.40, 2.70, 2.80, 2.90, 3.00,
                            3.10, 3.30, 3.50, 4.00, 5.00};
  // Mass bins from Run 1 (Theo's analysis)
  Double_t massBinsRun1arr[] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.47, 0.62, 0.7,
                                0.77, 0.8, 0.9, 0.95, 0.99, 1.02, 1.03, 1.1, 1.4, 1.7,
                                2, 2.3, 2.6, 2.8, 2.9, 3, 3.04, 3.08, 3.1, 3.12, 3.2,
                                3.5, 5.00};
  std::vector<Double_t> massBins;
  if(!useRun1binning){
    for(Int_t j = 0; j < sizeof(massBinsArr)/sizeof(massBinsArr[0]); j++){
      massBins.push_back(massBinsArr[j]);
    }
  }else{
    for(Int_t j = 0; j < sizeof(massBinsRun1arr)/sizeof(massBinsRun1arr[0]); j++){
      massBins.push_back(massBinsRun1arr[j]);
    }
  }
  task->SetMassBins(massBins);


  Double_t pairPtBinsArr[] = {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35,
                              0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,
                              0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30,
                              1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10,
                              2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90,
                              3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70,
                              3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50,
                              5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 10.0};
  std::vector<Double_t> pairPtBins;
  for(Int_t k = 0; k < sizeof(pairPtBinsArr)/sizeof(pairPtBinsArr[0]); k++){
    pairPtBins.push_back(pairPtBinsArr[k]);
  }
  task->SetPairPtBins(pairPtBins);

  // #########################################################
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
  // Add MCSignals
  // Add single track signal definition, and return vector used for calculating
  // efficiencies from non resonant dielectron pairs. E.g use all ULS pairs,
  // electrons from D or B decays etc
  std::vector<Bool_t> DielectronsPairNotFromSameMother = AddSingleLegMCSignal(task);
  task->AddMCSignalsWhereDielectronPairNotFromSameMother(DielectronsPairNotFromSameMother);
  // Add standard "same-mother" dielectron pair signal
  AddPairMCSignal(task);

  // #########################################################
  // Adding cutsettings
  for(Int_t iCut = 0; iCut < nDie; ++iCut){
    TString cutDefinition(arrNames->At(iCut)->GetName());
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(cutDefinition, SDDstatus);
    if(!filter){
      std::cout << "Invalid cut setting specified!!" << std::endl;
      return 0x0;
    }
    task->AddTrackCuts(filter);
    Printf("Successfully added task with cut set: %s\n", cutDefinition);
    // Apply PID post calibration to ITS(0) and TOF(1)
    if(applyPIDcorr){
      ApplyPIDpostCalibration(task, 0, SDDstatus);
      ApplyPIDpostCalibration(task, 1, SDDstatus);
    }
  }

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("efficiency%d", wagonnr), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;
}
