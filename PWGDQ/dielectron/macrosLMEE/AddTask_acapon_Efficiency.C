// Names should contain a comma seperated list of cut settings
// Current options: all, electrons, kCutSet1, TTreeCuts, V0_TPCcorr, V0_ITScorr
AliAnalysisTaskElectronEfficiencyV2* AddTask_acapon_Efficiency(TString names                = "kCutSet1",
                                                               Int_t whichGen               = 0, // 0=all sources, 1=Jpsi, 2=HS
                                                               Int_t wagonnr                = 0,
                                                               TString centrality           = "0;100;V0A", // min;max;estimator
                                                               Bool_t SDDstatus             = kTRUE,
                                                               Bool_t applyPIDcorr          = kTRUE,
                                                               std::string resoFilename     = "", // Leave blank to not use resolution files
                                                               Bool_t DoPairing             = kTRUE,
                                                               Bool_t DoULSLS               = kTRUE,
                                                               Bool_t DeactivateLS          = kTRUE,
                                                               Bool_t DoCocktailWeighting   = kFALSE,
                                                               Bool_t GetCocktailFromAlien  = kFALSE,
                                                               std::string CocktailFilename = "",
                                                               Bool_t cutlibPreloaded       = kFALSE,
                                                               Bool_t getFromAlien          = kFALSE)
{

  std::cout << "########################################\nADDTASK of ANALYSIS started\n########################################" << std::endl;

  TObjArray* arrNames = names.Tokenize(";");
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
  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(Form("%s%d",names.Data(), wagonnr));

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
  // Cut lib
  LMEECutLib* cutlib = new LMEECutLib(SDDstatus);
  // Event selection is the same for all cut settings
  task->SetEnablePhysicsSelection(kTRUE);
  task->SetTriggerMask(triggerNames);
  task->SetEventFilter(cutlib->GetEventCuts(kFALSE, kFALSE));

  // Set centrality requirements
  // There is certainly a better way to do this next section....
  TObjArray* centDetails = centrality.Tokenize(";");
  TString tempMinFloat = TString(centDetails->At(0)->GetName());
  TString tempMaxFloat = TString(centDetails->At(1)->GetName());
  Float_t minCent = tempMinFloat.Atoi();
  Float_t maxCent = tempMaxFloat.Atoi();
  std::cout << "CentMin = " <<  minCent << "  CentMax = " <<  maxCent << std::endl;
  task->SetCentrality(minCent, maxCent);
  TString whichCentEst = TString(centDetails->At(2)->GetName());
  task->SetCentralityEstimator(whichCentEst);

  // #########################################################
  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  // Do not set here your analysis pt-cuts
  task->SetMinPtGen(minGenPt);
  task->SetMaxPtGen(maxGenPt);
  task->SetMinEtaGen(minGenEta);
  task->SetMaxEtaGen(maxGenEta);


  // #########################################################
  // Set kinematic cuts for pairing
  task->SetKinematicCuts(minPtCut, maxPtCut, minEtaCut, maxEtaCut);

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

  // Mass and pairPt binning should match kMee and kPtee in Config_acapon.C
  Double_t massBinsArr[] = {0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
                            0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,
                            0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,
                            0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,
                            0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,
                            0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,
                            0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,
                            0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,
                            0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,
                            0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,
                            1.00,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,
                            1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,
                            2.10,2.20,2.30,2.40,2.50,2.60,2.70,2.80,2.90,3.00,
                            3.01,3.02,3.03,3.04,3.05,3.06,3.07,3.08,3.09,3.10,
                            3.11,3.12,3.30,3.50,4.00,4.50,5.00};
  std::vector<Double_t> massBins;
  for(Int_t j = 0; j < sizeof(massBinsArr)/sizeof(massBinsArr[0]); j++){
    massBins.push_back(massBinsArr[j]);
  }
  task->SetMassBins(massBins);

  Double_t pairPtBinsArr[] = {0.000,0.025,0.050,0.075,0.100,0.125,0.150,0.175,0.200,0.225,0.250,
                              0.275,0.300,0.325,0.350,0.375,0.400,0.425,0.450,0.475,0.500,0.550,
                              0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,1.000,1.050,1.100,
                              1.150,1.200,1.250,1.300,1.350,1.400,1.450,1.500,1.550,1.600,1.650,
                              1.700,1.750,1.800,1.850,1.900,1.950,2.000,2.050,2.100,2.150,2.200,
                              2.250,2.300,2.350,2.400,2.450,2.500,2.600,2.700,2.800,2.900,3.000,
                              3.100,3.200,3.300,3.400,3.500,3.600,3.700,3.800,3.900,4.000,4.100,
                              4.200,4.300,4.400,4.500,5.000,5.500,6.000,6.500,7.000,8.000,10.00};
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
  // Adding cut settings
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
