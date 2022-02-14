AliAnalysisTaskElectronEfficiencyV2* AddTask_slehner_ElectronEfficiency(
                                                                Double_t centMin=0.,
                                                                Double_t centMax=100.,
                                                                Bool_t PIDCorr=kFALSE,
                                                                Bool_t useAODFilterCuts=kFALSE,
                                                                TString TMVAweight,
                                                                Int_t genGroup=0,
                                                                Bool_t fromAlien,
                                                                TString date="ddmmyy",
                                                                Int_t wagonnr=0,
                                                                Bool_t purej=kTRUE
        ) {

  std::cout << "########################################\nADDTASK of ANALYSIS started\n########################################" << std::endl;
  TString configFile="Config_slehner_Efficiency.C";
  // #########################################################
  // #########################################################
  // Configuring Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  TString fileName        = AliAnalysisManager::GetCommonFileName();
  fileName = "AnalysisResults.root"; // create a subfolder in the file

  // #########################################################
  // #########################################################
  // Loading individual config file either local or from Alien

  
  if(fromAlien) TString configBasePath(TString("alien:///alice/cern.ch/user/s/selehner/configs/")+date+TString("/"));
  else TString configBasePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");

  TString configFile("Config_slehner_Efficiency.C");
  TString configFilePath(configBasePath+configFile);
  
  TString myConfig =TString::Format("alien_cp %s .",configFilePath.Data());
  gSystem->Exec(myConfig);
  
  gSystem->Exec(TString("alien_cp alien:///alice/cern.ch/user/s/selehner/cutlibs/LMEECutLib_slehner.C ."));

  configBasePath=Form("%s/",gSystem->pwd());

  TString configLMEECutLib("LMEECutLib_slehner.C");
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);

  Bool_t err = kFALSE;
  err |= gROOT->LoadMacro(configLMEECutLibPath.Data());
  err |= gROOT->LoadMacro(configFile.Data());
//  if (err) { Error("AddTask_slehner_ElectronEfficiency","Config(s) could not be loaded!"); }

  // #########################################################
  // #########################################################
  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(TString::Format("MaxTrCuts_MaxPIDCuts").Data());

  // #########################################################
  // #########################################################
  // Possibility to set generator. If nothing set all generators are taken into account
  // task->SetGeneratorName(generatorName);
  if(setGens){
    TString generators="";
    if(genGroup&1<<0) generators+= "Hijing_0;";
    if(genGroup&1<<1) generators+= "pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7;";
    if(genGroup&1<<2) generators+= "Pythia CC_8;Pythia BB_8;Pythia B_8";
    if(genGroup&1<<3) generators+= "Starlight_0;";
    if(genGroup&1<<4) generators+= "Hijing_1;";
    cout<<"Efficiency based on MC generators: "<<generators<<endl;
    TString generatorsPair=generators;
    task->SetGeneratorMCSignalName(generatorsPair);
    task->SetGeneratorULSSignalName(generators);
  }
  // #########################################################
  // #########################################################
	// Cut lib
  LMEECutLib* cutlib = new LMEECutLib();
  // Event selection. Is the same for all the different cutsettings
//  task->SetEnablePhysicsSelection(kTRUE);
//  task->SetTriggerMask(triggerNames);
  task->SelectCollisionCandidates(triggerNames);
  task->SetEventFilter(cutlib->GetEventCuts(centMin, centMax,purej)); // All cut sets have same event cuts
  
//  maybe redundant since already set in eventcuts above
  std::cout << "CentMin = " << centMin << "  CentMax = " << centMax << std::endl;
//  task->SetCentrality(centMin, centMax);
  
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
  
//  double mbinsarr[] = { 0.00, 0.02 ,0.04 ,0.08 ,0.14 ,0.22 ,0.38 ,0.54 ,1.1 ,1.7 ,2.5 ,2.9 ,3.0 ,3.1 ,3.3 ,3.5 ,4.0 ,5.0}; //Carsten's current
  vector<double>mbins;
  for(int i=0; i< sizeof(mbinsarr) / sizeof(mbinsarr[0]); i++){ mbins.push_back(mbinsarr[i]); }
  task->SetMassBins(mbins);

//  double ptbinsarr[]= {0.0,0.4,0.6,1,2.5,8};
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
//  AddSingleLegMCSignal(task);
  AddPairMCSignal(task);
  std::vector<Bool_t> DielectronsPairNotFromSameMother = AddSingleLegMCSignal(task);
  task->AddMCSignalsWhereDielectronPairNotFromSameMother(DielectronsPairNotFromSameMother);
  
  Config_slehner_Efficiency(task, useAODFilterCuts, TMVAweight);
  
//    Int_t trackCut = 0;
//    Int_t PIDCut = 0;
//    Int_t MVACut=0;
//    for(trackCut = 0; trackCut <=maxTrCuts; ++trackCut){
//        std::cout << "CutTr: "<<trackCut<<" CutPID: "<<PIDCut<<" MVA Cut: "<<-1+MVACut*0.2<<" added"<< std::endl;
//        AliAnalysisFilter* filter = SetupTrackCutsAndSettings(trackCut, PIDCut, MVACut, useAODFilterCuts,TMVAweight);
//        task->AddTrackCuts(filter);
//    }
//    trackCut = 0;
//    for(PIDCut = 0; PIDCut <=maxPIDCuts; ++PIDCut){  
//        std::cout << "CutTr: "<<trackCut<<" CutPID: "<<PIDCut<<" MVA Cut: "<<-1+MVACut*0.2<<" added"<< std::endl;
//        AliAnalysisFilter* filter = SetupTrackCutsAndSettings(trackCut, PIDCut, MVACut, useAODFilterCuts,TMVAweight);
//        task->AddTrackCuts(filter);
//    }
//    PIDCut = 0;
//    for(MVACut = 0; MVACut <=maxMVACut; ++MVACut){  
//        std::cout << "CutTr: "<<trackCut<<" CutPID: "<<PIDCut<<" MVA Cut: "<<-1+MVACut*0.2<<" added"<< std::endl;
//        AliAnalysisFilter* filter = SetupTrackCutsAndSettings(trackCut, PIDCut, MVACut, useAODFilterCuts,TMVAweight);
//        task->AddTrackCuts(filter);
//    }
    
  if(PIDCorr) setPIDCorrections(task);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("efficiency%d", wagonnr), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;
}
