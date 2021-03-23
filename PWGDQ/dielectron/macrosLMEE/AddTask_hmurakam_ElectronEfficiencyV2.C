AliAnalysisTaskElectronEfficiencyV2* AddTask_hmurakam_ElectronEfficiencyV2(TString name        = "test",
                                                                           Int_t whichGen      = 1, // 0=all sources, 1=HS, 2=Jpsi
                                                                           Bool_t isAOD        = kFALSE,
                                                                           Bool_t getFromAlien = kFALSE,
                                                                           TString configFile  = "./Config_hmurakam_ElectronEfficiencyV2.C",
                                                                           Bool_t tofcor       = kFALSE,
                                                                           TString year        = "16",
                                                                           Bool_t usePhiV      = kTRUE,
                                                                           Double_t maxMee     = 0.14,
                                                                           Double_t minphiv    = 2.0,
                                                                           Bool_t DeactivateLS = kFALSE)

{

  std::cout << "########################################\nADDTASK of ANALYSIS started\n########################################" << std::endl;

  // Configuring Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName = "LMEE.root"; // create a subfolder in the file

  // Loading individual config file either local or from Alien
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  //TString configBasePath= Form("%s/",gSystem->pwd());
  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/macrosLMEE/%s .",configFile.Data()))))
    {
      configBasePath=Form("%s/",gSystem->pwd());
    }
  TString configFilePath(configBasePath+configFile);
  
  // Loading config and cutlib
  Bool_t err=kFALSE;
  err |= gROOT->LoadMacro(configFilePath.Data());
  if (err) { Error("AddTask_hmurakam_ElectronEfficiencyV2","Config(s) could not be loaded!"); return 0x0; }

  std::cout << "Configpath:  " << configFilePath << std::endl;
 
  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(name.Data());

  //  Possibility to set generator. If nothing set all generators are taken into account
  if(whichGen == 0){
    std::cout << "No generator specified. Looking at all sources" << std::endl;
    task->SetGeneratorMCSignalName("");
    task->SetGeneratorULSSignalName("");
  } else if(whichGen == 1){
    if(year == "16"){
      std::cout << "2016 Sample Generator names specified -> Pythia CC_1, Pythia BB_1 and Pythia B_1" << std::endl;
      task->SetGeneratorMCSignalName("Pythia CC_1;Pythia BB_1;Pythia B_1");
      task->SetGeneratorULSSignalName("Pythia CC_1;Pythia BB_1;Pythia B_1");
    }else{
      std::cout << "2017 and 2018 Generator names specified -> Pythia CC_0, Pythia BB_0 and Pythia B_0" << std::endl;
      task->SetGeneratorMCSignalName("Pythia CC_0;Pythia BB_0;Pythia B_0");
      task->SetGeneratorULSSignalName("Pythia CC_0;Pythia BB_0;Pythia B_0");
    }
  } else if(whichGen == 2){
    //not used so far
    std::cout << "Generator names specified -> Jpsi2ee_1 and B2Jpsi2ee_1" << std::endl;
    task->SetGeneratorMCSignalName("Jpsi2ee_1;B2Jpsi2ee_1");
    task->SetGeneratorULSSignalName("Jpsi2ee_1;B2Jpsi2ee_1");
  }

  // Set TOF correction
  if(tofcor){
    SetTOFSigmaEleCorrection(task, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, year.Data());
  }

  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE);
  task->SetTriggerMask(triggerNames);
  task->SetEventFilter(SetupEventCuts(isAOD)); //returns eventCuts from Config.
  task->SetCentrality(centMin, centMax);
  
  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  // Do not set here your analysis pt-cuts
  task->SetMinPtGen(minGenPt);
  task->SetMaxPtGen(maxGenPt);
  task->SetMinEtaGen(minGenEta);
  task->SetMaxEtaGen(maxGenEta);

  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  task->SetKinematicCuts(minPtCut, maxPtCut, minEtaCut, maxEtaCut);

  // Set Binning
  // Pt
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
  // Mass
  if (useMassVector == true) {
    std::vector<double> massBinsVec;
    for (unsigned int i = 0; i < nBinsPt+1; ++i){
      massBinsVec.push_back(massBins[i]);
    }
    task->SetMassBins(massBinsVec);
  }
  else task->SetMassBinsLinear (minMassBin, maxMassBin, stepsMassBin);
  // PairPt
  if(usePairPtVector == true){
    std::vector<double> pairptBinsVec;
    for (unsigned int i = 0; i < nBinsPairPt+1; ++i){
      pairptBinsVec.push_back(pairptBins[i]);
    }
    task->SetPairPtBins(pairptBinsVec);
  }
  else task->SetPairPtBinsLinear(minPairPtBin, maxPairPtBin, stepsPairPtBin);

  // Resolution File, If resoFilename = "" no correction is applied
  SetResolutionFile(year);
  cout << resoFilename << endl;
  //  task->SetResolutionFile(resoFilename);
  //  task->SetResolutionFileFromAlien(resoFilenameFromAlien);
  task->SetResolutionFile(resoFilename,resoFilenameFromAlien);
  task->SetResolutionDeltaPtBinsLinear   (DeltaMomMin, DeltaMomMax, NbinsDeltaMom);
  task->SetResolutionRelPtBinsLinear   (RelMomMin, RelMomMax, NbinsRelMom);
  task->SetResolutionEtaBinsLinear  (DeltaEtaMin, DeltaEtaMax, NbinsDeltaEta);
  task->SetResolutionPhiBinsLinear  (DeltaPhiMin, DeltaPhiMax, NbinsDeltaPhi);
  task->SetResolutionThetaBinsLinear(DeltaThetaMin, DeltaThetaMax, NbinsDeltaTheta);

  // Set centrality correction. If resoFilename = "" no correction is applied
  //task->SetCentralityFile(centralityFilename);

  // Pairing related config
  task->SetDoPairing(DoPairing);
  task->SetULSandLS(DoULSLS);
  task->SetDeactivateLS(DeactivateLS);

  task->SetPhiVBinsLinear(minPhiVBin, maxPhiVBin, stepsPhiVBin);
  task->SetFillPhiV(kFALSE);
  task->SetPhiVCut(usePhiV,maxMee,minphiv);

  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  AddSingleLegMCSignal(task);
  AddPairMCSignal(task);

  // Adding cutsettings
  TObjArray*  arrNames=names.Tokenize(";");
  const Int_t nDie=arrNames->GetEntriesFast();

  printf("Add %d cuts\n",nDie);
  for (int iCut = 0; iCut < nDie; ++iCut){
    //TString cutDefinition(arrNames->At(iCut)->GetName());
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(iCut, isAOD);
    task->AddTrackCuts(filter);
  }

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("efficiency", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;
}
