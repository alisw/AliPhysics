
AliAnalysisTaskElectronEfficiencyV2* AddTask_jjung_efficiency(
		TString name = "name",
 		TString outputname = "AnalysisResults.root",	
		Bool_t isAOD=0, 
		Bool_t getFromAlien = kFALSE, 
		TString configFile="Config_jjung_lowmass.C",
                TString generatorNameForMCSignal  =  "Pythia",
                TString generatorNameForULSSignal =  "Pythia",
		Int_t triggerNames = AliVEvent::kINT7,
		Double_t centMin = 0.,
                Double_t centMax = 10.,
		Float_t PtMin =  0.2,
                Float_t PtMax = 10.0,
                Float_t EtaMin = -0.8,
                Float_t EtaMax = +0.8,
		Bool_t UsePtVec = kTRUE,
                Bool_t UsePairVec = kTRUE,
		Bool_t DoULSLS = kTRUE,
		std::string resoFilename   = "",
                Bool_t DoCocktailWeighting = kTRUE,
		std::string cocktailFilename = "",
		Bool_t DoCentralityCorrection = kFALSE,
		std::string centralityFilename = "",
		std::string calibFileName = "",
		Int_t wagonnr = 0) {

  std::cout << "########################################\nADDTASK of ANALYSIS started\n########################################" << std::endl;

  // #########################################################
  // #########################################################
  // Configuring Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    Error("AddTask_jjung_ElectronEfficiencyV2", "No analysis manager found.");
    return 0;
  }
  else { std::cout << "Analysis manager found!" << std::endl;}
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName = outputname; // create a subfolder in the file

  // #########################################################
  // #########################################################
  // Loading individual config file either local or from Alien

  // TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  TString configBasePath= "/data4/jung/localLegotrainNewEfficiency/";
  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/j/jjung/%s file:./",configFile.Data())))
      ) {
    configBasePath=Form("%s/",gSystem->pwd());
  }
  TString configFilePath(configBasePath+configFile);
  std::cout << configFilePath << std::endl;

  // Loading config and cutlib
  Bool_t err=kFALSE;
  err |= gROOT->LoadMacro(configFilePath.Data());
  //if (err) { Error("AddTask_jjung_ElectronEfficiency_v2","Config(s) could not be loaded!"); return 0x0; }

  // Download resolution file (configured in your config.C)
  //if (GetResolutionFromAlien == kTRUE){
  //  std::cout << "Trying to download resolution file" << std::endl;
  //  gSystem->Exec(Form("alien_cp alien://%s .",resoFilenameFromAlien.c_str()));
  //  std::cout << "Load resolution file from AliEn" << std::endl;
  //}
 
  //// Download centrality file (configured in your config.C)
  //if (GetCentralityFromAlien == kTRUE && !gSystem->Exec(Form("alien_cp alien://%s .",CentralityFilenameFromAlien.c_str()))){
  //  std::cout << "Load centrality file from AliEn" << std::endl;
  //}

  // #########################################################
  // #########################################################
  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(name.Data());
  gROOT->GetListOfSpecials()->Add(task);

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
  task->SetEventFilter((reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("GetEventCuts(%i)",wagonnr)))));
  if(wagonnr!=0) task->SetCentrality(centMin, centMax);
  // #########################################################
  // #########################################################
  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  // Do not set here your analysis pt-cuts
  task->SetMinPtGen(0.1);
  task->SetMaxPtGen(1e+10);
  task->SetMinEtaGen(-1.5);
  task->SetMaxEtaGen(+1.5);


  // #########################################################
  // #########################################################
  // 4D single efficiency from pairs
  //task->SetWriteLegsFromPair(WriteLegsFromPair);
  //task->SetPtMinLegsFromPair(ptMinLegsFromPair);
  //task->SetPtMaxLegsFromPair(ptMaxLegsFromPair);
  //task->SetPtNBinsLegsFromPair(ptNBinsLegsFromPair);
  //task->SetEtaMinLegsFromPair(etaMinLegsFromPair);
  //task->SetEtaMaxLegsFromPair(etaMaxLegsFromPair);
  //task->SetEtaNBinsLegsFromPair(etaNBinsLegsFromPair);
  //task->SetPhiMinLegsFromPair(phiMinLegsFromPair);
  //task->SetPhiMaxLegsFromPair(phiMaxLegsFromPair);
  //task->SetPhiNBinsLegsFromPair(phiNBinsLegsFromPair);
  //task->SetOpAngleMinLegsFromPair(opAngleMinLegsFromPair);
  //task->SetOpAngleMaxLegsFromPair(opAngleMaxLegsFromPair);
  //task->SetOpAngleNBinsLegsFromPair(opAngleNBinsLegsFromPair);


  // #########################################################
  // #########################################################
  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  task->SetKinematicCuts(PtMin, PtMax, EtaMin, EtaMax);

  // #########################################################
  // #########################################################
  // Set Binning
  if (UsePtVec == true) {

    const Int_t    nBinsPt = 73; 
    const Double_t PtBins[] = {
    0.00,0.050, 0.075, 0.080, 0.085, 0.090, 0.095, 0.10, 0.105, 0.11,
    0.115, 0.12, 0.125, 0.13,0.135, 0.14, 0.145, 0.15, 0.155, 0.16,
    0.165, 0.170, 0.175, 0.180, 0.185, 0.190, 0.195, 0.20,0.21,0.22,
    0.23,0.24,0.25, 0.26,0.27,0.28,0.29,0.30,0.32,0.34,
    0.36,0.38,0.40,0.42, 0.44,0.46,0.48, 0.49,0.52,0.54,
    0.56,0.58,0.60,0.65,0.70,0.75,0.80,0.90,1.00,1.10,
    1.20,1.40,1.60,1.80,2.00,2.40,2.80,3.20,3.70,4.50,
    6.00,8.00,10.0
    };
    std::vector<double> ptBinsVec;
    for (unsigned int i = 0; i < nBinsPt+1; ++i){
      ptBinsVec.push_back(PtBins[i]);
    }
    task->SetPtBins(ptBinsVec);
  }
  else task->SetPtBinsLinear   (0,  10, (Int_t)gROOT->ProcessLine("GetNbinsPt()"));
  task->SetEtaBinsLinear  (-1, 1, (Int_t)gROOT->ProcessLine("GetNbinsEta()"));
  task->SetPhiBinsLinear  (0, TMath::TwoPi(), (Int_t)gROOT->ProcessLine("GetNbinsPhi()"));
  task->SetThetaBinsLinear(0, TMath::TwoPi(), (Int_t)gROOT->ProcessLine("GetNbinsTheta()"));

  if (UsePairVec == true) {

    const Double_t MeeBins[] = { 0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
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
                               2.20,2.40,2.60,2.80,3.00,3.01,3.02,3.03,3.04,3.05,
                               3.06,3.07,3.08,3.09,3.10,3.30,3.50,4.00
                             };
    const Int_t nBinsMee = ( sizeof(MeeBins) / sizeof(MeeBins[0]) )-1;

    std::vector<double> meeBinsVec;
    for (unsigned int i = 0; i < nBinsMee+1; ++i){
      meeBinsVec.push_back(MeeBins[i]);
    }
    task->SetMassBins(meeBinsVec);
  }
  else task->SetMassBinsLinear (0., 4.0, (Int_t)gROOT->ProcessLine("GetNbinsMee()"));
  if (UsePairVec == true) {

    const Double_t PteeBins[] = {
                                0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
                                0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,
                                0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,
                                0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,
                                0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,
                                0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,
                                1.00,1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,
                                1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,
                                2.00,2.05,2.10,2.15,2.20,2.25,2.30,2.35,2.40,2.45,
                                2.50,2.60,2.70,2.80,2.90,3.00,3.10,3.20,3.30,3.40,
                                3.50,3.60,3.70,3.80,3.90,4.00,4.10,4.20,4.30,4.40,
                                4.50,5.00,5.50,6.00,6.50,7.00,8.00,10.00
                              };
    const Int_t nBinsPtee = ( sizeof(PteeBins) / sizeof(PteeBins[0]) )-1;    

    std::vector<double> pteeBinsVec;
    for (unsigned int i = 0; i < nBinsPtee+1; ++i){
      pteeBinsVec.push_back(PteeBins[i]);
    }
    task->SetPairPtBins(pteeBinsVec);
  }
  else task->SetPairPtBinsLinear(0., 10., (Int_t)gROOT->ProcessLine("GetNbinsPtee()"));

  // #########################################################
  // #########################################################
  // Resolution File, If resoFilename = "" no correction is applied
  task->SetResolutionFile(resoFilename);
  task->SetResolutionFileFromAlien("/alice/cern.ch/user/j/jjung/supportFiles/" + resoFilename);
  task->SetResolutionFile(resoFilename,"/alice/cern.ch/user/j/jjung/supportFiles/" + resoFilename);
  //task->SetSmearGenerated(SetGeneratedSmearingHistos); // cross check smearing the MC at single level and filling resolution maps
  
  task->SetResolutionDeltaPtBinsLinear   (-10., 2., (Int_t)gROOT->ProcessLine("GetNbinsDeltaMom()"));
  task->SetResolutionRelPtBinsLinear   (0., 2.,  (Int_t)gROOT->ProcessLine("GetNbinsRelMom()"));
  task->SetResolutionEtaBinsLinear  (-0.4, 0.4, (Int_t)gROOT->ProcessLine("GetNbinsDeltaEta()"));
  task->SetResolutionPhiBinsLinear  (-0.4, 0.4, (Int_t)gROOT->ProcessLine("GetNbinsDeltaPhi()"));
  task->SetResolutionThetaBinsLinear(-0.4, 0.4, (Int_t)gROOT->ProcessLine("GetNbinsDeltaTheta()"));


  // #########################################################
  // #########################################################
  // Set centrality correction. If resoFilename = "" no correction is applied
  if(wagonnr!=0){
    task->SetCentralityFile(centralityFilename);
    task->SetCentralityFile(centralityFilename,"/alice/cern.ch/user/j/jjung/supportFiles/" + centralityFilename);
  }
  // #########################################################
  // #########################################################
  // Set MCSignal and Cutsetting to fill the support histograms
  task->SetSupportHistoMCSignalAndCutsetting(0, 0);


  // #########################################################
  // #########################################################
  // Set Cocktail weighting
  task->SetDoCocktailWeighting(DoCocktailWeighting);
  task->SetCocktailWeighting(cocktailFilename);
  task->SetCocktailWeightingFromAlien("/alice/cern.ch/user/j/jjung/supportFiles/" + cocktailFilename);
  task->SetCocktailWeighting(cocktailFilename,"/alice/cern.ch/user/j/jjung/supportFiles/" + cocktailFilename);

  // #########################################################
  // #########################################################
  // Pairing related config
  Bool_t DoPairing = kTRUE;
  task->SetDoPairing(DoPairing);
  task->SetULSandLS(DoULSLS);
  //task->SetDeactivateLS(DeactivateLS);

  // #########################################################
  // #########################################################
  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  gROOT->ProcessLine(Form("AddSingleLegMCSignal(%s)",task->GetName()));//not task itself, task name
  gROOT->ProcessLine(Form("AddPairMCSignal(%s)"     ,task->GetName()));//not task itself, task name
  //Done in config
  //std::vector<bool> DielectronsPairNotFromSameMother = AddSingleLegMCSignal(task);
  //task->AddMCSignalsWhereDielectronPairNotFromSameMother(DielectronsPairNotFromSameMother);

  // #########################################################
  // #########################################################
  // Set mean and width correction for ITS, TPC and TOF
  //set PID map for ITS TOF in MC.
  TFile *rootfile = 0x0;
  std::string calibFileNameFromAlien = "/alice/cern.ch/user/j/jjung/supportFiles/" + calibFileName;

  if(calibFileName != "") rootfile = TFile::Open(calibFileName.c_str(),"READ");
  if(calibFileNameFromAlien != "" && !rootfile && getFromAlien){
    std::cout << "Location in AliEN: " << calibFileNameFromAlien << std::endl;
    gSystem->Exec(Form("alien_cp alien://%s file:./", calibFileNameFromAlien.c_str()));
    std::cout << "Copy resolution from Alien" << std::endl;
    rootfile = TFile::Open(calibFileName.c_str(), "READ");

    if (!rootfile) { 
      std::cout << "Could not open file: " << calibFileNameFromAlien << std::endl;
    }
  }
  if(rootfile && rootfile->IsOpen()){
    TH3D *h3mean_ITS  = (TH3D*)rootfile->Get("h3mean_ITS");
    TH3D *h3width_ITS = (TH3D*)rootfile->Get("h3width_ITS");
    h3mean_ITS ->SetDirectory(0);
    h3width_ITS->SetDirectory(0);
    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kITS, h3mean_ITS, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);
    task->SetWidthCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kITS, h3width_ITS, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta );

    TH3D *h3mean_TOF = (TH3D*)rootfile->Get("h3mean_TOF");
    TH3D *h3width_TOF = (TH3D*)rootfile->Get("h3width_TOF");
    h3mean_TOF ->SetDirectory(0);
    h3width_TOF->SetDirectory(0);
    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, h3mean_TOF, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);
    task->SetWidthCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, h3width_TOF, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);

    rootfile->Close();
  }



  // #########################################################
  // #########################################################
  // Adding cutsettings
  //TObjArray*  arrNames=names.Tokenize(";");
  //const Int_t nDie=arrNames->GetEntriesFast();
  const Int_t nDie = (Int_t)gROOT->ProcessLine("GetN()");
  for (int iCut = 0; iCut < nDie; ++iCut){
    //TString cutDefinition(arrNames->At(iCut)->GetName());
    //AliAnalysisFilter* filter = SetupTrackCutsAndSettings(iCut);
    AliAnalysisFilter *filter = reinterpret_cast<AliAnalysisFilter*>(gROOT->ProcessLine(Form("SetupTrackCutsAndSettings(%d)",iCut)));
    task->AddTrackCuts(filter);
  }


  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("efficiency%d",wagonnr), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;
}
