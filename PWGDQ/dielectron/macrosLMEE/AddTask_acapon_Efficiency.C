//Names should contain a comma seperated list of cut settings
//Current options: all, electrons, kCutSet1, TTreeCuts, V0_TPCcorr, V0_ITScorr
AliAnalysisTaskElectronEfficiencyV2* AddTask_acapon_Efficiency(TString names = "kCutSet1",
                                                               Int_t whichGen = 0, // 0=gen. purpose, 1=Jpsi, 2=HF, 3=all
                                                               Int_t wagonnr = 0,
                                                               Int_t centrality = 0,
                                                               Bool_t cutlibPreloaded = kFALSE,
                                                               Bool_t getFromAlien = kFALSE
																																) {

  std::cout << "########################################\nADDTASK of ANALYSIS started\n########################################" << std::endl;
	
  TString configFile  = "Config_acapon_Efficiency.C";
	TObjArray *arrNames = names.Tokenize(";");
	Int_t nDie          = arrNames->GetEntries();
	Printf("Number of implemented cuts: %i", nDie);

	Bool_t SDDstatus = kTRUE;

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
  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/a/acapon/PWGDQ/dielectron/macrosLMEE/%s .",configFile.Data())))
      && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PWGDQ/dielectron/macrosLMEE/LMEECutLib_acapon.C ."))
      ) {
    configBasePath=Form("%s/",gSystem->pwd());
  }
  TString configFilePath(configBasePath+configFile);
  TString configLMEECutLib("LMEECutLib_acapon.C");
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);


  Bool_t err = kFALSE;
  if(!cutlibPreloaded){ // should not be needed but seems to be...
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
  task->SetEventFilter(cutlib->GetEventCuts(LMEECutLib::kAllSpecies)); // All cut sets have same event cuts

  Double_t centMin = 0.;
  Double_t centMax = 90.;
  GetCentrality(centrality, centMin, centMax);
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
  task->SetKinematicCuts(minPtCut, maxPtCut, minEtaCut, maxEtaCut);

  // #########################################################
  // #########################################################
  // Set Binning
  if(usePtVector == kTRUE){
    std::vector<Double_t> ptBinsVec;
    for (UInt_t i = 0; i < nBinsPt+1; ++i){
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
  // Adding cutsettings
  for(Int_t iCut = 0; iCut < nDie; ++iCut){
    TString cutDefinition(arrNames->At(iCut)->GetName());
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(cutDefinition);
    task->AddTrackCuts(filter);
		Printf("Successfully added task with cut set: %s\n", cutDefinition);
    //DoAdditionalWork(task);
  }


  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("efficiency%d", wagonnr), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;
}
