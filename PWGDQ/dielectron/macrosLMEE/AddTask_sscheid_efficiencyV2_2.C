AliAnalysisTaskElectronEfficiencyV2* AddTask_sscheid_efficiencyV2_2(Bool_t getFromAlien = kFALSE,
                                                                                TString cFileName = "Config_sscheid_efficiencyV2_raph.C",
                                                                                UInt_t trigger = AliVEvent::kINT7,
                                                                                Int_t rejpileup = 1,
                                                                                const Int_t CenMin = 0,
                                                                                const Int_t CenMax = 10,
                                                                                const Float_t PtMin = 0.2,
                                                                                const Float_t PtMax = 10.0,
                                                                                const Float_t EtaMin = -0.8,
                                                                                const Float_t EtaMax = +0.8,
                                                                                const Bool_t UsePtVec = kTRUE,
                                                                                const Bool_t DoULSLS = kTRUE,
                                                                                const TString generators = "",
                                                                                const std::string resolutionFilename = "",
                                                                                const std::string cocktailFilename = "",
                                                                                const std::string centralityFilename = "",
                                                                                TString calibFileName = "",
                                                                                const TString outname = "AnalysisResults.root",
                                                                                Int_t version = 0,
                                                                                TString CentralityEstimator = "")
{

  //get the current analysis manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_rbailhac_test", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath = "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  if (!gSystem->AccessPathName(cFileName)) {
    printf("Configfile already present\n");
    configBasePath = Form("%s/", gSystem->pwd());
  } else if (getFromAlien && (!gSystem->Exec(Form("alien_cp alien:/alice/cern.ch/user/h/hscheid/PWGDQ/dielectron/macrosLMEE/%s file:.", cFileName.Data())))) {
    printf("Copy Configfile from alien\n");
    configBasePath = Form("%s/", gSystem->pwd());
  }

  TString configFilePath(configBasePath + cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;

  if (!gROOT->GetListOfGlobalFunctions()->FindObject("SetupTrackCutsAndSettings")) {
    printf("Load macro now\n");
    gROOT->LoadMacro(configFilePath.Data());
  }

  // trigger
  TString triggername = "NULL";
  if (trigger == (UInt_t)AliVEvent::kINT7)
    triggername = "kINT7";
  else if (trigger == (UInt_t)AliVEvent::kCentral)
    triggername = "kCentral";
  else if (trigger == (UInt_t)AliVEvent::kSemiCentral)
    triggername = "kSemiCentral";
  else if (trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral))
    triggername = "kCombinedCentralityTriggers";
  else if (trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral))
    triggername = "kCombinedCentral";
  else if (trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kSemiCentral))
    triggername = "kCombinedSemiCentral";

  // generators
  TString suffixgen = "";
  if (generators.Contains("Pythia CC") && (generators.Contains("Pythia BB") || generators.Contains("Pythia B")))
    suffixgen = "_CC_BB";
  else if (generators.Contains("Pythia CC"))
    suffixgen = "_CC";
  else if (generators.Contains("Pythia BB") || generators.Contains("Pythia B"))
    suffixgen = "_BB";
  else if (generators.Contains("pizero_0"))
    suffixgen = "_LF";
  else
    suffixgen = "";

  //create task and add it to the manager (MB)
  TString appendix;
  appendix += TString::Format("Cen%d_%d_%s_%s_Pileup%d_%d", CenMin, CenMax, triggername.Data(), suffixgen.Data(), rejpileup, version);
  if (CentralityEstimator.Contains("V0")) appendix += TString::Format("_Cen%s", CentralityEstimator.Data());
  printf("appendix %s\n", appendix.Data());

  //##########################################################
  //############################################################
  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(Form("TaskElectronEfficiencyV2_%s", appendix.Data()));
  gROOT->GetListOfSpecials()->Add(task); //this is only for ProcessLine(AddMCSignal);

  // #########################################################
  // #########################################################
  // Set TOF correction
  // if(tofcor){
  //  SetEtaCorrectionTOFRMS(task, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
  //  SetEtaCorrectionTOFMean(task, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
  //}

  // #########################################################
  // #########################################################
  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE); //always ON in Run2 analyses for both data and MC.
  task->SetTriggerMask(trigger);
  if (CentralityEstimator.Contains("V0")) {
    printf("Use %s centrality estimator for event selection\n", CentralityEstimator.Data());
    task->SetEventFilter((reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("GetEventCuts(%f,%f,%d,\"%s\")", (Float_t)CenMin, (Float_t)CenMax, rejpileup, CentralityEstimator.Data())))));
  } else {
    printf("Use default old centrality estimator V0M for event selection\n");
    task->SetEventFilter((reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("GetEventCuts(%f,%f,%d,\"%s\")", (Float_t)CenMin, (Float_t)CenMax, rejpileup, "V0M")))));
  }

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
  // Set minimum and maximum values for pairing
  task->SetKinematicCuts(PtMin, PtMax, EtaMin, EtaMax);

  // #########################################################
  // #########################################################
  // Set Binning single variables
  if (UsePtVec == true) {
    const Int_t Npt = 68;
    Double_t pte[Npt] = {0.00, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.20, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, 0.295, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00, 2.40, 2.80, 3.20, 3.70, 4.50, 6.00, 8.00, 12.0};
    std::vector<double> v_pte(pte, std::end(pte));
    task->SetPtBins(v_pte);
  } else
    task->SetPtBinsLinear(0, 10, 100);
  task->SetEtaBinsLinear(-1., 1., 20);           // 40 before
  task->SetPhiBinsLinear(0, TMath::TwoPi(), 36); // 90 before
  task->SetThetaBinsLinear(0, TMath::TwoPi(), 60);

  // pair variables
  const Int_t Nmee = 150;
  Double_t mee[Nmee] = {};
  for (Int_t i = 0; i < 110; i++)
    mee[i] = 0.01 * (i - 0) + 0.0; //from 0 to 1.09 GeV/c2, every 0.01 GeV/c2
  for (Int_t i = 110; i < Nmee; i++)
    mee[i] = 0.1 * (i - 110) + 1.1; //from 1.1 to 5 GeV/c2, evety 0.1 GeV/c2
  std::vector<double> v_mee(mee, std::end(mee));

  const Int_t NpTee = 130;
  Double_t pTee[NpTee] = {};
  for (Int_t i = 0; i < 20; i++)
    pTee[i] = 0.005 * (i - 0) + 0.0; //from 0 to 0.095 GeV/c, every 0.005 GeV/c
  for (Int_t i = 20; i < 119; i++)
    pTee[i] = 0.1 * (i - 20) + 0.1; //from 0.1 to 9.9 GeV/c, evety 0.1 GeV/c
  for (Int_t i = 119; i < NpTee; i++)
    pTee[i] = 1.0 * (i - 119) + 10.0; //from 10 to 20 GeV/c, evety 1.0 GeV/c
  std::vector<double> v_pTee(pTee, std::end(pTee));
  task->SetMassBins(v_mee);
  task->SetPairPtBins(v_pTee);
  //task->SetPhiVBinsLinear(0, TMath::Pi(), 100);
  //task->SetFillPhiV(kTRUE);

  // #########################################################
  // #########################################################
  //task->SetSmearGenerated(kFALSE); // cross check smearing the MC at single level and filling resolution maps
  // Resolution File, If resoFilename = "" no correction is applied
  task->SetResolutionFile(resolutionFilename, "/alice/cern.ch/user/h/hscheid/supportFiles/" + resolutionFilename);
  task->SetResolutionDeltaPtBinsLinear(-10., 2., (Int_t)gROOT->ProcessLine("GetNbinsDeltaMom()"));
  task->SetResolutionRelPtBinsLinear(0., 2., (Int_t)gROOT->ProcessLine("GetNbinsRelMom()"));
  task->SetResolutionEtaBinsLinear(-0.4, 0.4, (Int_t)gROOT->ProcessLine("GetNbinsDeltaEta()"));
  task->SetResolutionPhiBinsLinear(-0.4, 0.4, (Int_t)gROOT->ProcessLine("GetNbinsDeltaPhi()"));
  task->SetResolutionThetaBinsLinear(-0.4, 0.4, (Int_t)gROOT->ProcessLine("GetNbinsDeltaTheta()"));

  //###########################################################
  //############################################################
  // Set MCSignal and Cutsetting to fill the support histograms
  task->SetSupportHistoMCSignalAndCutsetting(0, 0); //fill support histograms for first MCsignal and first cutsetting

  // #########################################################
  // #########################################################
  // Set centrality correction. If resoFilename = "" no correction is applied
  task->SetCentralityFile(centralityFilename, "/alice/cern.ch/user/h/hscheid/supportFiles/" + centralityFilename);
  if (CentralityEstimator.Contains("V0")) {
    printf("Set %s centrality estimator \n", CentralityEstimator.Data());
    task->SetCentralityEstimator(CentralityEstimator.Data());
  } else {
    printf("Use default old centrality estimator V0M for the task\n");
  }

  // #########################################################
  // #########################################################
  // Pairing related config
  task->SetDoPairing(kTRUE);
  task->SetULSandLS(DoULSLS);

  //#####################################################
  //######################################################
  // Generators
  cout << "Efficiency based on MC generators: " << generators << endl;
  TString generatorsPair = generators;
  task->SetGeneratorMCSignalName(generatorsPair);
  task->SetGeneratorULSSignalName(generators);

  //###############################################
  //##############################################
  if (cocktailFilename != "") task->SetDoCocktailWeighting(kTRUE);
  task->SetCocktailWeighting(cocktailFilename, "/alice/cern.ch/user/h/hscheid/supportFiles/" + cocktailFilename);

  // #########################################################
  // #########################################################
  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  gROOT->ProcessLine(Form("AddSingleLegMCSignal(%s)", task->GetName())); //not task itself, task name
  gROOT->ProcessLine(Form("AddPairMCSignal(%s)", task->GetName()));      //not task itself, task name

  // #########################################################
  // #########################################################
  // Adding cutsettings
  // Cuts
  // Number of cuts
  const Int_t nDie = (Int_t)gROOT->ProcessLine("GetN()");
  //add dielectron analysis with different cuts to the task
  for (Int_t i = 0; i < nDie; ++i) { //nDie defined in config file
    AliAnalysisFilter* filter = reinterpret_cast<AliAnalysisFilter*>(gROOT->ProcessLine(Form("SetupTrackCutsAndSettings(%d)", i)));
    task->AddTrackCuts(filter);
  }

  //########################################
  //########################################
  TString outlistname = Form("efficiency_%s", appendix.Data());
  const TString fileName = outname;
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(outlistname, TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", fileName.Data())));

  return task;
}
