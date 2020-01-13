AliAnalysisTaskElectronEfficiencyV2* AddTask_dsekihat_ElectronEfficiencyV2_PbPb(
    Bool_t getFromAlien = kFALSE,
    TString configFile="Config_dsekihat_ElectronEfficiencyV2_PbPb.C",
    UInt_t trigger = AliVEvent::kINT7,
    const Int_t CenMin =  0,
    const Int_t CenMax = 10,
    const Float_t PtMin =  0.2,
    const Float_t PtMax = 10.0,
    const Float_t EtaMin = -0.8,
    const Float_t EtaMax = +0.8,
    const TString generators = "pizero_0;eta_1;etaprime_2;rho_3;omega_4;phi_5;jpsi_6;Pythia CC_0;Pythia BB_0;Pythia B_0;",
    const Bool_t isLHC19f2 = kTRUE,
    const std::string resolutionAlien ="",
    const std::string cocktailAlien   ="",
    const std::string centralityAlien =""
    ){

  // Configuring Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  Bool_t isAOD = mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();
  printf("isAOD = %d\n",isAOD);

  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(Form("TaskElectronEfficiencyV2_Cen%d_%d_kINT7",CenMin,CenMax));
  gROOT->GetListOfSpecials()->Add(task);//this is only for ProcessLine(AddMCSignal);

  TString configBasePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");
  //TString configBasePath("./");
  if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/macrosLMEE/%s .",configFile.Data())))){
    configBasePath=Form("%s/",gSystem->pwd());
  }
  TString configFilePath(configBasePath + configFile);

  //load cut library first
  TString libFilePath(configBasePath + "LMEECutLib_dsekihat.C");
  std::cout << "Configpath:  " << configFilePath << std::endl;

  gROOT->LoadMacro(libFilePath.Data());//library first
  gROOT->LoadMacro(configFilePath.Data());

  // Adding cutsettings
  Int_t nCut = gROOT->ProcessLine("GetN()");
  for (int iCut = 0; iCut < nCut; ++iCut){
    AliAnalysisFilter *filter = reinterpret_cast<AliAnalysisFilter*>(gROOT->ProcessLine(Form("Config_dsekihat_ElectronEfficiencyV2_PbPb(%d,%d,%f,%f,%f,%f)",iCut,isAOD,PtMin,PtMax,EtaMin,EtaMax)));
    task->AddTrackCuts(filter);
  }

  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE);//always ON in Run2 analyses for both data and MC.
  task->SetTriggerMask(trigger);
  task->SetEventFilter(reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("LMEECutLib::SetupEventCuts(%f,%f,%d,\"%s\")",(Float_t)CenMin,(Float_t)CenMax,kTRUE,"V0M"))));//kTRUE is for Run2
  //task->SetCentralityEstimator("V0M");
  //task->SetCentrality(CenMin, CenMax);

  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  // Do not set here your analysis pt-cuts
  task->SetMinPtGen(0.1);
  task->SetMaxPtGen(1e+10);
  task->SetMinEtaGen(-1.5);
  task->SetMaxEtaGen(+1.5);

  // Set minimum and maximum values for pairing
  task->SetKinematicCuts(PtMin, PtMax, EtaMin, EtaMax);

  // Set Binning
  task->SetPtBinsLinear    (0, 10, 1000);
  task->SetEtaBinsLinear   (-1, +1, 20);
  task->SetPhiBinsLinear   (0, TMath::TwoPi(), 90);
  task->SetThetaBinsLinear (0, TMath::TwoPi(), 60);
  task->SetMassBinsLinear  (0, 5, 500);
  task->SetPairPtBinsLinear(0, 20, 200);
  task->SetPhiVBinsLinear  (0, TMath::Pi(), 72);
  task->SetFillPhiV(kTRUE);

  task->SetSmearGenerated(kFALSE);
  task->SetResolutionDeltaPtBinsLinear( -10,  +10, 2000);
  task->SetResolutionRelPtBinsLinear  (   0,    2,  400);
  task->SetResolutionEtaBinsLinear    (-0.4, +0.4,  200);
  task->SetResolutionPhiBinsLinear    (-0.4, +0.4,  200);
  task->SetResolutionThetaBinsLinear  (-0.4, +0.4,  200);

  // Set MCSignal and Cutsetting to fill the support histograms
  task->SetSupportHistoMCSignalAndCutsetting(0,0);//fill support histograms for first MCsignal and first cutsetting

  // Pairing related config
  task->SetDoPairing(kTRUE);
  //task->SetULSandLS(kTRUE);
  //task->SetDeactivateLS(kTRUE);

  //TString generators = "pizero_0;eta_1;etaprime_2;rho_3;omega_4;phi_5;jpsi_6;";
  //TString generators = "pizero_0;eta_1;etaprime_2;rho_3;omega_4;phi_5;jpsi_6;Pythia CC_0;Pythia BB_0;Pythia B_0;";
  //TString generators = "Pythia CC_0;Pythia BB_0;Pythia B_0;";
  //TString generators = "Pythia CC_0;";
  //TString generators = "Pythia BB_0;Pythia B_0;";
  //TString generators = "Pythia BB_0;";
  //TString generators = "Pythia B_0";

  cout<<"Efficiency based on MC generators: " << generators <<endl;
  TString generatorsPair=generators;
  task->SetGeneratorMCSignalName(generatorsPair);
  task->SetGeneratorULSSignalName(generators);
  task->SetLHC19f2MC(isLHC19f2);

  // Resolution File, If resoFilename = "" no correction is applied
  task->SetCentralityFile(centralityAlien);
  task->SetResolutionFileFromAlien(resolutionAlien);
  task->SetCocktailWeightingFromAlien(cocktailAlien);

  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  gROOT->ProcessLine(Form("AddSingleLegMCSignal(%s)",task->GetName()));//not task itself, task name
  gROOT->ProcessLine(Form("AddPairMCSignal(%s)"     ,task->GetName()));//not task itself, task name

  const TString fileName = AliAnalysisManager::GetCommonFileName();
	const TString dirname = Form("PWGDQ_LMEE_ElectronEfficiencyV2_Cen%d_%d_kINT7",CenMin,CenMax);
  mgr->AddTask(task);
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("Efficiency_Cen%d_%d_kINT7",CenMin,CenMax), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",fileName.Data(),dirname.Data())));
  return task;
}

