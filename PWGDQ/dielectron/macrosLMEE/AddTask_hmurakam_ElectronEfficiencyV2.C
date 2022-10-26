AliAnalysisTaskElectronEfficiencyV2* AddTask_hmurakam_ElectronEfficiencyV2(
TString name             = "test",
// 0=all sources, 1=HS, 2=Jpsi
Int_t whichGen           = 0,
Bool_t getFromAlien      = kFALSE,
TString configFile       = "Config_hmurakam_ElectronEfficiencyV2.C",
TString calibFileNameTOF = "alien:///alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/calibLMEE/calMaps_TOF_mc.root",
//Binning of Pt,Mass,Ptee
Bool_t usePtVector       = kTRUE,
Bool_t useMassVector     = kFALSE,
Bool_t usePteeVector     = kFALSE,
Bool_t DoULSLS           = kTRUE,
Bool_t DeactivateLS      = kFALSE,
TString year             = "16",
Bool_t usePhiV           = kFALSE,
Double_t maxMee          = 0.04,
Double_t minphiv         = 2.0,
TString suffix           = "",
TString outputFileName   = "LMEE.root",
// Binning of resolution histograms
Int_t NbinsDeltaMom      = 1,
Int_t NbinsRelMom        = 1,
Int_t NbinsDeltaEta      = 1,
Int_t NbinsDeltaTheta    = 1,
Int_t NbinsDeltaPhi      = 1
)
{

  //Fiducial cut
  Float_t minPtCut         = 0.2;
  Float_t maxPtCut         = 15;
  Float_t minEtaCut        = -0.8;
  Float_t maxEtaCut        = +0.8;

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_hmurakam_ElectronEfficiencyV2", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  if (!gSystem->AccessPathName(configFile)) {
    printf("Configfile already present\n");
    configBasePath=Form("%s/",gSystem->pwd());
  }
  else if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/macrosLMEE/%s file:./",configFile.Data()))) ){
    printf("Copy Configfile from alien\n");
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+configFile);
  std::cout << "Configpath:  " << configFilePath << std::endl;

  // Loading config and cutlib
  if (!gROOT->GetListOfGlobalFunctions()->FindObject("Config_hmurakam_ElectronEfficiencyV2.C")) {
    printf("Load macro now\n");
    gROOT->LoadMacro(configFilePath.Data());
  }

  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(name.Data());
  gROOT->GetListOfSpecials()->Add(task);//this is only for ProcessLine(AddMCSignal);

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
    std::cout << "Generator names specified -> Jpsi2ee_1 and B2Jpsi2ee_1" << std::endl;
    task->SetGeneratorMCSignalName("Jpsi2ee_1;B2Jpsi2ee_1");
    task->SetGeneratorULSSignalName("Jpsi2ee_1;B2Jpsi2ee_1");
  }

  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE);
  task->SetTriggerMask(AliVEvent::kINT7);
  task->SetEventFilter((reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine("GetEventCuts()"))));

  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  // Do not set here your analysis pt-cuts
  task->SetMinPtGen(0.1);
  task->SetMaxPtGen(100);
  task->SetMinEtaGen(-1.5);
  task->SetMaxEtaGen(+1.5);

  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  task->SetKinematicCuts(minPtCut, maxPtCut, minEtaCut, maxEtaCut);

  // Set Binning
  // Pt
  const Double_t ptBins[] = {
    //    0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
    //    1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,4.00,5.0,6.0,7.0,10.0,20.0};
    0.0,
    0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
    0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4,
    0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6,
    0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8,
    0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0,
    1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2,
    1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4,
    1.42, 1.44, 1.46, 1.48, 1.5, 1.52, 1.54, 1.56, 1.58, 1.6,
    1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74, 1.76, 1.78, 1.8,
    1.82, 1.84, 1.86, 1.88, 1.9, 1.92, 1.94, 1.96, 1.98, 2.0,
    2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5.0, 10.0, 20.0};

    const Int_t nBinsPt =  ( sizeof(ptBins) / sizeof(ptBins[0]) )-1;
  if (usePtVector == true) {
    std::vector<double> ptBinsVec;
    for (Int_t i = 0; i < nBinsPt+1; ++i){
      ptBinsVec.push_back(ptBins[i]);
    }
    task->SetPtBins(ptBinsVec);
  }
  else task->SetPtBinsLinear(0,  8, 800);
  task->SetEtaBinsLinear(-1.0, 1.0, 40);
  task->SetPhiBinsLinear(0, TMath::TwoPi(), 90);
  task->SetThetaBinsLinear(0, TMath::TwoPi(), 60);

  // Mass
  const Double_t massBins[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  const Int_t nBinsMass =  ( sizeof(massBins) / sizeof(massBins[0]) )-1;
  if (useMassVector == true) {
    std::vector<double> massBinsVec;
    for (Int_t i = 0; i < nBinsMass+1; ++i){
      massBinsVec.push_back(massBins[i]);
    }
    task->SetMassBins(massBinsVec);
  }
  else task->SetMassBinsLinear   (0,  4, 800);//Default

  // Pair ptee  2022/03/17
  const Double_t pteeBins[] = {0.0,
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
    2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
    3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
    5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0,
    7.0, 8.0, 9.0, 10.0};
  if (usePteeVector == true) {
    const Int_t nBinsPtee = ( sizeof(pteeBins) / sizeof(pteeBins[0]) )-1;
    std::vector<double> pteeBinsVec;
    for (Int_t i = 0; i < nBinsPtee+1; ++i){
      pteeBinsVec.push_back(pteeBins[i]);
    }
    task->SetPairPtBins(pteeBinsVec);
  }
  else task->SetPairPtBinsLinear   (0,10,50);

  TGrid::Connect("alien://");
  //  Resolution File, If resoFilename = "" no correction is applied
  std::string resoFilename = Form("%s.root",year.Data());
  task->SetResolutionFile(resoFilename,"/alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/resolution/" + resoFilename);
  task->SetResolutionDeltaPtBinsLinear   (-10.0, 2.0, 1);
  task->SetResolutionRelPtBinsLinear   (0., 2.0, NbinsRelMom);
  task->SetResolutionEtaBinsLinear  (-0.4, 0.4, NbinsDeltaEta);
  task->SetResolutionPhiBinsLinear  (-0.4, 0.4, NbinsDeltaPhi);
  task->SetResolutionThetaBinsLinear(-0.4, 0.4, NbinsDeltaTheta);

  // Set centrality correction. If resoFilename = "" no correction is applied
  // task->SetCentralityFile(centralityFilename);
  // Pairing related config
  task->SetDoPairing(kTRUE);
  task->SetULSandLS(DoULSLS);
  task->SetDeactivateLS(DeactivateLS);
  task->SetPhiVBinsLinear(0, TMath::Pi(), 180);
  task->SetFillPhiV(kFALSE);

  //Set Phiv Cut
  task->SetPhiVCut(usePhiV,maxMee,minphiv);

  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  //  AddSingleLegMCSignal(task);
  AliDielectronSignalMC eleFinalState("eleFinalState","eleFinalState");
  eleFinalState.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  task->AddSingleLegMCSignal(eleFinalState);

  //
  AliDielectronSignalMC eleFinalStateFromD("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD.SetMotherPDGs(402, 402); // open charm mesons and baryons together
  eleFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  task->AddSingleLegMCSignal(eleFinalStateFromD);
  //
  AliDielectronSignalMC eleFinalStateFromB("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB.SetMotherPDGs(502, 502); // open charm mesons and baryons together
  eleFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);
  task->AddSingleLegMCSignal(eleFinalStateFromB);

  // this is used to get electrons from charmed mesons in a environment where GEANT is doing the decay of D mesons, like in LHC18b5a
  // ordering is according to MCSignals of single legs
  printf("Init the DielectronsPairNotFromSameMother vector\n");
  std::vector<bool> DielectronsPairNotFromSameMother;
  DielectronsPairNotFromSameMother.push_back(false);
  DielectronsPairNotFromSameMother.push_back(false);
  DielectronsPairNotFromSameMother.push_back(false);
  task->AddMCSignalsWhereDielectronPairNotFromSameMother(DielectronsPairNotFromSameMother);

  if(whichGen == 0 || whichGen == 2){
    //    AddPairMCSignalLFJPsi(task);
    AliDielectronSignalMC pair_sameMother("sameMother","sameMother");
    pair_sameMother.SetLegPDGs(11,-11);
    pair_sameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    pair_sameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother.SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons. should have no effect on final state ele.
    AliDielectronSignalMC eleFromJPsi("eleFromJPsi", "eleFromJPsi");
    eleFromJPsi.SetLegPDGs(11,-11);
    eleFromJPsi.SetCheckBothChargesLegs(kTRUE,kTRUE);
    eleFromJPsi.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    eleFromJPsi.SetMotherPDGs(443, 443);
    eleFromJPsi.SetMothersRelation(AliDielectronSignalMC::kSame);
    eleFromJPsi.SetCheckBothChargesMothers(kTRUE,kTRUE);
    task->AddPairMCSignal(pair_sameMother);
    task->AddPairMCSignal(eleFromJPsi);
  }else if(whichGen == 1){
    //    AddPairMCSignalHF(task);
    AliDielectronSignalMC pair_sameMother("sameMother","sameMother");
    pair_sameMother.SetLegPDGs(11,-11);
    pair_sameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    pair_sameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother.SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons. should have no effect on final state ele.
    task->AddPairMCSignal(pair_sameMother);
  }else {
    printf("no PairMCSignal added\n");
  };


  // PID postcalibration
  // TOF
  gSystem->Exec(Form("alien_cp %s file:./",calibFileNameTOF.Data()));
  TH2F* histMean2DTOF  = 0x0;
  TH2F* histWidth2DTOF = 0x0;
  TFile *rootfileTOF   = 0x0;
  printf("reading : %s for TOF PID calibration\n",calibFileNameTOF.Data());
  rootfileTOF = TFile::Open(calibFileNameTOF,"READ");
  histMean2DTOF  = (TH2F*)rootfileTOF->Get(Form("m%s",year.Data()));
  histWidth2DTOF = (TH2F*)rootfileTOF->Get(Form("w%s",year.Data()));
  printf("%s and %s\n",Form("m%s",year.Data()),Form("w%s",year.Data()));

  for (Int_t i = 0; i <= histMean2DTOF->GetNbinsX()+1; i++){
    for (Int_t k = 0; k <= histMean2DTOF->GetNbinsY()+1; k++){
      if ( (i == 0) || (k == 0) || (i > histMean2DTOF->GetNbinsX()) || (k > histMean2DTOF->GetNbinsY())) { // under/overflows
	histMean2DTOF->SetBinContent(i, k, 0.0 );
	histWidth2DTOF->SetBinContent(i, k, 1.0 );
      }
    }
  }

  cout<<"Adding mean & width TOF PID correction" <<endl;
  task->SetWidthCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, histWidth2DTOF, AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, histMean2DTOF, AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);


  // Adding cutsettings
  // Cuts
  // Number of cuts
  const Int_t nDie = (Int_t)gROOT->ProcessLine("GetN()");
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliAnalysisFilter *filter = reinterpret_cast<AliAnalysisFilter*>(gROOT->ProcessLine(Form("Config_hmurakam_ElectronEfficiencyV2(%d)",i)));
    task->AddTrackCuts(filter);
  }

  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName = outputFileName; // create a subfolder in the file
  TString outlistname = Form("efficiency%s",suffix.Data());
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(outlistname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;

}
