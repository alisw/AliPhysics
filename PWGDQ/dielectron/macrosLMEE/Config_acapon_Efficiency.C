Bool_t SetGeneratedSmearingHistos = kFALSE;


Bool_t GetCentralityFromAlien = kFALSE;
std::string centralityFilename = "";
std::string centralityFilenameFromAlien = "/alice/cern.ch/user/a/acapon/.root";

const Int_t triggerNames = AliVEvent::kINT7;

const Int_t nMCSignal   = 0;
const Int_t nCutsetting = 0;

const Double_t minGenPt  = 0.05;
const Double_t maxGenPt  = 20;
const Double_t minGenEta = -1.5;
const Double_t maxGenEta = 1.5;

const Double_t minPtCut  = 0.2;
const Double_t maxPtCut  = 15.0;
const Double_t minEtaCut = -0.8;
const Double_t maxEtaCut = 0.8;
// const Double_t minPtCut = 0.2;
// const Double_t maxPtCut = 8.0;
// const Double_t minEtaCut = -0.8;
// const Double_t maxEtaCut = 0.8;


// binning of single leg histograms
Bool_t usePtVector = kTRUE;
Double_t ptBins[] = {0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,
                     0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,
                     0.90,0.95,1.00,1.10,1.20,1.30,1.40,1.50,1.60,
                     1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,
                     4.00,5.00,6.00,7.00,8.00,10.0,15.0};
const Int_t nBinsPt =  ( sizeof(ptBins) / sizeof(ptBins[0]) )-1;

const Double_t minPtBin   = 0;
const Double_t maxPtBin   = 20;
const Int_t    stepsPtBin = 800;

const Double_t minEtaBin   = -1.0;
const Double_t maxEtaBin   = 1.0;
const Int_t    stepsEtaBin = 20;

const Double_t minPhiBin   = 0;
const Double_t maxPhiBin   = 6.3;
const Int_t    stepsPhiBin = 20;

const Double_t minThetaBin   = 0;
const Double_t maxThetaBin   = TMath::TwoPi();
const Int_t    stepsThetaBin = 60;

const Double_t minMassBin     = 0;
const Double_t maxMassBin     = 5;
const Int_t    stepsMassBin   = 500;
const Double_t minPairPtBin   = 0;
const Double_t maxPairPtBin   = 10;
const Int_t    stepsPairPtBin = 100;

// Binning of resolution histograms
const Int_t    NbinsDeltaMom   = 2000;
const Double_t DeltaMomMin     = -10.0;
const Double_t DeltaMomMax     = 10.0;
const Int_t    NbinsRelMom     = 400;
const Double_t RelMomMin       = 0.0;
const Double_t RelMomMax       = 2.0;
const Int_t    NbinsDeltaEta   = 200;
const Double_t DeltaEtaMin     = -0.4;
const Double_t DeltaEtaMax     = 0.4;
const Int_t    NbinsDeltaTheta = 200;
const Double_t DeltaThetaMin   = -0.4;
const Double_t DeltaThetaMax   = 0.4;
const Int_t    NbinsDeltaPhi   = 200;
const Double_t DeltaPhiMin     = -0.4;
const Double_t DeltaPhiMax     = 0.4;

void GetCentrality(const Int_t centrality, Double_t& CentMin, Double_t& CentMax){
  std::cout << "GetCentrality with centrality " << centrality << std::endl;
  if     (centrality == 0){CentMin = 0;  CentMax = 100;}
  else if(centrality == 1){CentMin = 0;  CentMax = 20;}
  else if(centrality == 2){CentMin = 20; CentMax = 40;}
  else if(centrality == 3){CentMin = 40; CentMax = 60;}
  else if(centrality == 4){CentMin = 60; CentMax = 100;}
  else if(centrality == 5){CentMin = 60; CentMax = 80;}
  else if(centrality == 6){CentMin = 80; CentMax = 100;}
  else if(centrality == 7){CentMin = 0;  CentMax = 5;}
  else if(centrality == 8){CentMin = -1; CentMax = -1;}
  else                      {std::cout << "WARNING::Centrality range not found....." std::endl;}
  return;
}

void ApplyPIDpostCalibration(AliAnalysisTaskElectronEfficiencyV2* task, Int_t whichDet, Bool_t wSDD){
  std::cout << task << std::endl;

  std::cout << "starting ApplyPIDpostCalibration()\n";
  if(whichDet == 0 && wSDD){// ITS
    std::cout << "Loading ITS correction" << std::endl;
    TString localPath = "/home/aaron/Data/diElec_framework_output/PIDcalibration/";
    TString fileName  = "outputITS_MC.root";

    TFile* inFile = TFile::Open(localPath+fileName);
    if(!inFile){
      gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/"+fileName+" .");
      std::cout << "Copy ITS correction from Alien" << std::endl;
      inFile = TFile::Open(fileName);
    }
    else {
      std::cout << "Correction loaded" << std::endl;
    }

    TH3D* mean = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));

    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kITS, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction   (AliAnalysisTaskElectronEfficiencyV2::kITS, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  }
  if(whichDet == 1){// TOF
    std::cout << "Loading TOF correction" << std::endl;
    TString localPath = "/home/aaron/Data/diElec_framework_output/PIDcalibration/";
    TString fileName = "outputTOF";
    if(wSDD == kTRUE){
      fileName.Append("_MC.root");
    }else{
      fileName.Append("_woSDD_MC.root");
    }

    TFile* inFile = TFile::Open(localPath+fileName);
    if(!inFile){
      gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/"+fileName+" .");
      std::cout << "Copy TOF correction from Alien" << std::endl;
      inFile = TFile::Open(fileName);
    }
    else {
      std::cout << "Correction loaded" << std::endl;
    }

    TH3D* mean = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));

    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction   (AliAnalysisTaskElectronEfficiencyV2::kTOF, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  }
}

// #########################################################
// #########################################################

AliAnalysisFilter* SetupTrackCutsAndSettings(TString cutDefinition, Bool_t wSDD)
{
  std::cout << "SetupTrackCutsAndSettings( cutInstance = " << cutDefinition << " )" <<std::endl;
  AliAnalysisFilter* anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); 

  LMEECutLib* LMcutlib = new LMEECutLib(wSDD);
  if(cutDefinition == "kResolutionCuts"){
    std::cout << "Resolution Track Cuts being set" << std::endl;
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kResolutionTrackCuts, LMEECutLib::kResolutionTrackCuts));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kCutSet1"){ // Typical analysis cut
    std::cout << "Setting up cut set 1" << std::endl;
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kCutSet1));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kCutSet2"){ // kCutSet1 w/o ITSshared cut
    std::cout << "Setting up cut set 2" << std::endl;
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet2, LMEECutLib::kCutSet1));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kCutSet3"){ // kCutSet1 with pT <  10 GeV
    std::cout << "Setting up cut set 3" << std::endl;
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet3, LMEECutLib::kCutSet1));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kCutSet4"){ // Hadron band rejection for wSDD data
    std::cout << "Setting up kCutSet4" << std::endl;
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kScheidCuts));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kTheoPID"){ // PID cut set from a Run 1 pPb analysis. Standard track cuts
    std::cout << "Setting up Theo PID. Standard track cuts." << std::endl;
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kTheoPID));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kScheidCuts"){
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kScheidCuts, LMEECutLib::kScheidCuts));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  // Cut sets to to vary fITSshared cut
  else if(cutDefinition == "kITSshared1"){
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared1, LMEECutLib::kScheidCuts));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kITSshared2"){
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared2, LMEECutLib::kScheidCuts));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kITSshared3"){
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared3, LMEECutLib::kScheidCuts));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kITSshared4"){
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared4, LMEECutLib::kScheidCuts));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kITSshared5"){
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared5, LMEECutLib::kScheidCuts));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kITSshared6"){
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared6, LMEECutLib::kScheidCuts));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kITSmin"){
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSmin, LMEECutLib::kScheidCuts));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else if(cutDefinition == "kITSmax"){
    anaFilter->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSmax, LMEECutLib::kScheidCuts));
    anaFilter->SetName(cutDefinition);
    anaFilter->Print();
  }
  else{
    std::cout << "Undefined cut definition...." << std::endl;
    return 0x0;
  }

  return anaFilter;
}

// #########################################################
// #########################################################
std::vector<Bool_t> AddSingleLegMCSignal(AliAnalysisTaskElectronEfficiencyV2* task){

  // SetLegPDGs() requires two pdg codes. For single tracks a dummy value is
  // passed, "1".


  // All final state electrons (excluding conversion electrons)
  AliDielectronSignalMC eleFinalState("eleFinalState","eleFinalState");
  eleFinalState.SetLegPDGs(11,1);
  eleFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalState.SetMotherPDGs(22, 22, kTRUE, kTRUE); // Exclude conversion electrons

  // Electrons from open charm mesons and baryons
  AliDielectronSignalMC eleFinalStateFromD("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD.SetLegPDGs(11,1);
  eleFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD.SetMotherPDGs(402, 402);
  eleFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);

  // Electrons from open beauty mesons and baryons
  AliDielectronSignalMC eleFinalStateFromB("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB.SetLegPDGs(11,1);
  eleFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB.SetMotherPDGs(502, 502);
  eleFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);

  // Add signals
  task->AddSingleLegMCSignal(eleFinalState);
  task->AddSingleLegMCSignal(eleFinalStateFromD);
  task->AddSingleLegMCSignal(eleFinalStateFromB);

  // This is used to get electrons not from same mother for pair efficiency.
  // Needed to look at D and B meson electrons as functionality to pair those is
  // not implemented in the framework. Instead, use all final start electrons
  // from D or B decays for efficiency correction, for example.
  // The ordering must match the ordering of the added signals above*.
  std::vector<Bool_t> DielectronsPairNotFromSameMother;
  DielectronsPairNotFromSameMother.push_back(kFALSE);
  DielectronsPairNotFromSameMother.push_back(kTRUE);
  DielectronsPairNotFromSameMother.push_back(kTRUE);

  return DielectronsPairNotFromSameMother;
}


// #########################################################
// #########################################################
void AddPairMCSignal(AliAnalysisTaskElectronEfficiencyV2* task){

  // Dielectron pairs from same mother (excluding conversions)
  AliDielectronSignalMC pair_sameMother("sameMother","sameMother");
  pair_sameMother.SetLegPDGs(11,-11);
  pair_sameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_sameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  // Set mother properties
  pair_sameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
  pair_sameMother.SetMotherPDGs(22,22,kTRUE,kTRUE); // Exclude conversion

  // Used pdg codes (defined in AliDielectronMC::ComparePDG)
  // 401: open charm meson
  // 404: charged open charmed mesons NO s quark
  // 405: neutral open charmed mesons
  // 406: charged open charmed mesons with s quark
  // 501: open beauty mesons
  // 503: all beauty hadrons
  // 504: charged open beauty mesons NO s quark
  // 505: neutral open beauty mesons
  // 506: charged open beauty mesons with s quark
  // all D mesons

  // decay channels
  // (1) D -> e X
  // (1) B -> e X
  // (2) B -> D X -> e X Y
  // (3) B -> e D X -> ee X Y always produces ULS pair

  // Electrons from open beauty mesons and baryons
  AliDielectronSignalMC eleFinalStateFromB("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB.SetLegPDGs(11,-11);
  eleFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB.SetMotherPDGs(502, 502);
  eleFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);
  eleFinalStateFromB.SetCheckCorrelatedHF(kTRUE);

  // Electrons from open charm mesons and baryons
  AliDielectronSignalMC eleFinalStateFromD("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD.SetLegPDGs(11,-11);
  eleFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD.SetMotherPDGs(402, 402);
  eleFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  eleFinalStateFromD.SetCheckCorrelatedHF(kTRUE);

  AliDielectronSignalMC eleFromJPsi("eleFromJPsi", "eleFromJPsi");
  eleFromJPsi.SetLegPDGs(11,-11);
  eleFromJPsi.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFromJPsi.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFromJPsi.SetMotherPDGs(443, 443);
  eleFromJPsi.SetMothersRelation(AliDielectronSignalMC::kSame);
  eleFromJPsi.SetCheckBothChargesMothers(kTRUE,kTRUE);

  task->AddPairMCSignal(pair_sameMother);
  task->AddPairMCSignal(eleFinalStateFromD);
  task->AddPairMCSignal(eleFinalStateFromB);
  task->AddPairMCSignal(eleFromJPsi);
  // task->AddPairMCSignal(pair_sameMother_pion);
  // task->AddPairMCSignal(pair_sameMother_eta);
  // task->AddPairMCSignal(pair_sameMother_etaP);
  // task->AddPairMCSignal(pair_sameMother_rho);
  // task->AddPairMCSignal(pair_sameMother_omega);
  // task->AddPairMCSignal(pair_sameMother_phi);
  // task->AddPairMCSignal(pair_sameMother_jpsi);
}
