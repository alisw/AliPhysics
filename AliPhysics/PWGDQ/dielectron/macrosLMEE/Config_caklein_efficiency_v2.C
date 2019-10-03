
// TString names=("kPbPb2015_Pt400_looseTOFif;kPbPb2015_Pt400_tightTOFreq;NewPID;kPbPb2015_Pt400_looseTOFif_0SITSCl;kPbPb2015_Pt400_tightTOFreq_0SITSCl;NewPID_0SITSCl;kPbPb2015_Pt400_looseTOFif_ITSMAP;kPbPb2015_Pt400_tightTOFreq_ITSMAP;NewPID_ITSMAP");
// TString names=("cut1_pt400;cut2_pt400;cut3_pt400;cut4_pt400;cut5_pt400;cut6_pt400;cut7_pt400;cut8_pt400;cut9_pt400;cut10_pt400;cut11_pt400;cut12_pt400;cut13_pt400;cut14_pt400;cut15_pt400;cut16_pt400;cut17_pt400;cut18_pt400;cut19_pt400;cut20_pt400");
// TString names=("cut5_pt200;kPbPb2015_Pt200_PID_cutoff_pion_kaon_proton;cut5_pt200_TOFreq");
TString names=("cut5_pt200_noPIDcuts;cut5_pt200;cut5_pt200_woSCcut");

// TString names=("cut1_pt200;cut2_pt200;cut3_pt200;cut4_pt200;cut5_pt200;cut6_pt200;cut7_pt200;cut8_pt200;cut9_pt200;cut10_pt200;cut11_pt200;cut12_pt200;cut13_pt200;cut14_pt200;cut15_pt200;cut16_pt200;cut17_pt200;cut18_pt200;cut19_pt200;cut20_pt200");
// TString names=("cut1_pt200_woSCcut;cut2_pt200_woSCcut;cut3_pt200_woSCcut;cut4_pt200_woSCcut;cut5_pt200_woSCcut;cut6_pt200_woSCcut;cut7_pt200_woSCcut;cut8_pt200_woSCcut;cut9_pt200_woSCcut;cut10_pt200_woSCcut;cut11_pt200_woSCcut;cut12_pt200_woSCcut;cut13_pt200_woSCcut;cut14_pt200_woSCcut;cut15_pt200_woSCcut;cut16_pt200_woSCcut;cut17_pt200_woSCcut;cut18_pt200_woSCcut;cut19_pt200_woSCcut;cut20_pt200_woSCcut");
// TString names=("cut1_pt200_wSCcut;cut2_pt200_wSCcut;cut3_pt200_wSCcut;cut4_pt200_wSCcut;cut5_pt200_wSCcut;cut6_pt200_wSCcut;cut7_pt200_wSCcut;cut8_pt200_wSCcut;cut9_pt200_wSCcut;cut10_pt200_wSCcut;cut11_pt200_wSCcut;cut12_pt200_wSCcut;cut13_pt200_wSCcut;cut14_pt200_wSCcut;cut15_pt200_wSCcut;cut16_pt200_wSCcut;cut17_pt200_wSCcut;cut18_pt200_wSCcut;cut19_pt200_wSCcut;cut20_pt200_wSCcut");


// TString names=("cut5_openKinematicCuts");
// TString names=("cut5_openKinematicCuts_without_CrossedOverFindableCut");
// TString names=("kPbPb2015_Pt100_ResolutionCuts");

// TString generatorName = "Hijing";

// TString generatorName = "";
// TString generatorName = "Hijing_0;pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7;Pythia CC_8;Pythia B_8;Pythia BB_8";
// TString generatorName = "pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7;Pythia CC_8;Pythia B_8;Pythia BB_8";
// TString generatorName = "Hijing_0;Hijing";
// TString generatorName = "pizero_1";
// TString generatorName = "eta_2";
// TString generatorName = "etaprime_3";
// TString generatorName = "rho_4";
// TString generatorName = "omega_5";
// TString generatorName = "phi_6";
// TString generatorName = "jpsi_7";
// TString generatorName = "Pythia_CC_";

// TString generatorName = "Hijing_0";
// TString generatorName = "Jpsi2ee_1";

TString generatorNameForMCSignal  = "pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7";
// TString generatorNameForMCSignal  = "Hijing_0";
TString generatorNameForULSSignal = "Hijing;pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7;Hijing_0;Pythia CC_8;Pythia B_8;Pythia BB_8";


Bool_t SetTPCCorrection = kFALSE;
Bool_t SetITSCorrection = kTRUE;
Bool_t SetTOFCorrection = kTRUE;

bool SetGeneratedSmearingHistos = false;

bool DoPairing = true;
bool DoULSLS   = true;
bool DeactivateLS   = true;

bool GetResolutionFromAlien = kTRUE;
// std::string resoFilename = "resolution_PbPb2015_0080_deltaXvsP_cut5_noKinematicCuts.root";
// std::string resoFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/resolution_PbPb2015_0080_deltaXvsP_cut5_noKinematicCuts.root";
std::string resoFilename = "resolution_PbPb2015_0080_deltaXvsP.root";
std::string resoFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/resolution_PbPb2015_0080_deltaXvsP.root";
// std::string resoFilename = "new_resolution_PbPb2015_0080_deltaXvsP.root";
// std::string resoFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/new_resolution_PbPb2015_0080_deltaXvsP.root";

bool DoCocktailWeighting = kTRUE;
bool GetCocktailFromAlien = kTRUE;
std::string CocktailFilename = "Cocktail_PbPb0080_5TeV_MotherPt.root";
std::string CocktailFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/Cocktail_PbPb0080_5TeV_MotherPt.root";
// std::string CocktailFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/cocktail_PbPb0080_5TeV.root";

bool GetCentralityFromAlien = kFALSE;
std::string centralityFilename = "";
std::string centralityFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/centralityLHC16g1.root";

const Int_t triggerNames = AliVEvent::kMB;

const int nMCSignal = 0;
const int nCutsetting = 0;


const double minGenPt = 0.05;
const double maxGenPt = 100;
const double minGenEta = -1.5;
const double maxGenEta =  1.5;


const double minPtCut = 0.2;
const double maxPtCut = 8.0;
const double minEtaCut = -0.8;
const double maxEtaCut =  0.8;
// const double minPtCut = 0.2;
// const double maxPtCut = 8.0;
// const double minEtaCut = -0.8;
// const double maxEtaCut = 0.8;



// binning of single leg histograms
bool usePtVector = true;
double ptBins[] = {0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
  1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,
  4.00,5.0,6.0,7.0,8.0
  };
const Int_t nBinsPt =  ( sizeof(ptBins) / sizeof(ptBins[0]) )-1;

const double minPtBin = 0;
const double maxPtBin = 8;
const int    stepsPtBin = 400;

const double minEtaBin = -1.0;
const double maxEtaBin =  1.0;
const int    stepsEtaBin = 20;

const double minPhiBin = 0;
const double maxPhiBin =  6.3;
const int    stepsPhiBin = 20;

const double minThetaBin = 0;
const double maxThetaBin =  TMath::TwoPi();
const int    stepsThetaBin = 60;

const double minMassBin = 0;
const double maxMassBin =  5;
const int    stepsMassBin = 250;
const double minPairPtBin = 0;
const double maxPairPtBin =  8;
const int    stepsPairPtBin = 160;

// Binning of resolution histograms
const int    NbinsDeltaMom    = 2000;
const double DeltaMomMin   =-10.0;
const double DeltaMomMax   = 10.0;
const int    NbinsRelMom      = 400;
const double RelMomMin     =  0.0;
const double RelMomMax     =  2.0;
const int    NbinsDeltaEta    = 200;
const double DeltaEtaMin   = -0.4;
const double DeltaEtaMax   =  0.4;
const int    NbinsDeltaTheta  = 200;
const double DeltaThetaMin = -0.4;
const double DeltaThetaMax =  0.4;
const int    NbinsDeltaPhi    = 200;
const double DeltaPhiMin   = -0.4;
const double DeltaPhiMax   =  0.4;

Int_t centrality = LMEECutLib::kPbPb0080; // DUMMY!!!!!!! will be overwritten

void GetCentrality(const int centrality, double& CentMin, double& CentMax){
  std::cout << "GetCentrality with centrality " << centrality << std::endl;
  if      (centrality == 0) {CentMin = 0; CentMax = 10;}
  else if (centrality == 1) {CentMin = 10; CentMax = 50;}
  else if (centrality == 2) {CentMin = 50; CentMax = 80;}
  else if (centrality == 3) {CentMin = 0;  CentMax = 80;}
  else if (centrality == 4) {CentMin = 0;  CentMax = 10;}
  else if (centrality == 5) {CentMin = 10; CentMax = 20;}
  else if (centrality == 6) {CentMin = 20; CentMax = 40;}
  else if (centrality == 7) {CentMin = 40; CentMax = 50;}
  else if (centrality == 8) {CentMin = 50; CentMax = 60;}
  else if (centrality == 9) {CentMin = 60; CentMax = 80;}

  else if (centrality == 10) {CentMin = 0; CentMax = 10;}
  else if (centrality == 11) {CentMin = 10; CentMax = 20;}
  else if (centrality == 12) {CentMin = 20; CentMax = 30;}
  else if (centrality == 13) {CentMin = 30; CentMax = 40;}
  else if (centrality == 14) {CentMin = 40; CentMax = 50;}
  else if (centrality == 15) {CentMin = 50; CentMax = 60;}
  else if (centrality == 16) {CentMin = 60; CentMax = 70;}
  else if (centrality == 17) {CentMin = 70; CentMax = 80;}
  else if (centrality == 18) {CentMin = 80; CentMax = 90;}

  else if (centrality == 20) {CentMin = 0; CentMax = 1;}
  else if (centrality == 21) {CentMin = 1; CentMax = 2;}
  else if (centrality == 22) {CentMin = 2; CentMax = 3;}
  else if (centrality == 23) {CentMin = 3; CentMax = 4;}
  else if (centrality == 24) {CentMin = 4; CentMax = 5;}
  else if (centrality == 25) {CentMin = 5; CentMax = 6;}
  else if (centrality == 26) {CentMin = 6; CentMax = 7;}
  else if (centrality == 27) {CentMin = 7; CentMax = 8;}
  else if (centrality == 28) {CentMin = 8; CentMax = 9;}
  else if (centrality == 29) {CentMin = 9; CentMax = 10;}

  else if (centrality == 30) {CentMin = 0; CentMax = 5;}
  else if (centrality == 31) {CentMin = 5; CentMax = 10;}
  else if (centrality == 32) {CentMin = 10; CentMax = 15;}
  else if (centrality == 33) {CentMin = 15; CentMax = 20;}

  else                      {CentMin = 0;  CentMax = 0;}
  return;
}

void DoAdditionalWork(AliAnalysisTaskElectronEfficiencyV2* task){
  std::cout << task << std::endl;

  std::cout << "starting DoAdditionalWork()\n";
  if (SetTPCCorrection == true){
    std::cout << "Loading TPC correction" << std::endl;
    std::string file_name = "recalib_mc_tpc_nsigmaele.root";
    TFile* _file = TFile::Open(file_name.c_str());

    if (!_file){
      gSystem->Exec(("alien_cp alien:///alice/cern.ch/user/c/cklein/data/recalibration/" + file_name + " .").c_str());
      std::cout << "Copy TPC correction from Alien" << std::endl;
      _file = TFile::Open(file_name.c_str());
    }
    else {
      std::cout << "Correction loaded" << std::endl;
    }

    TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));

    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTPC, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction   (AliAnalysisTaskElectronEfficiencyV2::kTPC, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  }
  if (SetITSCorrection == true){
    std::cout << "Loading ITS correction" << std::endl;
    std::string file_name = "recalib_mc_its_nsigmaele.root";
    TFile* _file = TFile::Open(file_name.c_str());

    if (!_file){
      gSystem->Exec(("alien_cp alien:///alice/cern.ch/user/c/cklein/data/recalibration/" + file_name + " .").c_str());
      std::cout << "Copy ITS correction from Alien" << std::endl;
      _file = TFile::Open(file_name.c_str());
    }
    else {
      std::cout << "Correction loaded" << std::endl;
    }

    TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));

    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kITS, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction   (AliAnalysisTaskElectronEfficiencyV2::kITS, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  }
  if (SetTOFCorrection == true){
    std::cout << "Loading TOF correction" << std::endl;
    std::string file_name = "recalib_mc_tof_nsigmaele.root";
    TFile* _file = TFile::Open(file_name.c_str());

    if (!_file){
      gSystem->Exec(("alien_cp alien:///alice/cern.ch/user/c/cklein/data/recalibration/" + file_name + " .").c_str());
      std::cout << "Copy TOF correction from Alien" << std::endl;
      _file = TFile::Open(file_name.c_str());
    }
    else {
      std::cout << "Correction loaded" << std::endl;
    }

    TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));

    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction   (AliAnalysisTaskElectronEfficiencyV2::kTOF, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  }
}

// #########################################################
// #########################################################

AliAnalysisFilter* SetupTrackCutsAndSettings(TString cutDefinition, Bool_t isAOD)
{
  std::cout << "TEST" << std::endl;
  std::cout << "SetupTrackCutsAndSettings( cutInstance = " << cutDefinition << " )" <<std::endl;
  AliAnalysisFilter *anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); // named constructor seems mandatory!

  AnalysisCut AnaCut;


  LMEECutLib* LMcutlib = new LMEECutLib();
  if (cutDefinition == "kPbPb2015_Pt100_ResolutionCuts"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt100_ResolutionCuts);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kResolutionTrackCuts);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5_pt200_noPIDcuts"){
    AnaCut.SetPIDAna(LMEECutLib::kPDGelectron);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5_openKinematicCuts"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt100_ResolutionCuts);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5_openKinematicCuts_without_CrossedOverFindableCut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt100_ResolutionCuts);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5_without_CrossedOverFindableCut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt200_PID_cutoff_pion_kaon_proton"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt200_PID_cutoff_pion_kaon_proton);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut1_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_1_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut2_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_2_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut3_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_3_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_3);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut4_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_4_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_4);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_5_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut6_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_6_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_6);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut7_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_7_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_7);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut8_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_8_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_8);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut9_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_9_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_9);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut10_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_10_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_10);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut11_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_11_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_11);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "cut12_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_12_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_12);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut13_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_13_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_13);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut14_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_14_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_14);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut15_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_15_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_15);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut16_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_16_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_16);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut17_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_17_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_17);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut18_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_18_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_18);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut19_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_19_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_19);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut20_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_20_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_20);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut1_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_1_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut2_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_2_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut3_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_3_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_3);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut4_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_4_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_4);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_5_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut6_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_6_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_6);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut7_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_7_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_7);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut8_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_8_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_8);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut9_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_9_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_9);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut10_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_10_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_10);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut11_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_11_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_11);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "cut12_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_12_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_12);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut13_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_13_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_13);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut14_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_14_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_14);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut15_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_15_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_15);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut16_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_16_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_16);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut17_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_17_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_17);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut18_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_18_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_18);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut19_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_19_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_19);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut20_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_20_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_20);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5_pt200_looserPionRejection"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_5_pt200_looserPionRejection);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5_pt200_looserPionRejection_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_5_pt200_looserPionRejection);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut1_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_1_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut2_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_2_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut3_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_3_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_3_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut4_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_4_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_4_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_5_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut6_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_6_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_6_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut7_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_7_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_7_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut8_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_8_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_8_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut9_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_9_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_9_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut10_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_10_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_10_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut11_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_11_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_11_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut12_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_12_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_12_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut13_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_13_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_13_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut14_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_14_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_14_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut15_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_15_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_15_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut16_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_16_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_16_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut17_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_17_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_17_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut18_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_18_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_18_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut19_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_19_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_19_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut20_pt200_woSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_20_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_20_woSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut1_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_1_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut2_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_2_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut3_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_3_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_3);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut4_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_4_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_4);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_5_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut6_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_6_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_6);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut7_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_7_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_7);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut8_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_8_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_8);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut9_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_9_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_9);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut10_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_10_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_10);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut11_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_11_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_11);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut12_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_12_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_12);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut13_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_13_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_13);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut14_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_14_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_14);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut15_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_15_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_15);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut16_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_16_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_16);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut17_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_17_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_17);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut18_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_18_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_18);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut19_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_19_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_19);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut20_pt200_TOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_20_pt200_TOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_20);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut1_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_1_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut2_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_2_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut3_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_3_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_3_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut4_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_4_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_4_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_5_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut6_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_6_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_6_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut7_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_7_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_7_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut8_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_8_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_8_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut9_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_9_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_9_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut10_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_10_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_10_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut11_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_11_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_11_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut12_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_12_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_12_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut13_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_13_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_13_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut14_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_14_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_14_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut15_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_15_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_15_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut16_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_16_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_16_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut17_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_17_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_17_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut18_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_18_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_18_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut19_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_19_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_19_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut20_pt200_wSCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_20_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_20_wSCcut);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  if (!isAOD) anaFilter->AddCuts( LMcutlib->GetESDTrackCutsAna(AnaCut) );
  anaFilter->AddCuts( LMcutlib->GetPIDCutsAna(AnaCut) );
  anaFilter->SetName(cutDefinition);
  return anaFilter;
}


// #########################################################
// #########################################################
AliAnalysisCuts* SetupEventCuts(Bool_t isAOD)
{
  std::cout << "Setup Event Cuts" << std::endl;
  // event cuts are identical for all analysis 'cutInstance's that run together!
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  if (isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
  eventCuts->Print();
  return eventCuts;
}



// #########################################################
// #########################################################
std::vector<bool> AddSingleLegMCSignal(AliAnalysisTaskElectronEfficiencyV2* task){
  AliDielectronSignalMC partFinalState("partFinalState","partFinalState");
  partFinalState.SetLegPDGs(0,1);//dummy second leg (never MCtrue)\n"
  // partFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  partFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);

  AliDielectronSignalMC eleFinalState("eleFinalState","eleFinalState");
  eleFinalState.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);

  AliDielectronSignalMC eleFinalStateFromSameMotherMeson("eleFinalStateFromSameMotherMeson","eleFinalStateFromSameMotherMeson");
  eleFinalStateFromSameMotherMeson.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromSameMotherMeson.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromSameMotherMeson.SetMotherPDGs(600, 600); // open charm mesons and baryons together
  eleFinalStateFromSameMotherMeson.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //
  AliDielectronSignalMC eleFinalStateFromD("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD.SetMotherPDGs(402, 402); // open charm mesons and baryons together
  eleFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC eleFinalStateFromB("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB.SetMotherPDGs(502, 502); // open charm mesons and baryons together
  eleFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);

  // AliDielectronSignalMC eleSecondary("eleSecondary","eleSecondary");
  // eleSecondary.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  // eleSecondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
  // eleSecondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  // task->AddSingleLegMCSignal(eleSecondary);
  //
  // AliDielectronSignalMC eleDontCare("eleDontCare","eleDontCare");
  // eleDontCare.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  // eleDontCare.SetCheckBothChargesLegs(kTRUE,kTRUE);
  // eleDontCare.SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // task->AddSingleLegMCSignal(eleDontCare);

  // task->AddSingleLegMCSignal(partFinalState);
  task->AddSingleLegMCSignal(eleFinalState);
  // task->AddSingleLegMCSignal(eleFinalStateFromPion);
  task->AddSingleLegMCSignal(eleFinalStateFromD);
  task->AddSingleLegMCSignal(eleFinalStateFromB);
  // task->AddSingleLegMCSignal(eleFinalStateFromSameMotherMeson);

  // this is used to get electrons from charmed mesons in a environment where GEANT is doing the decay of D mesons, like in LHC18b5a
  // ordering is according to MCSignals of single legs
 std::vector<bool> DielectronsPairNotFromSameMother;
 DielectronsPairNotFromSameMother.push_back(false);
 DielectronsPairNotFromSameMother.push_back(true);
 DielectronsPairNotFromSameMother.push_back(true);
 // DielectronsPairNotFromSameMother.push_back(false);
 return DielectronsPairNotFromSameMother;
}


// #########################################################
// #########################################################
void AddPairMCSignal(AliAnalysisTaskElectronEfficiencyV2* task){

    AliDielectronSignalMC pair_sameMother("sameMother","sameMother");
    pair_sameMother.SetLegPDGs(11,-11);
    pair_sameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother.SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons. should have no effect on final state ele.

    AliDielectronSignalMC pair_sameMother_pion("sameMother_pion","sameMother_pion");
    pair_sameMother_pion.SetLegPDGs(11,-11);
    pair_sameMother_pion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_pion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_pion.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_pion.SetMotherPDGs(111,111); //

    AliDielectronSignalMC pair_sameMother_eta("sameMother_eta","sameMother_eta");
    pair_sameMother_eta.SetLegPDGs(11,-11);
    pair_sameMother_eta.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_eta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_eta.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_eta.SetMotherPDGs(221,221); //

    AliDielectronSignalMC pair_sameMother_etaP("sameMother_etaP","sameMother_etaP");
    pair_sameMother_etaP.SetLegPDGs(11,-11);
    pair_sameMother_etaP.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_etaP.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_etaP.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_etaP.SetMotherPDGs(331,331); //

    AliDielectronSignalMC pair_sameMother_rho("sameMother_rho","sameMother_rho");
    pair_sameMother_rho.SetLegPDGs(11,-11);
    pair_sameMother_rho.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_rho.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_rho.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_rho.SetMotherPDGs(113, 113); //

    AliDielectronSignalMC pair_sameMother_omega("sameMother_omega","sameMother_omega");
    pair_sameMother_omega.SetLegPDGs(11,-11);
    pair_sameMother_omega.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_omega.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_omega.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_omega.SetMotherPDGs(223, 223); //

    AliDielectronSignalMC pair_sameMother_phi("sameMother_phi","sameMother_phi");
    pair_sameMother_phi.SetLegPDGs(11,-11);
    pair_sameMother_phi.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_phi.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_phi.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_phi.SetMotherPDGs(333, 333); //

    AliDielectronSignalMC pair_sameMother_jpsi("sameMother_jpsi","sameMother_jpsi");
    pair_sameMother_jpsi.SetLegPDGs(11,-11);
    pair_sameMother_jpsi.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_jpsi.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_jpsi.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_jpsi.SetMotherPDGs(443, 443); //

    AliDielectronSignalMC pair_sameMother_CharmedMesonsWithSameMother("CharmedMesonsWithSameMother","CharmedMesonsWithSameMother");
    pair_sameMother_CharmedMesonsWithSameMother.SetLegPDGs(11,-11);
    pair_sameMother_CharmedMesonsWithSameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_CharmedMesonsWithSameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_CharmedMesonsWithSameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_CharmedMesonsWithSameMother.SetMotherPDGs(402, 402); //

    AliDielectronSignalMC pair_conversion("pair_conversion","pair_conversion");
    pair_conversion.SetLegPDGs(11,-11);
    pair_conversion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_conversion.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_conversion.SetMothersRelation(AliDielectronSignalMC::kSame);



    task->AddPairMCSignal(pair_sameMother);
    // task->AddPairMCSignal(pair_sameMother_pion);
    // task->AddPairMCSignal(pair_sameMother_eta);
    // task->AddPairMCSignal(pair_sameMother_etaP);
    // task->AddPairMCSignal(pair_sameMother_rho);
    // task->AddPairMCSignal(pair_sameMother_omega);
    // task->AddPairMCSignal(pair_sameMother_phi);
    // task->AddPairMCSignal(pair_sameMother_jpsi);
    // task->AddPairMCSignal(pair_sameMother_CharmedMesonsWithSameMother);
    // task->AddPairMCSignal(pair_conversion);
}
