// TRACK CUTS
// ################################################################
// ################# PreFilter Track Cut Primary ##################
// ################################################################
TString names_Prim_Track_PreFilter_Cuts=("JPID_sum_pt75_PreFilter;JPID_sum1_pt75_sec_kV0_PreFilter");
// TString names_Prim_Track_PreFilter_Cuts=("JPID_sum_pt75;JPID_sum1_pt75_sec_kV0");

// ################################################################
// ################# PreFilter Track Cut Secondary ################
// ################################################################
//         !!!!!!!    actually not in use     !!!!!!
TString names_Sec_Track_PreFilter_Cuts=("noPID_V0_PreFilter;track_V0_PreFilter");

// ################################################################
// ################# Standard Track Cut Primary ###################
// ################################################################
// TString names_Prim_Cuts=("noPID");     // still has kin cuts (pt 75 MeV/c)
// TString names_Prim_Cuts=("JPID_sum_pt75");
// TString names_Prim_Cuts=("JPID_sum_pt75;JPID_sum2_pt75_sec_OnlyCosOpenAngle;JPID_sum7_pt75_sec_OnlyM;JPID_sum8_pt75_sec_OnlyArmPt;JPID_sum9_pt75_sec_OnlyArmAlpha;JPID_sum1_pt75_sec_kV0");
// TString names_Prim_Cuts=("JPID_sum_pt75;JPID_sum6_pt75_sec_wCosChi2LegDistRPsiPairM;JPID_sum7_pt75_sec_wCosChi2LegDistRPsiPairMArmPt;JPID_sum8_pt75_sec_wCosChi2LegDistRPsiPairMArmPtAlpha");
TString names_Prim_Track_standard_Cuts=("JPID_sum_pt75;JPID_sum1_pt75_sec_kV0");
// TString names_Prim_Track_standard_Cuts=("JPID_sum_pt75_PreFilter;JPID_sum1_pt75_sec_kV0_PreFilter");

// ################################################################
// ############### Standard Track Cut Secondary ###################
// ################################################################
TString names_Sec_Track_standard_Cuts=("noPID_V0_standard;track_V0_standard");



// PAIR CUTS
// ################################################################
// ################# PreFilter Pair Cut Primary ###################
// ################################################################
//                        *** MISSING ***

// ################################################################
// ################# PreFilter Pair Cut Secondary #################
// ################################################################
// TString names_Sec_Pair_PreFilter_Cuts=("noPID;kV0");
TString names_Sec_Pair_PreFilter_Cuts=("noPID;pairkV0_PreFilter");

// ################################################################
// ################# Standard Pair Cut Primary ####################
// ################################################################
TString names_Prim_Pair_Cuts=("pairJPID_sum_pt20;pairJPID_sum1_pt20_sec_kV0");

// ################################################################
// ################# Standard Pair Cut Secondary ##################
// ################################################################
// TString names_Sec_Pair_standard_Cuts=("noPID");      // still has kin cuts (pt 75 MeV/c)
// TString names_Sec_Pair_standard_Cuts=("JPID_sum_pt75_secondary");
// TString names_Sec_Pair_standard_Cuts=("kV0");
// TString names_Sec_Pair_standard_Cuts=("noPID;kV0OnlyCosOpenAngle;kV0OnlyM;kV0OnlyArmPt;kV0OnlyArmAlpha;kV0");
// TString names_Sec_Pair_standard_Cuts=("noPID;kV0wPairM;kV0wPairMArmPt;kV0wPairMArmPtAlpha");
TString names_Sec_Pair_standard_Cuts=("noPID;pairkV0");
// TString names_Sec_Pair_standard_Cuts=("noPID;kV0_PreFilter");


                                                                                // TString names_Sec_Pair_standard_Cuts=("onlyPIDcut1_noTrackCuts");
                                                                                // TString names_Sec_Pair_standard_Cuts=("cut1_pt75_secondary");
                                                                                // TString names_Sec_Pair_standard_Cuts=("cut1_pt75");
                                                                                // TString names_Sec_Pair_standard_Cuts=("noPID;cut1_pt75_secondary");


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!
// 400 MEV CUT !!!!!!!!!!!!!!!
// TString names=("cut5_pt400;cut11_pt400;cut23_pt400");

// TString names=("cut5_openKinematicCuts");
// TString names=("cut5_openKinematicCuts_without_CrossedOverFindableCut");
// TString names=("kPbPb2015_Pt100_ResolutionCuts");

// TString generatorNameForMCSignal  = "Pythia CC_0;Pythia CC_1;Pythia B_0;Pythia B_1";
// TString generatorNameForULSSignal = "Pythia CC_0;Pythia CC_1;Pythia B_0;Pythia B_1";
TString generatorNameForMCSignal  = "";
TString generatorNameForULSSignal = "";
// TString generatorNameForMCSignal  = "pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7";
// TString generatorNameForULSSignal = "Hijing;pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7;Hijing_0;Pythia CC_8;Pythia B_8;Pythia BB_8";


Bool_t SetTPCCorrection = kFALSE;
Bool_t SetITSCorrection = kFALSE;
Bool_t SetTOFCorrection = kFALSE;


bool debug = true;

bool DoPairing         = true;
bool DoFourPairing     = true;
bool UsePreFilter      = true;
bool UseSecPreFilter   = true;
bool DoMassCut         = true;
bool V0OnFlyStatus     = true; // true stands for OnFlyStatus:aktive ; false means deaktivated
// bool DoULSLS   = true;

bool UseMCDataSig   = true; // if it is selected true the running time is increasing drastically, Reducing time for example by mass cut.

bool GetResolutionFromAlien = kTRUE;
// std::string resoFilename = "resolution_PbPb2015_0080_deltaXvsP_cut5_noKinematicCuts.root";
// std::string resoFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/resolution_PbPb2015_0080_deltaXvsP_cut5_noKinematicCuts.root";
std::string resoFilename = "resolution_pp2016f_deltaXvsP_LHC17d1_vsPt.root";
std::string resoFilenameFromAlien = "/alice/cern.ch/user/f/feisenhu/data/resolution_pp2016f_deltaXvsP_LHC17d1_vsPt.root";
// std::string resoFilename = "resolution_PbPb2015_0080_deltaXvsP.root";
// std::string resoFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/resolution_PbPb2015_0080_deltaXvsP.root";
// std::string resoFilename = "new_resolution_PbPb2015_0080_deltaXvsP.root";
// std::string resoFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/new_resolution_PbPb2015_0080_deltaXvsP.root";

bool DoCocktailWeighting = kFALSE;
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


const double minGenPt = 0.01;
const double maxGenPt = 100;
const double minGenEta = -1.5;
const double maxGenEta =  1.5;


// const double minPtCut = 0.2;
// const double minPtCut = 0.;
// const double maxPtCut = 100.0;
// const double minEtaCut = -100;
// const double maxEtaCut =  100;
// const double minPtCut = 0.200;
const double minPtCut = 0.075;
const double maxPtCut = 8.0;
const double minEtaCut = -0.8;
const double maxEtaCut = 0.8;


// const double upperMassCutPrimaries = 0.547862;
// const double lowerMassCutPrimaries = 0.1349766;
// const double upperMassCutPrimaries = 1.;
const double lowerMassCutPrimaries = 0.1;
const double upperMassCutPrimaries = 0.2;
// const double lowerPrimSecPreFilterMass = 0.1;
// const double upperPrimSecPreFilterMass = 0.165;
const double lowerPrimSecPreFilterMass = 0.03;
const double upperPrimSecPreFilterMass = 0.25;
const double lowerSecSecPreFilterMass = 0.1;
const double upperSecSecPreFilterMass = 0.2;
const double massCutSecondaries = 0.01;
const double photonMass  = 0.0;




// binning of single leg histograms
bool usePtVector = false;
double ptBins[] = {0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
  1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,
  4.00,5.0,6.0,7.0,8.0
  };

// double ptBins[] = {0.000,0.015,0.030,0.045,0.060,0.075,0.090,0.100,0.110,0.120,0.130,0.140,0.150,0.160,0.170,0.180,0.190,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
  // 1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,
  // 4.00,5.0,6.0,7.0,8.0
  // };

const Int_t nBinsPt =  ( sizeof(ptBins) / sizeof(ptBins[0]) )-1;

const double minPtBin = 0;
const double maxPtBin = 8;
//const int    stepsPtBin = 320;
const int    stepsPtBin = 80;
// const int    stepsPtBin = 800;

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
const double maxMassBin =  1;
const int    stepsMassBin = 1000;
const double minPairPtBin = 0;
const double maxPairPtBin =  4;
// const double maxPairPtBin =  8;
 const int    stepsPairPtBin = 80;
// const int    stepsPairPtBin = 400;

// Binning of resolution histograms
const int    NbinsDeltaMom    = 1000;
const double DeltaMomMin   =-10.0;
const double DeltaMomMax   = 10.0;
const int    NbinsRelMom      = 200;
const double RelMomMin     =  0.0;
const double RelMomMax     =  2.0;
const int    NbinsDeltaEta    = 100;
const double DeltaEtaMin   = -0.4;
const double DeltaEtaMax   =  0.4;
const int    NbinsDeltaTheta  = 100;
const double DeltaThetaMin = -0.4;
const double DeltaThetaMax =  0.4;
const int    NbinsDeltaPhi    = 100;
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

void DoAdditionalWork(AliAnalysisTaskEtaReconstruction* task){
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

    task->SetCentroidCorrFunction(AliAnalysisTaskEtaReconstruction::kTPC, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction   (AliAnalysisTaskEtaReconstruction::kTPC, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
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

    task->SetCentroidCorrFunction(AliAnalysisTaskEtaReconstruction::kITS, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction   (AliAnalysisTaskEtaReconstruction::kITS, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
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

    task->SetCentroidCorrFunction(AliAnalysisTaskEtaReconstruction::kTOF, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction   (AliAnalysisTaskEtaReconstruction::kTOF, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  }
}

// #########################################################
// #########################################################

AliAnalysisFilter* SetupTrackCutsAndSettings(TString cutDefinition, Bool_t isAOD)
{
  std::cout << "SetupTrackCutsAndSettings( cutInstance = " << cutDefinition << " )" <<std::endl;
  AliAnalysisFilter *anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); // named constructor seems mandatory!

  AnalysisCut AnaCut;


  LMEECutLib* LMcutlib = new LMEECutLib();

  if (cutDefinition == "noPID"){
    // AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt75);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);

    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }


  else if (cutDefinition == "noPID_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt200_noPID);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "onlyJPID_sum_noTrackCuts"){
    AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  /////////////////////////////////////////////////////////
  //             primary track cut settings              //
  /////////////////////////////////////////////////////////
  else if (cutDefinition == "JPID_sum_pt75" || cutDefinition == "JPID_sum1_pt75_sec_kV0"){
    AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "JPID_sum_pt75_PreFilter" || cutDefinition == "JPID_sum1_pt75_sec_kV0_PreFilter"){
    AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01_PreFilter);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2_PreFilter);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "JPID_sum_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }


  // else if (cutDefinition == "JPID_sum1_pt75_sec_kV0"){
  //   AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01);
  //   // AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }


  /////////////////////////////////////////////////////////
  //             secondary track cut settings            //
  /////////////////////////////////////////////////////////
  else if (cutDefinition == "noPID_V0_PreFilter" || cutDefinition == "noPID_V0_standard" || cutDefinition == "track_V0_PreFilter"){
    // AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetPIDAna(LMEECutLib::kNoKinPIDCuts);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  // else if (cutDefinition == "trackkV0"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   // AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kV0track);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }

  else if (cutDefinition == "track_V0_standard"){
    AnaCut.SetPIDAna(LMEECutLib::kPID_V0_TPC_Pt75);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }


  // else if (cutDefinition == "JPID_sum_pt75_secondary"){
  //   AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01);
  //   // AnaCut.SetTrackSelectionAna(LMEECutLib::kV0);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1_secondary);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  // else if (cutDefinition == "JPID_sum_pt200_secondary"){
  //   AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01_pt200);
  //   // AnaCut.SetTrackSelectionAna(LMEECutLib::kV0);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1_secondary);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }

  /////////////////////////////////////////////////////////
  //              primary pair cut settings              //
  /////////////////////////////////////////////////////////
  else if (cutDefinition == "pairJPID_sum_pt20" || cutDefinition == "pairJPID_sum1_pt20_sec_kV0"){
    // AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt75);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  // else if (cutDefinition == "pairJPID_sum1_pt20_sec_kV0"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }

  /////////////////////////////////////////////////////////
  //             secondary pair cut settings             //
  /////////////////////////////////////////////////////////
  else if (cutDefinition == "pairkV0"){
    // AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt75);
    // AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1_secondary);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0pair);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "pairkV0_PreFilter"){
    // AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt75);
    // AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kV0OnTheFly);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1_secondary);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0pair_PreFilter);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  // else if (cutDefinition == "kV0OnlyCosOpenAngle"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyCos);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0OnlyChi2NDF"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyChi2NDF);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0OnlyLegDist"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyLegDist);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0OnlyR"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyR);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0OnlyPsiPair"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyPsiPair);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0OnlyM"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyM);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0OnlyArmPt"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyArmPt);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0OnlyArmAlpha"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyArmAlpha);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  //
  // else if (cutDefinition == "kV0wCosOpenAngle"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyCos);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0wChi2NDF"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0wCosChi2LegDist"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDist);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "wCosChi2LegDistR"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDistR);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "wCosChi2LegDistRPsiPair"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDistRPsiPair);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0wCosChi2LegDistRPsiPairM"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDistRPsiPairM);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0wCosChi2LegDistRPsiPairMArmPt"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDistRPsiPairMArmPt);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0wCosChi2LegDistRPsiPairMArmPtAlpha"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDistRPsiPairMArmPtAlpha);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0wPairM"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_wPairM);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0wPairMArmPt"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_wPairMArmPt);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }
  //
  // else if (cutDefinition == "kV0wPairMArmPtAlpha"){
  //   AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
  //   AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
  //   AnaCut.SetPairCutsAna(LMEECutLib::kV0_wPairMArmPtAlpha);
  //   AnaCut.SetCentrality(centrality);
  //   AnaCut.SetStandardCut();
  // }


  // Old Cut settings from CaKlein
  else if (cutDefinition == "cut1_pt75"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_1_pt75);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "onlyPIDcut1_noTrackCuts"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_1_pt75);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kNoTrackCuts);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut1_pt75_secondary"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_1_pt75);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kV0);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1_secondary);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  if (!isAOD) anaFilter->AddCuts( LMcutlib->GetESDTrackCutsAna(AnaCut) );

                                                                                // AliDielectronVarCuts *massCuts = new AliDielectronVarCuts("massCuts","massCuts");
                                                                                // massCuts->AddCut(AliDielectronVarManager::kM, 0.0,   0.250, kTRUE);
  // AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("gammaV0Cuts","gammaV0Cuts");
  // gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);  // kAll(default), kOffline or kOnTheFly
  // gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kOffline);  // kAll(default), kOffline or kOnTheFly
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.02, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // should increase purity...
  // gammaV0Cuts->SetExcludeTracks(kTRUE);

  anaFilter->AddCuts( LMcutlib->GetPIDCutsAna(AnaCut) );
  // anaFilter->AddCuts(gammaV0Cuts);
  // anaFilter->AddCuts(massCuts);
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
std::vector<bool> AddSinglePrimaryLegMCSignal(AliAnalysisTaskEtaReconstruction* task){
  AliDielectronSignalMC anyPartFinalState("anyPartFinalState","anyPartFinalState");
  anyPartFinalState.SetLegPDGs(0,1);//dummy second leg (never MCtrue)\n"
  anyPartFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  anyPartFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);

  AliDielectronSignalMC eleFinalState("eleFinalState","eleFinalState");
  eleFinalState.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);

  AliDielectronSignalMC eleFinalStateFromSameMotherMeson("eleFinalStateFromSameMotherMeson","eleFinalStateFromSameMotherMeson");
  eleFinalStateFromSameMotherMeson.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalStateFromSameMotherMeson.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromSameMotherMeson.SetMotherPDGs(600, 600); // open charm mesons and baryons together
  eleFinalStateFromSameMotherMeson.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //
  AliDielectronSignalMC eleFinalStateFromD("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD.SetMotherPDGs(402, 402); // open charm mesons and baryons together
  eleFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC eleFinalStateFromB("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB.SetMotherPDGs(502, 502); // open charm mesons and baryons together
  eleFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC eleFinalStateFromEta("eleFinalStateFromEta","eleFinalStateFromEta");
  eleFinalStateFromEta.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalStateFromEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromEta.SetMotherPDGs(221, 221); // open charm mesons and baryons together
  eleFinalStateFromEta.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC eleFinalStateFromPion("eleFinalStateFromPion","eleFinalStateFromPion");
  eleFinalStateFromPion.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalStateFromPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromPion.SetMotherPDGs(111, 111); // open charm mesons and baryons together
  eleFinalStateFromPion.SetCheckBothChargesMothers(kTRUE,kTRUE);

  AliDielectronSignalMC PhotonFinalState("PhotonFinalState","PhotonFinalState");
  PhotonFinalState.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);

  AliDielectronSignalMC PhotonFinalStateFromSameMotherMeson("PhotonFinalStateFromSameMotherMeson","PhotonFinalStateFromSameMotherMeson");
  PhotonFinalStateFromSameMotherMeson.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalStateFromSameMotherMeson.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalStateFromSameMotherMeson.SetMotherPDGs(600, 600); // open charm mesons and baryons together
  PhotonFinalStateFromSameMotherMeson.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //
  AliDielectronSignalMC PhotonFinalStateFromD("PhotonFinalStateFromD","PhotonFinalStateFromD");
  PhotonFinalStateFromD.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  PhotonFinalStateFromD.SetMotherPDGs(402, 402); // open charm mesons and baryons together
  PhotonFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC PhotonFinalStateFromB("PhotonFinalStateFromB","PhotonFinalStateFromB");
  PhotonFinalStateFromB.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  PhotonFinalStateFromB.SetMotherPDGs(502, 502); // open charm mesons and baryons together
  PhotonFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC PhotonFinalStateFromEta("PhotonFinalStateFromEta","PhotonFinalStateFromEta");
  PhotonFinalStateFromEta.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalStateFromEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalStateFromEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  PhotonFinalStateFromEta.SetMotherPDGs(221, 221); // open charm mesons and baryons together
  PhotonFinalStateFromEta.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC PhotonFinalStateFromPion("PhotonFinalStateFromPion","PhotonFinalStateFromPion");
  PhotonFinalStateFromPion.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalStateFromPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalStateFromPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  PhotonFinalStateFromPion.SetMotherPDGs(111, 111); // open charm mesons and baryons together
  PhotonFinalStateFromPion.SetCheckBothChargesMothers(kTRUE,kTRUE);




  if(UseMCDataSig) task->AddSinglePrimaryLegMCSignal(anyPartFinalState);

  task->AddSinglePrimaryLegMCSignal(eleFinalState);
  task->AddSinglePrimaryLegMCSignal(eleFinalStateFromPion);
  // task->AddSinglePrimaryLegMCSignal(eleFinalStateFromD);
  // task->AddSinglePrimaryLegMCSignal(eleFinalStateFromB);
  // task->AddSinglePrimaryLegMCSignal(eleFinalStateFromSameMotherMeson);
  task->AddSinglePrimaryLegMCSignal(eleFinalStateFromEta);
  // //
  // task->AddSinglePrimaryLegMCSignal(PhotonFinalState);
  // task->AddSinglePrimaryLegMCSignal(PhotonFinalStateFromPion);
  // task->AddSinglePrimaryLegMCSignal(PhotonFinalStateFromD);
  // task->AddSinglePrimaryLegMCSignal(PhotonFinalStateFromB);
  // task->AddSinglePrimaryLegMCSignal(PhotonFinalStateFromSameMotherMeson);
  // task->AddSinglePrimaryLegMCSignal(PhotonFinalStateFromEta);

  // this is used to get electrons from charmed mesons in a environment where GEANT is doing the decay of D mesons, like in LHC18b5a
  // ordering is according to MCSignals of single legs
 std::vector<bool> DielectronsPairNotFromSameMother;
 DielectronsPairNotFromSameMother.push_back(true);
 DielectronsPairNotFromSameMother.push_back(true);
 DielectronsPairNotFromSameMother.push_back(true);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 if(UseMCDataSig) DielectronsPairNotFromSameMother.push_back(false);
 return DielectronsPairNotFromSameMother;
}


std::vector<bool> AddSingleSecondaryLegMCSignal(AliAnalysisTaskEtaReconstruction* task){

  AliDielectronSignalMC eleSecondary("eleSecondary","eleSecondary");
  eleSecondary.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleSecondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleSecondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);

  AliDielectronSignalMC eleDontCare("eleDontCare","eleDontCare");
  eleDontCare.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleDontCare.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleDontCare.SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);

  AliDielectronSignalMC eleSecondaryFromPhoton("eleSecondaryFromPhoton","eleSecondaryFromPhoton");
  eleSecondaryFromPhoton.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleSecondaryFromPhoton.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleSecondaryFromPhoton.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  eleSecondaryFromPhoton.SetMotherPDGs(22, 1); // open charm mesons and baryons together
  eleSecondaryFromPhoton.SetCheckBothChargesMothers(kTRUE,kTRUE);
  // eleSecondaryFromPhoton.SetMotherSource(AliDielectronSignalMC::kSecondary);

  AliDielectronSignalMC eleSecondaryFromPhotonFromEta("eleSecondaryFromPhotonFromEta","eleSecondaryFromPhotonFromEta");
  eleSecondaryFromPhotonFromEta.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleSecondaryFromPhotonFromEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleSecondaryFromPhotonFromEta.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  eleSecondaryFromPhotonFromEta.SetMotherPDGs(22, 1); // open charm mesons and baryons together
  eleSecondaryFromPhotonFromEta.SetCheckBothChargesMothers(kTRUE,kTRUE);
  // eleSecondaryFromPhotonFromEta.SetMotherSource(AliDielectronSignalMC::kSecondary);
  eleSecondaryFromPhotonFromEta.SetGrandMotherPDGs(221, 1); // open charm mesons and baryons together
  eleSecondaryFromPhotonFromEta.SetCheckBothChargesMothers(kTRUE,kTRUE);


  AliDielectronSignalMC PhotonSecondary("PhotonSecondary","PhotonSecondary");
  PhotonSecondary.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonSecondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonSecondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);

  AliDielectronSignalMC PhotonDontCare("PhotonDontCare","PhotonDontCare");
  eleDontCare.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  eleDontCare.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleDontCare.SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);

  AliDielectronSignalMC anyPartSecondary("anyPartSecondary","anyPartSecondary");
  anyPartSecondary.SetLegPDGs(0,1);
  anyPartSecondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
  anyPartSecondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);


  if(UseMCDataSig) task->AddSingleSecondaryLegMCSignal(anyPartSecondary);
  //
  // task->AddSingleSecondaryLegMCSignal(eleDontCare);
  // task->AddSingleSecondaryLegMCSignal(eleSecondary);
  task->AddSingleSecondaryLegMCSignal(eleSecondaryFromPhoton);
  // task->AddSingleSecondaryLegMCSignal(eleSecondaryFromPhotonFromEta);
  // task->AddSingleSecondaryLegMCSignal(PhotonSecondary);
  // task->AddSingleSecondaryLegMCSignal(PhotonDontCare);

  // this is used to get electrons from charmed mesons in a environment where GEANT is doing the decay of D mesons, like in LHC18b5a
  // ordering is according to MCSignals of single legs
 std::vector<bool> DielectronsPairNotFromSameMother;
 DielectronsPairNotFromSameMother.push_back(true);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 if(UseMCDataSig) DielectronsPairNotFromSameMother.push_back(false);
 return DielectronsPairNotFromSameMother;
}

// #########################################################
// #########################################################
void AddPrimaryPairMCSignal(AliAnalysisTaskEtaReconstruction* task){

    AliDielectronSignalMC elePair_sameMother_finalstate("elePair_sameMother_finalstate","elePair_sameMother_finalstate");
    elePair_sameMother_finalstate.SetLegPDGs(11,-11);
    elePair_sameMother_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_sameMother_finalstate.SetMothersRelation(AliDielectronSignalMC::kSame);

    AliDielectronSignalMC elePair_DifferentMother_finalstate("elePair_DifferentMother_finalstate","elePair_DifferentMother_finalstate");
    elePair_DifferentMother_finalstate.SetLegPDGs(11,-11);
    elePair_DifferentMother_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_DifferentMother_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_DifferentMother_finalstate.SetMothersRelation(AliDielectronSignalMC::kDifferent);

    AliDielectronSignalMC elePair_UndefinedMother_finalstate("elePair_UndefinedMother_finalstate","elePair_UndefinedMother_finalstate");
    elePair_UndefinedMother_finalstate.SetLegPDGs(11,-11);
    elePair_UndefinedMother_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_UndefinedMother_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_UndefinedMother_finalstate.SetMothersRelation(AliDielectronSignalMC::kUndefined);

    AliDielectronSignalMC elePair_sameMother_photon_finalstate("elePair_sameMother_photon_finalstate","elePair_sameMother_photon_finalstate");
    elePair_sameMother_photon_finalstate.SetLegPDGs(11,-11);
    elePair_sameMother_photon_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_photon_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_sameMother_photon_finalstate.SetMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_photon_finalstate.SetMotherPDGs(22,22);

    AliDielectronSignalMC elePair_sameMother_eta_finalstate("elePair_sameMother_eta_finalstate","elePair_sameMother_eta_finalstate");
    elePair_sameMother_eta_finalstate.SetLegPDGs(11,-11);
    elePair_sameMother_eta_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_eta_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_sameMother_eta_finalstate.SetMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_eta_finalstate.SetMotherPDGs(221,221);

    AliDielectronSignalMC elePair_DifferentMother_eta_finalstate("elePair_DifferentMother_eta_finalstate","elePair_DifferentMother_eta_finalstate");
    elePair_DifferentMother_eta_finalstate.SetLegPDGs(11,-11);
    elePair_DifferentMother_eta_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_DifferentMother_eta_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_DifferentMother_eta_finalstate.SetMothersRelation(AliDielectronSignalMC::kDifferent);
    elePair_DifferentMother_eta_finalstate.SetMotherPDGs(221,221);

    AliDielectronSignalMC elePair_UndefinedMother_eta_finalstate("elePair_UndefinedMother_eta_finalstate","elePair_UndefinedMother_eta_finalstate");
    elePair_UndefinedMother_eta_finalstate.SetLegPDGs(11,-11);
    elePair_UndefinedMother_eta_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_UndefinedMother_eta_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_UndefinedMother_eta_finalstate.SetMothersRelation(AliDielectronSignalMC::kUndefined);
    elePair_UndefinedMother_eta_finalstate.SetMotherPDGs(221,221);

    AliDielectronSignalMC elePair_sameMother_pion_finalstate("elePair_sameMother_pion_finalstate","elePair_sameMother_pion_finalstate");
    elePair_sameMother_pion_finalstate.SetLegPDGs(11,-11);
    elePair_sameMother_pion_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_pion_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_sameMother_pion_finalstate.SetMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_pion_finalstate.SetMotherPDGs(111,111); //

    AliDielectronSignalMC elePair_DifferentMother_pion_finalstate("elePair_DifferentMother_pion_finalstate","elePair_DifferentMother_pion_finalstate");
    elePair_DifferentMother_pion_finalstate.SetLegPDGs(11,-11);
    elePair_DifferentMother_pion_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_DifferentMother_pion_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_DifferentMother_pion_finalstate.SetMothersRelation(AliDielectronSignalMC::kDifferent);
    elePair_DifferentMother_pion_finalstate.SetMotherPDGs(111,111);

    AliDielectronSignalMC elePair_UndefinedMother_pion_finalstate("elePair_UndefinedMother_pion_finalstate","elePair_UndefinedMother_pion_finalstate");
    elePair_UndefinedMother_pion_finalstate.SetLegPDGs(11,-11);
    elePair_UndefinedMother_pion_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_UndefinedMother_pion_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_UndefinedMother_pion_finalstate.SetMothersRelation(AliDielectronSignalMC::kUndefined);
    elePair_UndefinedMother_pion_finalstate.SetMotherPDGs(111,111);

    AliDielectronSignalMC elePair_sameMother_CharmedMesonsWithSameMother("CharmedMesonsWithSameMother","CharmedMesonsWithSameMother");
    elePair_sameMother_CharmedMesonsWithSameMother.SetLegPDGs(11,-11);
    elePair_sameMother_CharmedMesonsWithSameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_CharmedMesonsWithSameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_sameMother_CharmedMesonsWithSameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_CharmedMesonsWithSameMother.SetMotherPDGs(402, 402); //

    AliDielectronSignalMC elePair_sameMother_BeautyMesonsWithSameMother("BeautyMesonsWithSameMother","BeautyMesonsWithSameMother");
    elePair_sameMother_BeautyMesonsWithSameMother.SetLegPDGs(11,-11);
    elePair_sameMother_BeautyMesonsWithSameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_BeautyMesonsWithSameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    elePair_sameMother_BeautyMesonsWithSameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_BeautyMesonsWithSameMother.SetMotherPDGs(502, 502); //

    AliDielectronSignalMC anyPair_anyPart_UndefinedMother_finalstate("anyPair_anyPart_UndefinedMother_finalstate","anyPair_anyPart_UndefinedMother_finalstate");
    anyPair_anyPart_UndefinedMother_finalstate.SetLegPDGs(0,0);
    anyPair_anyPart_UndefinedMother_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    anyPair_anyPart_UndefinedMother_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    anyPair_anyPart_UndefinedMother_finalstate.SetMothersRelation(AliDielectronSignalMC::kUndefined);


    if(UseMCDataSig) task->AddPrimaryPairMCSignal(anyPair_anyPart_UndefinedMother_finalstate);

    task->AddPrimaryPairMCSignal(elePair_sameMother_finalstate);
    task->AddPrimaryPairMCSignal(elePair_DifferentMother_finalstate);
    task->AddPrimaryPairMCSignal(elePair_UndefinedMother_finalstate);
    task->AddPrimaryPairMCSignal(elePair_sameMother_photon_finalstate);
    task->AddPrimaryPairMCSignal(elePair_sameMother_eta_finalstate);
    task->AddPrimaryPairMCSignal(elePair_DifferentMother_eta_finalstate);
    task->AddPrimaryPairMCSignal(elePair_UndefinedMother_eta_finalstate);
    task->AddPrimaryPairMCSignal(elePair_sameMother_pion_finalstate);
    task->AddPrimaryPairMCSignal(elePair_DifferentMother_pion_finalstate);
    task->AddPrimaryPairMCSignal(elePair_UndefinedMother_pion_finalstate);
    // task->AddPrimaryPairMCSignal(elePair_sameMother_CharmedMesonsWithSameMother);
    // task->AddPrimaryPairMCSignal(elePair_sameMother_BeautyMesonsWithSameMother);
}

void AddSecondaryPairMCSignal(AliAnalysisTaskEtaReconstruction* task){
    AliDielectronSignalMC elePair_sameMother_secondary("elePair_sameMother_secondary","elePair_sameMother_secondary");
    elePair_sameMother_secondary.SetLegPDGs(11,-11);
    elePair_sameMother_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    elePair_sameMother_secondary.SetMothersRelation(AliDielectronSignalMC::kSame);

    AliDielectronSignalMC elePair_UndefinedMother_secondary("elePair_UndefinedMother_secondary","elePair_UndefinedMother_secondary");
    elePair_UndefinedMother_secondary.SetLegPDGs(11,-11);
    elePair_UndefinedMother_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_UndefinedMother_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    elePair_UndefinedMother_secondary.SetMothersRelation(AliDielectronSignalMC::kUndefined);

    AliDielectronSignalMC elePair_DifferentMother_secondary("elePair_DifferentMother_secondary","elePair_DifferentMother_secondary");
    elePair_DifferentMother_secondary.SetLegPDGs(11,-11);
    elePair_DifferentMother_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_DifferentMother_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    elePair_DifferentMother_secondary.SetMothersRelation(AliDielectronSignalMC::kDifferent);

    AliDielectronSignalMC elePair_sameMother_photon_secondary("elePair_sameMother_photon_secondary","elePair_sameMother_photon_secondary");
    elePair_sameMother_photon_secondary.SetLegPDGs(11,-11);
    elePair_sameMother_photon_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_photon_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    elePair_sameMother_photon_secondary.SetMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_photon_secondary.SetMotherPDGs(22,22);

    AliDielectronSignalMC elePair_sameMother_photon_secondaryfromMaterial("elePair_sameMother_photon_secondaryfromMaterial","elePair_sameMother_photon_secondaryfromMaterial");
    elePair_sameMother_photon_secondaryfromMaterial.SetLegPDGs(11,-11);
    elePair_sameMother_photon_secondaryfromMaterial.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_photon_secondaryfromMaterial.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
    //mother
    elePair_sameMother_photon_secondaryfromMaterial.SetMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_photon_secondaryfromMaterial.SetMotherPDGs(22,22);

    AliDielectronSignalMC elePair_sameMother_photon_secondaryfromWD("elePair_sameMother_photon_secondaryfromWD","elePair_sameMother_photon_secondaryfromWD");
    elePair_sameMother_photon_secondaryfromWD.SetLegPDGs(11,-11);
    elePair_sameMother_photon_secondaryfromWD.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_photon_secondaryfromWD.SetLegSources(AliDielectronSignalMC::kSecondaryFromWeakDecay, AliDielectronSignalMC::kSecondaryFromWeakDecay);
    //mother
    elePair_sameMother_photon_secondaryfromWD.SetMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_photon_secondaryfromWD.SetMotherPDGs(22,22);

    AliDielectronSignalMC elePair_UndefinedMother_photon_secondary("elePair_UndefinedMother_photon_secondary","elePair_UndefinedMother_photon_secondary");
    elePair_UndefinedMother_photon_secondary.SetLegPDGs(11,-11);
    elePair_UndefinedMother_photon_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_UndefinedMother_photon_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    elePair_UndefinedMother_photon_secondary.SetMothersRelation(AliDielectronSignalMC::kUndefined);
    elePair_UndefinedMother_photon_secondary.SetMotherPDGs(22,22);

    AliDielectronSignalMC elePair_DifferentMother_photon_secondary("elePair_DifferentMother_photon_secondary","elePair_DifferentMother_photon_secondary");
    elePair_DifferentMother_photon_secondary.SetLegPDGs(11,-11);
    elePair_DifferentMother_photon_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_DifferentMother_photon_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    elePair_DifferentMother_photon_secondary.SetMothersRelation(AliDielectronSignalMC::kDifferent);
    elePair_DifferentMother_photon_secondary.SetMotherPDGs(22,22);

    AliDielectronSignalMC elePair_sameMother_photon_secondary_pion("elePair_sameMother_photon_secondary_pion","elePair_sameMother_photon_secondary_pion");
    elePair_sameMother_photon_secondary_pion.SetLegPDGs(11,-11);
    elePair_sameMother_photon_secondary_pion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_photon_secondary_pion.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    elePair_sameMother_photon_secondary_pion.SetMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_photon_secondary_pion.SetMotherPDGs(22,22);
    //grand-mother
    // elePair_sameMother_photon_secondary_pion.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_photon_secondary_pion.SetGrandMotherPDGs(111,111);

    AliDielectronSignalMC elePair_sameMother_photon_secondary_eta("elePair_sameMother_photon_secondary_eta","elePair_sameMother_photon_secondary_eta");
    elePair_sameMother_photon_secondary_eta.SetLegPDGs(11,-11);
    elePair_sameMother_photon_secondary_eta.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_sameMother_photon_secondary_eta.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    elePair_sameMother_photon_secondary_eta.SetMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_photon_secondary_eta.SetMotherPDGs(22,22);
    //grand-mother
    // elePair_sameMother_photon_secondary_eta.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
    elePair_sameMother_photon_secondary_eta.SetGrandMotherPDGs(221,221);



    AliDielectronSignalMC elePair_conversion_secondary("elePair_conversion_secondary","elePair_conversion_secondary");
    elePair_conversion_secondary.SetLegPDGs(11,-11);
    elePair_conversion_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_conversion_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    elePair_conversion_secondary.SetMothersRelation(AliDielectronSignalMC::kSame);

    AliDielectronSignalMC elePair_random_secondary("elePair_random_secondary","elePair_random_secondary");
    elePair_random_secondary.SetLegPDGs(11,-11);
    elePair_random_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_random_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);


    AliDielectronSignalMC elePair_NotSameMother_secondary("elePair_NotSameMother_secondary","elePair_NotSameMother_secondary");
    elePair_NotSameMother_secondary.SetLegPDGs(11,-11);
    elePair_NotSameMother_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    elePair_NotSameMother_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    elePair_NotSameMother_secondary.SetMothersRelation(AliDielectronSignalMC::kDifferent);

    AliDielectronSignalMC anyPair_anyPart_random_secondary("anyPair_anyPart_random_secondary","anyPair_anyPart_random_secondary");
    anyPair_anyPart_random_secondary.SetLegPDGs(0,0);
    anyPair_anyPart_random_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    anyPair_anyPart_random_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);


    if(UseMCDataSig) task->AddSecondaryPairMCSignal(anyPair_anyPart_random_secondary);

    // task->AddSecondaryPairMCSignal(elePair_sameMother_secondary);
    // task->AddSecondaryPairMCSignal(elePair_UndefinedMother_secondary);
    // task->AddSecondaryPairMCSignal(elePair_DifferentMother_secondary);
task->AddSecondaryPairMCSignal(elePair_sameMother_photon_secondary);
    // task->AddSecondaryPairMCSignal(elePair_sameMother_photon_secondaryfromMaterial);
    // task->AddSecondaryPairMCSignal(elePair_sameMother_photon_secondaryfromWD);
    // task->AddSecondaryPairMCSignal(elePair_DifferentMother_photon_secondary);
    // task->AddSecondaryPairMCSignal(elePair_UndefinedMother_photon_secondary);
task->AddSecondaryPairMCSignal(elePair_sameMother_photon_secondary_pion);
task->AddSecondaryPairMCSignal(elePair_sameMother_photon_secondary_eta);
task->AddSecondaryPairMCSignal(elePair_conversion_secondary);
task->AddSecondaryPairMCSignal(elePair_random_secondary);
task->AddSecondaryPairMCSignal(elePair_NotSameMother_secondary);
}



void AddFourPairMCSignal_PrimSec(AliAnalysisTaskEtaReconstruction* task){
  // AliDielectronSignalMC for 4 particles
  // Seperate it into 2 pair MCSignals
  //
  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_FinalState("FourElePair1_FinalState","FourElePair1_FinalState");
        FourElePair1_FinalState.SetLegPDGs(11,-11);
        FourElePair1_FinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_FinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_FinalState.SetMothersRelation(AliDielectronSignalMC::kSame);


      // Second Pair
        AliDielectronSignalMC FourElePair2_Secondary_samePhoton("FourElePair2_Secondary_samePhoton","FourElePair2_Secondary_samePhoton");
        FourElePair2_Secondary_samePhoton.SetLegPDGs(11,-11);
        FourElePair2_Secondary_samePhoton.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_Secondary_samePhoton.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_Secondary_samePhoton.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_Secondary_samePhoton.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_Secondary_samePhoton.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_Secondary_samePhoton.SetMotherPDGs(22, 22); //


  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_FinalState_Dalitz("FourElePair1_FinalState_Dalitz","FourElePair1_FinalState_Dalitz");
        FourElePair1_FinalState_Dalitz.SetLegPDGs(11,-11);
        FourElePair1_FinalState_Dalitz.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_FinalState_Dalitz.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_FinalState_Dalitz.SetMothersRelation(AliDielectronSignalMC::kSame);
        // check if mother of first pair is grand mother of second pair or vice versa
        FourElePair1_FinalState_Dalitz.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)

      // Second Pair
        AliDielectronSignalMC FourElePair2_Secondary_samePhoton_Dalitz("FourElePair2_Secondary_samePhoton_Dalitz","FourElePair2_Secondary_samePhoton_Dalitz");
        FourElePair2_Secondary_samePhoton_Dalitz.SetLegPDGs(11,-11);
        FourElePair2_Secondary_samePhoton_Dalitz.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_Secondary_samePhoton_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_Secondary_samePhoton_Dalitz.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_Secondary_samePhoton_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_Secondary_samePhoton_Dalitz.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_Secondary_samePhoton_Dalitz.SetMotherPDGs(22, 22); //

        // check if mother of first pair is grand mother of second pair or vice versa
        FourElePair2_Secondary_samePhoton_Dalitz.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)

  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_samePion("FourElePair1_samePion","FourElePair1_samePion");
        FourElePair1_samePion.SetLegPDGs(11,-11);
        FourElePair1_samePion.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_samePion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_samePion.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_samePion.SetMotherPDGs(111, 111); //

      // Second Pair
        AliDielectronSignalMC FourElePair2_samePion("FourElePair2_samePion","FourElePair2_samePion");
        FourElePair2_samePion.SetLegPDGs(11,-11);
        FourElePair2_samePion.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_samePion.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_samePion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_samePion.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_samePion.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_samePion.SetMotherPDGs(22, 22); //
        // grand-mother
        // FourElePair2_samePion.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_samePion.SetGrandMotherPDGs(111, 111); //

  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_samePion_Dalitz("FourElePair1_samePion_Dalitz","FourElePair1_samePion_Dalitz");
        FourElePair1_samePion_Dalitz.SetLegPDGs(11,-11);
        FourElePair1_samePion_Dalitz.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_samePion_Dalitz.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_samePion_Dalitz.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_samePion_Dalitz.SetMotherPDGs(111, 111); //
        // check if mother of first pair is grand mother of second pair or vice versa
        FourElePair1_samePion_Dalitz.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)

      // Second Pair
        AliDielectronSignalMC FourElePair2_samePion_Dalitz("FourElePair2_samePion_Dalitz","FourElePair2_samePion_Dalitz");
        FourElePair2_samePion_Dalitz.SetLegPDGs(11,-11);
        FourElePair2_samePion_Dalitz.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_samePion_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_samePion_Dalitz.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_samePion_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_samePion_Dalitz.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_samePion_Dalitz.SetMotherPDGs(22, 22); //
        // grand-mother
        // FourElePair2_samePion_Dalitz.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_samePion_Dalitz.SetGrandMotherPDGs(111, 111); //
        // check if mother of first pair is grand mother of second pair or vice versa
        FourElePair2_samePion_Dalitz.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)

  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_sameEta("FourElePair1_sameEta","FourElePair1_sameEta");
        FourElePair1_sameEta.SetLegPDGs(11,-11);
        FourElePair1_sameEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_sameEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_sameEta.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_sameEta.SetMotherPDGs(221, 221); //

      // Second Pair
        AliDielectronSignalMC FourElePair2_sameEta("FourElePair2_sameEta","FourElePair2_sameEta");
        FourElePair2_sameEta.SetLegPDGs(11,-11);
        FourElePair2_sameEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_sameEta.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_sameEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_sameEta.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_sameEta.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_sameEta.SetMotherPDGs(22, 22); //
        // grand-mother
        // FourElePair2_sameEta.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_sameEta.SetGrandMotherPDGs(221, 221); //

  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_sameEta_Dalitz("FourElePair1_sameEta_Dalitz","FourElePair1_sameEta_Dalitz");
        FourElePair1_sameEta_Dalitz.SetLegPDGs(11,-11);
        FourElePair1_sameEta_Dalitz.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_sameEta_Dalitz.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_sameEta_Dalitz.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_sameEta_Dalitz.SetMotherPDGs(221, 221); //
        // check if mother of first pair is grand mother of second pair or vice versa
        FourElePair1_sameEta_Dalitz.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)

      // Second Pair
        AliDielectronSignalMC FourElePair2_sameEta_Dalitz("FourElePair2_sameEta_Dalitz","FourElePair2_sameEta_Dalitz");
        FourElePair2_sameEta_Dalitz.SetLegPDGs(11,-11);
        FourElePair2_sameEta_Dalitz.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_sameEta_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_sameEta_Dalitz.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_sameEta_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_sameEta_Dalitz.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_sameEta_Dalitz.SetMotherPDGs(22, 22); //
        // grand-mother
        // FourElePair2_sameEta_Dalitz.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_sameEta_Dalitz.SetGrandMotherPDGs(221, 221); //
        // check if mother of first pair is grand mother of second pair or vice versa
        FourElePair2_sameEta_Dalitz.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)


  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_MissMatchPionEta("FourElePair1_MissMatchPionEta","FourElePair1_MissMatchPionEta");
        FourElePair1_MissMatchPionEta.SetLegPDGs(11,-11);
        FourElePair1_MissMatchPionEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_MissMatchPionEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_MissMatchPionEta.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_MissMatchPionEta.SetMotherPDGs(111, 111); //
        // check if mother of first pair is grand mother of second pair or vice versa
        // FourElePair1_MissMatchPionEta.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)

      // Second Pair
        AliDielectronSignalMC FourElePair2_MissMatchPionEta("FourElePair2_MissMatchPionEta","FourElePair2_MissMatchPionEta");
        FourElePair2_MissMatchPionEta.SetLegPDGs(11,-11);
        FourElePair2_MissMatchPionEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_MissMatchPionEta.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_MissMatchPionEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_MissMatchPionEta.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_MissMatchPionEta.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_MissMatchPionEta.SetMotherPDGs(22, 22); //
        // grand-mother
        // FourElePair2_MissMatchPionEta.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_MissMatchPionEta.SetGrandMotherPDGs(221, 221); //
        // check if mother of first pair is grand mother of second pair or vice versa
        // FourElePair2_MissMatchPionEta.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)


  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_MissMatchEtaPion("FourElePair1_MissMatchEtaPion","FourElePair1_MissMatchEtaPion");
        FourElePair1_MissMatchEtaPion.SetLegPDGs(11,-11);
        FourElePair1_MissMatchEtaPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_MissMatchEtaPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_MissMatchEtaPion.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_MissMatchEtaPion.SetMotherPDGs(221, 221); //
        // check if mother of first pair is grand mother of second pair or vice versa
        // FourElePair1_MissMatchEtaPion.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)

      // Second Pair
        AliDielectronSignalMC FourElePair2_MissMatchEtaPion("FourElePair2_MissMatchEtaPion","FourElePair2_MissMatchEtaPion");
        FourElePair2_MissMatchEtaPion.SetLegPDGs(11,-11);
        FourElePair2_MissMatchEtaPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_MissMatchEtaPion.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_MissMatchEtaPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_MissMatchEtaPion.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_MissMatchEtaPion.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_MissMatchEtaPion.SetMotherPDGs(22, 22); //
        // grand-mother
        // FourElePair2_MissMatchEtaPion.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_MissMatchEtaPion.SetGrandMotherPDGs(111, 111); //
        // check if mother of first pair is grand mother of second pair or vice versa
        // FourElePair2_MissMatchEtaPion.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)


  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_FinalState_UndefinedMother("FourElePair1_FinalState_UndefinedMother","FourElePair1_FinalState_UndefinedMother");
        FourElePair1_FinalState_UndefinedMother.SetLegPDGs(11,-11);
        FourElePair1_FinalState_UndefinedMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_FinalState_UndefinedMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_FinalState_UndefinedMother.SetMothersRelation(AliDielectronSignalMC::kUndefined);


      // Second Pair
        AliDielectronSignalMC FourElePair2_Secondary_UndefinedMother("FourElePair2_Secondary_UndefinedMother","FourElePair2_Secondary_UndefinedMother");
        FourElePair2_Secondary_UndefinedMother.SetLegPDGs(11,-11);
        FourElePair2_Secondary_UndefinedMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_Secondary_UndefinedMother.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_Secondary_UndefinedMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_Secondary_UndefinedMother.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_Secondary_UndefinedMother.SetMothersRelation(AliDielectronSignalMC::kUndefined);
        // FourElePair2_Secondary_UndefinedMother.SetMotherPDGs(22, 22); //


  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourAnyPartPair1_FinalState_UndefinedMother("FourAnyPartPair1_FinalState_UndefinedMother","FourAnyPartPair1_FinalState_UndefinedMother");
        FourAnyPartPair1_FinalState_UndefinedMother.SetLegPDGs(0,0);
        FourAnyPartPair1_FinalState_UndefinedMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourAnyPartPair1_FinalState_UndefinedMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourAnyPartPair1_FinalState_UndefinedMother.SetMothersRelation(AliDielectronSignalMC::kUndefined);


      // Second Pair
        AliDielectronSignalMC FourAnyPartPair2_Secondary_UndefinedMother("FourAnyPartPair2_Secondary_UndefinedMother","FourAnyPartPair2_Secondary_UndefinedMother");
        FourAnyPartPair2_Secondary_UndefinedMother.SetLegPDGs(0,0);
        FourAnyPartPair2_Secondary_UndefinedMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourAnyPartPair2_Secondary_UndefinedMother.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // mother
        FourAnyPartPair2_Secondary_UndefinedMother.SetMothersRelation(AliDielectronSignalMC::kUndefined);


//________________________________________________________
    task->AddFourPairMCSignal_PrimSec(FourElePair1_FinalState);
    task->AddFourPairMCSignal_PrimSec(FourElePair2_Secondary_samePhoton);

    task->AddFourPairMCSignal_PrimSec(FourElePair1_FinalState_Dalitz);
    task->AddFourPairMCSignal_PrimSec(FourElePair2_Secondary_samePhoton_Dalitz);

    task->AddFourPairMCSignal_PrimSec(FourElePair1_samePion);
    task->AddFourPairMCSignal_PrimSec(FourElePair2_samePion);

    task->AddFourPairMCSignal_PrimSec(FourElePair1_samePion_Dalitz);
    task->AddFourPairMCSignal_PrimSec(FourElePair2_samePion_Dalitz);

    task->AddFourPairMCSignal_PrimSec(FourElePair1_sameEta);
    task->AddFourPairMCSignal_PrimSec(FourElePair2_sameEta);

    task->AddFourPairMCSignal_PrimSec(FourElePair1_sameEta_Dalitz);
    task->AddFourPairMCSignal_PrimSec(FourElePair2_sameEta_Dalitz);

    task->AddFourPairMCSignal_PrimSec(FourElePair1_MissMatchPionEta);
    task->AddFourPairMCSignal_PrimSec(FourElePair2_MissMatchPionEta);

    task->AddFourPairMCSignal_PrimSec(FourElePair1_MissMatchEtaPion);
    task->AddFourPairMCSignal_PrimSec(FourElePair2_MissMatchEtaPion);

    task->AddFourPairMCSignal_PrimSec(FourElePair1_FinalState_UndefinedMother);
    task->AddFourPairMCSignal_PrimSec(FourElePair2_Secondary_UndefinedMother);

    if(UseMCDataSig) task->AddFourPairMCSignal_PrimSec(FourAnyPartPair1_FinalState_UndefinedMother);
    if(UseMCDataSig) task->AddFourPairMCSignal_PrimSec(FourAnyPartPair2_Secondary_UndefinedMother);
}




void AddFourPairMCSignal_SecSec(AliAnalysisTaskEtaReconstruction* task){
  // AliDielectronSignalMC for 4 particles
  // Seperate it into 2 pair MCSignals
  //
  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_Secondary_samePhoton("FourElePair1_Secondary_samePhoton","FourElePair1_Secondary_samePhoton");
        FourElePair1_Secondary_samePhoton.SetLegPDGs(11,-11);
        FourElePair1_Secondary_samePhoton.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_Secondary_samePhoton.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // mother
        FourElePair1_Secondary_samePhoton.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_Secondary_samePhoton.SetMotherPDGs(22, 22); //


      // Second Pair
        AliDielectronSignalMC FourElePair2_Secondary_samePhoton("FourElePair2_Secondary_samePhoton","FourElePair2_Secondary_samePhoton");
        FourElePair2_Secondary_samePhoton.SetLegPDGs(11,-11);
        FourElePair2_Secondary_samePhoton.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_Secondary_samePhoton.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // mother
        FourElePair2_Secondary_samePhoton.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_Secondary_samePhoton.SetMotherPDGs(22, 22); //


  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_Secondary_samePhoton_sameEta("FourElePair1_Secondary_samePhoton_sameEta","FourElePair1_Secondary_samePhoton_sameEta");
        FourElePair1_Secondary_samePhoton_sameEta.SetLegPDGs(11,-11);
        FourElePair1_Secondary_samePhoton_sameEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_Secondary_samePhoton_sameEta.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // mother
        FourElePair1_Secondary_samePhoton_sameEta.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_Secondary_samePhoton_sameEta.SetMotherPDGs(22, 22); //
        // grand mother
        FourElePair1_Secondary_samePhoton_sameEta.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_Secondary_samePhoton_sameEta.SetGrandMotherPDGs(221, 221); //


      // Second Pair
        AliDielectronSignalMC FourElePair2_Secondary_samePhoton_sameEta("FourElePair2_Secondary_samePhoton_sameEta","FourElePair2_Secondary_samePhoton_sameEta");
        FourElePair2_Secondary_samePhoton_sameEta.SetLegPDGs(11,-11);
        FourElePair2_Secondary_samePhoton_sameEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_Secondary_samePhoton_sameEta.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // mother
        FourElePair2_Secondary_samePhoton_sameEta.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_Secondary_samePhoton_sameEta.SetMotherPDGs(22, 22); //
        // grand mother
        FourElePair2_Secondary_samePhoton_sameEta.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_Secondary_samePhoton_sameEta.SetGrandMotherPDGs(221, 221); //

  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourAnyPartPair1_Secondary_UndefinedMother("FourAnyPartPair1_Secondary_UndefinedMother","FourAnyPartPair1_Secondary_UndefinedMother");
        FourAnyPartPair1_Secondary_UndefinedMother.SetLegPDGs(0,0);
        FourAnyPartPair1_Secondary_UndefinedMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourAnyPartPair1_Secondary_UndefinedMother.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // mother
        FourAnyPartPair1_Secondary_UndefinedMother.SetMothersRelation(AliDielectronSignalMC::kUndefined);


      // Second Pair
        AliDielectronSignalMC FourAnyPartPair2_Secondary_UndefinedMother("FourAnyPartPair2_Secondary_UndefinedMother","FourAnyPartPair2_Secondary_UndefinedMother");
        FourAnyPartPair2_Secondary_UndefinedMother.SetLegPDGs(0,0);
        FourAnyPartPair2_Secondary_UndefinedMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourAnyPartPair2_Secondary_UndefinedMother.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // mother
        FourAnyPartPair2_Secondary_UndefinedMother.SetMothersRelation(AliDielectronSignalMC::kUndefined);

        //
        // //________________________________________________________
        //     // First Pair
        //       AliDielectronSignalMC FourPairMCSignalTest1("FourPairMCSignalTest1","FourPairMCSignalTest1");
        //       FourPairMCSignalTest1.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        //       FourPairMCSignalTest1.SetMothersRelation(AliDielectronSignalMC::kSame);
        //     // Second Pair
        //       AliDielectronSignalMC FourPairMCSignalTest2("FourPairMCSignalTest2","FourPairMCSignalTest2");
        //       FourPairMCSignalTest2.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        //       FourPairMCSignalTest2.SetMothersRelation(AliDielectronSignalMC::kSame);


// ________________________________________________________
    task->AddFourPairMCSignal_SecSec(FourElePair1_Secondary_samePhoton);
    task->AddFourPairMCSignal_SecSec(FourElePair2_Secondary_samePhoton);

    task->AddFourPairMCSignal_SecSec(FourElePair1_Secondary_samePhoton_sameEta);
    task->AddFourPairMCSignal_SecSec(FourElePair2_Secondary_samePhoton_sameEta);


    if(UseMCDataSig) task->AddFourPairMCSignal_SecSec(FourAnyPartPair1_Secondary_UndefinedMother);
    if(UseMCDataSig) task->AddFourPairMCSignal_SecSec(FourAnyPartPair2_Secondary_UndefinedMother);

    // task->AddFourPairMCSignal_SecSec(FourPairMCSignalTest1);
    // task->AddFourPairMCSignal_SecSec(FourPairMCSignalTest2);

}
