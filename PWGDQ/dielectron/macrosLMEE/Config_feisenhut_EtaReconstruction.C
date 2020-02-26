
// TString names=("kPbPb2015_Pt400_looseTOFif;kPbPb2015_Pt400_tightTOFreq;NewPID;kPbPb2015_Pt400_looseTOFif_0SITSCl;kPbPb2015_Pt400_tightTOFreq_0SITSCl;NewPID_0SITSCl;kPbPb2015_Pt400_looseTOFif_ITSMAP;kPbPb2015_Pt400_tightTOFreq_ITSMAP;NewPID_ITSMAP");
// TString names=("cut1_pt400;cut2_pt400;cut3_pt400;cut4_pt400;cut5_pt400;cut6_pt400;cut7_pt400;cut8_pt400;cut9_pt400;cut10_pt400;cut11_pt400;cut12_pt400;cut13_pt400;cut14_pt400;cut15_pt400;cut16_pt400;cut17_pt400;cut18_pt400;cut19_pt400;cut20_pt400");
// TString names=("noPID;cut5_pt75;kPbPb2015_Pt75_PID_cutoff_pion_kaon_proton;cut5_pt75_TOFreq");
// TString names=("cut1_pt75;cut2_pt75;cut3_pt75;cut4_pt75;cut5_pt75;cut6_pt75;cut7_pt75;cut8_pt75;cut9_pt75;cut10_pt75;cut11_pt75;cut12_pt75;cut13_pt75;cut14_pt75;cut15_pt75;cut16_pt75;cut17_pt75;cut18_pt75;cut19_pt75;cut20_pt75;cut21_pt75;cut22_pt75;cut23_pt75;cut24_pt75;cut25_pt75;cut26_pt75;cut27_pt75;cut28_pt75;cut29_pt75;cut30_pt75");
// TString names=("cut5_pt75;cut5_pt75_noPID_wSharedCluster;cut5_pt75_noPID_noSharedCluster");
// TString names=("cut1_pt75;cut2_pt75;cut3_pt75;cut4_pt75;cut5_pt75;cut6_pt75;cut7_pt75;cut8_pt75;cut9_pt75;cut10_pt75");

// +++++++++++++++++ Cuts for primary electrons +++++++++++++++++
// TString names_Prim_Cuts=("noPID");     // still has kin cuts (pt 75 MeV/c)
// TString names_Prim_Cuts=("noPID_pt200");  // still has kin cuts (pt 200 MeV/c)
// TString names_Prim_Cuts=("onlyJPID_sum_noTrackCuts");
// TString names_Prim_Cuts=("JPID_sum_pt75");
// TString names_Prim_Cuts=("JPID_sum_pt200");
// TString names_Prim_Cuts=("JPID_sum_pt75_secondary");
// TString names_Prim_Cuts=("JPID_sum1_pt75_sec_kV0");
// TString names_Prim_Cuts=("TESTCutskV0");

// TString names_Prim_Cuts=("noPID;JPID_sum1_pt75_sec_OnlyImpParXY;JPID_sum2_pt75_sec_OnlyImpParZ;JPID_sum3_pt75_sec_OnlyNclsTPC;JPID_sum4_pt75_sec_OnlyTPCchi2Cl;JPID_sum5_pt75_sec_OnlySharedCls;JPID_sum1_pt75_sec_kV0;JPID_sum_pt75");
// TString names_Prim_Cuts=("noPID;JPID_sum1_pt75_sec_wImpParXY;JPID_sum2_pt75_sec_wImpParXYZ;JPID_sum3_pt75_sec_wImpParXYZ_NclsTPC;JPID_sum4_pt75_sec_wImpParXYZ_NclsTPC_TPCchi2Cl;JPID_sum5_pt75_sec_wImpParXYZ_NclsTPC_TPCchi2Cl_SharedCls;JPID_sum_pt75");

// TString names_Prim_Cuts=("JPID_sum_pt75;JPID_sum2_pt75_sec_OnlyCosOpenAngle;JPID_sum3_pt75_sec_OnlyChi2NDF;JPID_sum4_pt75_sec_OnlyLegDist;JPID_sum5_pt75_sec_OnlyR;JPID_sum6_pt75_sec_OnlyPsiPair;JPID_sum7_pt75_sec_OnlyM;JPID_sum8_pt75_sec_OnlyArmPt;JPID_sum9_pt75_sec_OnlyArmAlpha;JPID_sum1_pt75_sec_kV0");
// TString names_Prim_Cuts=("JPID_sum_pt75;JPID_sum1_pt75_sec_wCos;JPID_sum2_pt75_sec_wCosChi2;JPID_sum3_pt75_sec_wCosChi2LegDist;JPID_sum4_pt75_sec_wCosChi2LegDistR;JPID_sum5_pt75_sec_wCosChi2LegDistRPsiPair;JPID_sum6_pt75_sec_wCosChi2LegDistRPsiPairM;JPID_sum7_pt75_sec_wCosChi2LegDistRPsiPairMArmPt;JPID_sum8_pt75_sec_wCosChi2LegDistRPsiPairMArmPtAlpha");
TString names_Prim_Cuts=("JPID_sum_pt75;JPID_sum2_pt75_sec_OnlyCosOpenAngle;JPID_sum7_pt75_sec_OnlyM;JPID_sum8_pt75_sec_OnlyArmPt;JPID_sum9_pt75_sec_OnlyArmAlpha;JPID_sum1_pt75_sec_kV0");
// TString names_Prim_Cuts=("JPID_sum_pt75;JPID_sum6_pt75_sec_wCosChi2LegDistRPsiPairM;JPID_sum7_pt75_sec_wCosChi2LegDistRPsiPairMArmPt;JPID_sum8_pt75_sec_wCosChi2LegDistRPsiPairMArmPtAlpha");

// +++++++++++++++++ Cuts for secondary electrons +++++++++++++++++
// TString names_Sec_Cuts=("noPID");      // still has kin cuts (pt 75 MeV/c)
// TString names_Sec_Cuts=("noPID_pt200");   // still has kin cuts (pt 200 MeV/c)
// TString names_Sec_Cuts=("onlyJPID_sum_noTrackCuts");
// TString names_Sec_Cuts=("JPID_sum_pt75");
// TString names_Sec_Cuts=("JPID_sum_pt75_secondary");
// TString names_Sec_Cuts=("JPID_sum_pt200_secondary");
// TString names_Sec_Cuts=("kV0");

// TString names_Sec_Cuts=("noPID;TESTCutsOnlyImpParXY;TESTCutsOnlyImpParZ;TESTCutsOnlyNclsTPC;TESTCutsOnlyTPCchi2Cl;TESTCutsOnlySharedCls;TESTCutskV0;JPID_sum_pt75_secondary");
// TString names_Sec_Cuts=("noPID;TESTCuts_wImpParXY;TESTCuts_wImpParXYZ;TESTCuts_wImpParXYZ_NclsTPC;TESTCuts_wImpParXYZ_NclsTPC_TPCchi2Cl;TESTCuts_wImpParXYZ_NclsTPC_TPCchi2Cl_SharedCls;JPID_sum_pt75_secondary");

// TString names_Sec_Cuts=("noPID;kV0OnlyCosOpenAngle;kV0OnlyChi2NDF;kV0OnlyLegDist;kV0OnlyR;kV0OnlyPsiPair;kV0OnlyM;kV0OnlyArmPt;kV0OnlyArmAlpha;kV0");
// TString names_Sec_Cuts=("noPID;kV0wCosOpenAngle;kV0wChi2NDF;kV0wCosChi2LegDist;wCosChi2LegDistR;wCosChi2LegDistRPsiPair;kV0wCosChi2LegDistRPsiPairM;kV0wCosChi2LegDistRPsiPairMArmPt;kV0wCosChi2LegDistRPsiPairMArmPtAlpha");
TString names_Sec_Cuts=("noPID;kV0OnlyCosOpenAngle;kV0OnlyM;kV0OnlyArmPt;kV0OnlyArmAlpha;kV0");
// TString names_Sec_Cuts=("noPID;kV0wPairM;kV0wPairMArmPt;kV0wPairMArmPtAlpha");


                                                                                // TString names_Sec_Cuts=("onlyPIDcut1_noTrackCuts");
                                                                                // TString names_Sec_Cuts=("cut1_pt75_secondary");
                                                                                // TString names_Sec_Cuts=("cut1_pt75");
                                                                                // TString names_Sec_Cuts=("noPID;cut1_pt75_secondary");


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

bool SetGeneratedSmearingHistos = false;

bool debug = false;

bool DoPairing      = true;
bool DoFourPairing  = true;
bool UsePreFilter   = true;
bool DoULSLS        = false;
bool DoMassCut      = false;

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


const double upperMassCutPrimaries = 0.547862;
const double lowerMassCutPrimaries = 0.1349766;
const double upperPreFilterMass = 0.2;
const double lowerPreFilterMass = 0.1;
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
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    // AnaCut.SetPIDAna(LMEECutLib::kNoPID_noKinCuts);
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

  // Primary Cut settings
  /////////////////////////////////////////////////////
  else if (cutDefinition == "JPID_sum_pt75"){
    AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2);
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

    // TEST Cut settings
  else if (cutDefinition == "JPID_sum1_pt75_sec_kV0"                   || cutDefinition == "JPID_sum1_pt75_sec_wCos"                               ||
           cutDefinition == "JPID_sum2_pt75_sec_OnlyCosOpenAngle"      || cutDefinition == "JPID_sum2_pt75_sec_wCosChi2"                           ||
           cutDefinition == "JPID_sum3_pt75_sec_OnlyChi2NDF"           || cutDefinition == "JPID_sum3_pt75_sec_wCosChi2LegDist"                    ||
           cutDefinition == "JPID_sum4_pt75_sec_OnlyLegDist"           || cutDefinition == "JPID_sum4_pt75_sec_wCosChi2LegDistR"                   ||
           cutDefinition == "JPID_sum5_pt75_sec_OnlyR"                 || cutDefinition == "JPID_sum5_pt75_sec_wCosChi2LegDistRPsiPair"            ||
           cutDefinition == "JPID_sum6_pt75_sec_OnlyPsiPair"           || cutDefinition == "JPID_sum6_pt75_sec_wCosChi2LegDistRPsiPairM"           ||
           cutDefinition == "JPID_sum7_pt75_sec_OnlyM"                 || cutDefinition == "JPID_sum7_pt75_sec_wCosChi2LegDistRPsiPairMArmPt"      ||
           cutDefinition == "JPID_sum8_pt75_sec_OnlyArmPt"             || cutDefinition == "JPID_sum8_pt75_sec_wCosChi2LegDistRPsiPairMArmPtAlpha" ||
           cutDefinition == "JPID_sum9_pt75_sec_OnlyArmAlpha"
          ){
    AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }


  // Secondary Cut settings
  /////////////////////////////////////////////////////
  else if (cutDefinition == "JPID_sum_pt75_secondary"){
    AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kV0);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1_secondary);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "JPID_sum_pt200_secondary"){
    AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01_pt200);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kV0);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1_secondary);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }



  else if (cutDefinition == "kV0"){
    // AnaCut.SetPIDAna(LMEECutLib::kPIDcut_TEST);
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    // AnaCut.SetPIDAna(LMEECutLib::kNoPID_noKinCuts);
    // AnaCut.SetPIDAna(LMEECutLib::kPID_Jeromian_01);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kV0);
    // AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1_secondary);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);

    AnaCut.SetPairCutsAna(LMEECutLib::kV0);

    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0OnlyCosOpenAngle"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyCos);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0OnlyChi2NDF"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyChi2NDF);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0OnlyLegDist"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyLegDist);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0OnlyR"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyR);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0OnlyPsiPair"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyPsiPair);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0OnlyM"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyM);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0OnlyArmPt"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyArmPt);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0OnlyArmAlpha"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyArmAlpha);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }


  else if (cutDefinition == "kV0wCosOpenAngle"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_onlyCos);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0wChi2NDF"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0wCosChi2LegDist"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDist);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "wCosChi2LegDistR"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDistR);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "wCosChi2LegDistRPsiPair"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDistRPsiPair);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0wCosChi2LegDistRPsiPairM"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDistRPsiPairM);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0wCosChi2LegDistRPsiPairMArmPt"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDistRPsiPairMArmPt);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0wCosChi2LegDistRPsiPairMArmPtAlpha"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_wCosChi2LegDistRPsiPairMArmPtAlpha);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0wPairM"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_wPairM);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0wPairMArmPt"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_wPairMArmPt);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kV0wPairMArmPtAlpha"){
    AnaCut.SetPIDAna(LMEECutLib::kNoPID_Pt20);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kDefaultNoTrackCuts);
    AnaCut.SetPairCutsAna(LMEECutLib::kV0_wPairMArmPtAlpha);
    AnaCut.SetCentrality(centrality);
    AnaCut.SetStandardCut();
  }


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
  /*AliDielectronSignalMC partFinalState("partFinalState","partFinalState");
  partFinalState.SetLegPDGs(0,1);//dummy second leg (never MCtrue)\n"
  partFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  partFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  */
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




  task->AddSinglePrimaryLegMCSignal(eleFinalState);
  // task->AddSinglePrimaryLegMCSignal(partFinalState);
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
 // DielectronsPairNotFromSameMother.push_back(false);
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
 // DielectronsPairNotFromSameMother.push_back(false);
 return DielectronsPairNotFromSameMother;
}

// #########################################################
// #########################################################
void AddPrimaryPairMCSignal(AliAnalysisTaskEtaReconstruction* task){

    AliDielectronSignalMC pair_sameMother_finalstate("pair_sameMother_finalstate","pair_sameMother_finalstate");
    pair_sameMother_finalstate.SetLegPDGs(11,-11);
    pair_sameMother_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_finalstate.SetMothersRelation(AliDielectronSignalMC::kSame);

    AliDielectronSignalMC pair_DifferentMother_finalstate("pair_DifferentMother_finalstate","pair_DifferentMother_finalstate");
    pair_DifferentMother_finalstate.SetLegPDGs(11,-11);
    pair_DifferentMother_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_DifferentMother_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_DifferentMother_finalstate.SetMothersRelation(AliDielectronSignalMC::kDifferent);

    AliDielectronSignalMC pair_UndefinedMother_finalstate("pair_UndefinedMother_finalstate","pair_UndefinedMother_finalstate");
    pair_UndefinedMother_finalstate.SetLegPDGs(11,-11);
    pair_UndefinedMother_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_UndefinedMother_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_UndefinedMother_finalstate.SetMothersRelation(AliDielectronSignalMC::kUndefined);

    AliDielectronSignalMC pair_sameMother_photon_finalstate("pair_sameMother_photon_finalstate","pair_sameMother_photon_finalstate");
    pair_sameMother_photon_finalstate.SetLegPDGs(11,-11);
    pair_sameMother_photon_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_photon_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_photon_finalstate.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_finalstate.SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons. should have no effect on final state ele.

    AliDielectronSignalMC pair_sameMother_eta_finalstate("pair_sameMother_eta_finalstate","pair_sameMother_eta_finalstate");
    pair_sameMother_eta_finalstate.SetLegPDGs(11,-11);
    pair_sameMother_eta_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_eta_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_eta_finalstate.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_eta_finalstate.SetMotherPDGs(221,221);

    AliDielectronSignalMC pair_DifferentMother_eta_finalstate("pair_DifferentMother_eta_finalstate","pair_DifferentMother_eta_finalstate");
    pair_DifferentMother_eta_finalstate.SetLegPDGs(11,-11);
    pair_DifferentMother_eta_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_DifferentMother_eta_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_DifferentMother_eta_finalstate.SetMothersRelation(AliDielectronSignalMC::kDifferent);
    pair_DifferentMother_eta_finalstate.SetMotherPDGs(221,221);

    AliDielectronSignalMC pair_UndefinedMother_eta_finalstate("pair_UndefinedMother_eta_finalstate","pair_UndefinedMother_eta_finalstate");
    pair_UndefinedMother_eta_finalstate.SetLegPDGs(11,-11);
    pair_UndefinedMother_eta_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_UndefinedMother_eta_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_UndefinedMother_eta_finalstate.SetMothersRelation(AliDielectronSignalMC::kUndefined);
    pair_UndefinedMother_eta_finalstate.SetMotherPDGs(221,221);

    AliDielectronSignalMC pair_sameMother_pion_finalstate("pair_sameMother_pion_finalstate","pair_sameMother_pion_finalstate");
    pair_sameMother_pion_finalstate.SetLegPDGs(11,-11);
    pair_sameMother_pion_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_pion_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_pion_finalstate.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_pion_finalstate.SetMotherPDGs(111,111); //

    AliDielectronSignalMC pair_DifferentMother_pion_finalstate("pair_DifferentMother_pion_finalstate","pair_DifferentMother_pion_finalstate");
    pair_DifferentMother_pion_finalstate.SetLegPDGs(11,-11);
    pair_DifferentMother_pion_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_DifferentMother_pion_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_DifferentMother_pion_finalstate.SetMothersRelation(AliDielectronSignalMC::kDifferent);
    pair_DifferentMother_pion_finalstate.SetMotherPDGs(111,111);

    AliDielectronSignalMC pair_UndefinedMother_pion_finalstate("pair_UndefinedMother_pion_finalstate","pair_UndefinedMother_pion_finalstate");
    pair_UndefinedMother_pion_finalstate.SetLegPDGs(11,-11);
    pair_UndefinedMother_pion_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_UndefinedMother_pion_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_UndefinedMother_pion_finalstate.SetMothersRelation(AliDielectronSignalMC::kUndefined);
    pair_UndefinedMother_pion_finalstate.SetMotherPDGs(111,111);

    AliDielectronSignalMC pair_sameMother_CharmedMesonsWithSameMother("CharmedMesonsWithSameMother","CharmedMesonsWithSameMother");
    pair_sameMother_CharmedMesonsWithSameMother.SetLegPDGs(11,-11);
    pair_sameMother_CharmedMesonsWithSameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_CharmedMesonsWithSameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_CharmedMesonsWithSameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_CharmedMesonsWithSameMother.SetMotherPDGs(402, 402); //

    AliDielectronSignalMC pair_sameMother_BeautyMesonsWithSameMother("BeautyMesonsWithSameMother","BeautyMesonsWithSameMother");
    pair_sameMother_BeautyMesonsWithSameMother.SetLegPDGs(11,-11);
    pair_sameMother_BeautyMesonsWithSameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_BeautyMesonsWithSameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_BeautyMesonsWithSameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_BeautyMesonsWithSameMother.SetMotherPDGs(502, 502); //


    task->AddPrimaryPairMCSignal(pair_sameMother_finalstate);
    task->AddPrimaryPairMCSignal(pair_DifferentMother_finalstate);
    task->AddPrimaryPairMCSignal(pair_UndefinedMother_finalstate);
    task->AddPrimaryPairMCSignal(pair_sameMother_photon_finalstate);
    task->AddPrimaryPairMCSignal(pair_sameMother_eta_finalstate);
    task->AddPrimaryPairMCSignal(pair_DifferentMother_eta_finalstate);
    task->AddPrimaryPairMCSignal(pair_UndefinedMother_eta_finalstate);
    // task->AddPrimaryPairMCSignal(pair_sameMother_pion_finalstate);
    // task->AddPrimaryPairMCSignal(pair_DifferentMother_pion_finalstate);
    // task->AddPrimaryPairMCSignal(pair_UndefinedMother_pion_finalstate);
    // task->AddPrimaryPairMCSignal(pair_sameMother_CharmedMesonsWithSameMother);
    // task->AddPrimaryPairMCSignal(pair_sameMother_BeautyMesonsWithSameMother);
}

void AddSecondaryPairMCSignal(AliAnalysisTaskEtaReconstruction* task){
    AliDielectronSignalMC pair_sameMother_secondary("pair_sameMother_secondary","pair_sameMother_secondary");
    pair_sameMother_secondary.SetLegPDGs(11,-11);
    pair_sameMother_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_sameMother_secondary.SetMothersRelation(AliDielectronSignalMC::kSame);

    AliDielectronSignalMC pair_UndefinedMother_secondary("pair_UndefinedMother_secondary","pair_UndefinedMother_secondary");
    pair_UndefinedMother_secondary.SetLegPDGs(11,-11);
    pair_UndefinedMother_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_UndefinedMother_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_UndefinedMother_secondary.SetMothersRelation(AliDielectronSignalMC::kUndefined);

    AliDielectronSignalMC pair_DifferentMother_secondary("pair_DifferentMother_secondary","pair_DifferentMother_secondary");
    pair_DifferentMother_secondary.SetLegPDGs(11,-11);
    pair_DifferentMother_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_DifferentMother_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_DifferentMother_secondary.SetMothersRelation(AliDielectronSignalMC::kDifferent);

    AliDielectronSignalMC pair_sameMother_photon_secondary("pair_sameMother_photon_secondary","pair_sameMother_photon_secondary");
    pair_sameMother_photon_secondary.SetLegPDGs(11,-11);
    pair_sameMother_photon_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_photon_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_sameMother_photon_secondary.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondary.SetMotherPDGs(22,22);

    AliDielectronSignalMC pair_sameMother_photon_secondaryfromMaterial("pair_sameMother_photon_secondaryfromMaterial","pair_sameMother_photon_secondaryfromMaterial");
    pair_sameMother_photon_secondaryfromMaterial.SetLegPDGs(11,-11);
    pair_sameMother_photon_secondaryfromMaterial.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_photon_secondaryfromMaterial.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
    //mother
    pair_sameMother_photon_secondaryfromMaterial.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondaryfromMaterial.SetMotherPDGs(22,22);

    AliDielectronSignalMC pair_sameMother_photon_secondaryfromWD("pair_sameMother_photon_secondaryfromWD","pair_sameMother_photon_secondaryfromWD");
    pair_sameMother_photon_secondaryfromWD.SetLegPDGs(11,-11);
    pair_sameMother_photon_secondaryfromWD.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_photon_secondaryfromWD.SetLegSources(AliDielectronSignalMC::kSecondaryFromWeakDecay, AliDielectronSignalMC::kSecondaryFromWeakDecay);
    //mother
    pair_sameMother_photon_secondaryfromWD.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondaryfromWD.SetMotherPDGs(22,22);

    AliDielectronSignalMC pair_UndefinedMother_photon_secondary("pair_UndefinedMother_photon_secondary","pair_UndefinedMother_photon_secondary");
    pair_UndefinedMother_photon_secondary.SetLegPDGs(11,-11);
    pair_UndefinedMother_photon_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_UndefinedMother_photon_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_UndefinedMother_photon_secondary.SetMothersRelation(AliDielectronSignalMC::kUndefined);
    pair_UndefinedMother_photon_secondary.SetMotherPDGs(22,22);

    AliDielectronSignalMC pair_DifferentMother_photon_secondary("pair_DifferentMother_photon_secondary","pair_DifferentMother_photon_secondary");
    pair_DifferentMother_photon_secondary.SetLegPDGs(11,-11);
    pair_DifferentMother_photon_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_DifferentMother_photon_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_DifferentMother_photon_secondary.SetMothersRelation(AliDielectronSignalMC::kDifferent);
    pair_DifferentMother_photon_secondary.SetMotherPDGs(22,22);

    AliDielectronSignalMC pair_sameMother_photon_secondary_pion("pair_sameMother_photon_secondary_pion","pair_sameMother_photon_secondary_pion");
    pair_sameMother_photon_secondary_pion.SetLegPDGs(11,-11);
    pair_sameMother_photon_secondary_pion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_photon_secondary_pion.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_sameMother_photon_secondary_pion.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondary_pion.SetMotherPDGs(22,22);
    //grand-mother
    pair_sameMother_photon_secondary_pion.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondary_pion.SetGrandMotherPDGs(111,111);

    AliDielectronSignalMC pair_sameMother_photon_secondary_eta("pair_sameMother_photon_secondary_eta","pair_sameMother_photon_secondary_eta");
    pair_sameMother_photon_secondary_eta.SetLegPDGs(11,-11);
    pair_sameMother_photon_secondary_eta.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_photon_secondary_eta.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_sameMother_photon_secondary_eta.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondary_eta.SetMotherPDGs(22,22);
    //grand-mother
    pair_sameMother_photon_secondary_eta.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondary_eta.SetGrandMotherPDGs(221,221);



    AliDielectronSignalMC pair_conversion_secondary("pair_conversion_secondary","pair_conversion_secondary");
    pair_conversion_secondary.SetLegPDGs(11,-11);
    pair_conversion_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_conversion_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_conversion_secondary.SetMothersRelation(AliDielectronSignalMC::kSame);

    AliDielectronSignalMC pair_random_secondary("pair_random_secondary","pair_random_secondary");
    pair_random_secondary.SetLegPDGs(11,-11);
    pair_random_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_random_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);


    AliDielectronSignalMC pair_NotSameMother_secondary("pair_NotSameMother_secondary","pair_NotSameMother_secondary");
    pair_NotSameMother_secondary.SetLegPDGs(11,-11);
    pair_NotSameMother_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_NotSameMother_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_NotSameMother_secondary.SetMothersRelation(AliDielectronSignalMC::kDifferent);

    // task->AddSecondaryPairMCSignal(pair_sameMother_secondary);
    // task->AddSecondaryPairMCSignal(pair_UndefinedMother_secondary);
    // task->AddSecondaryPairMCSignal(pair_DifferentMother_secondary);
task->AddSecondaryPairMCSignal(pair_sameMother_photon_secondary);
    // task->AddSecondaryPairMCSignal(pair_sameMother_photon_secondaryfromMaterial);
    // task->AddSecondaryPairMCSignal(pair_sameMother_photon_secondaryfromWD);
    // task->AddSecondaryPairMCSignal(pair_DifferentMother_photon_secondary);
    // task->AddSecondaryPairMCSignal(pair_UndefinedMother_photon_secondary);
    task->AddSecondaryPairMCSignal(pair_sameMother_photon_secondary_pion);
task->AddSecondaryPairMCSignal(pair_sameMother_photon_secondary_eta);
task->AddSecondaryPairMCSignal(pair_conversion_secondary);
task->AddSecondaryPairMCSignal(pair_random_secondary);
task->AddSecondaryPairMCSignal(pair_NotSameMother_secondary);
}

void AddFourPairMCSignal(AliAnalysisTaskEtaReconstruction* task){
  //#########################################
  //########### Unlike Signe Signals ########
  //#########################################
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
        AliDielectronSignalMC FourElePair2_Secondary_from_Photon("FourElePair2_Secondary_from_Photon","FourElePair2_Secondary_from_Photon");
        FourElePair2_Secondary_from_Photon.SetLegPDGs(11,-11);
        FourElePair2_Secondary_from_Photon.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_Secondary_from_Photon.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_Secondary_from_Photon.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_Secondary_from_Photon.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_Secondary_from_Photon.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_Secondary_from_Photon.SetMotherPDGs(22, 22); //


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
        AliDielectronSignalMC FourElePair2_Secondary_from_Photon_Dalitz("FourElePair2_Secondary_from_Photon_Dalitz","FourElePair2_Secondary_from_Photon_Dalitz");
        FourElePair2_Secondary_from_Photon_Dalitz.SetLegPDGs(11,-11);
        FourElePair2_Secondary_from_Photon_Dalitz.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_Secondary_from_Photon_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_Secondary_from_Photon_Dalitz.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_Secondary_from_Photon_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_Secondary_from_Photon_Dalitz.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_Secondary_from_Photon_Dalitz.SetMotherPDGs(22, 22); //

        // check if mother of first pair is grand mother of second pair or vice versa
        FourElePair2_Secondary_from_Photon_Dalitz.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)

  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_FromPion("FourElePair1_FromPion","FourElePair1_FromPion");
        FourElePair1_FromPion.SetLegPDGs(11,-11);
        FourElePair1_FromPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_FromPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_FromPion.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_FromPion.SetMotherPDGs(111, 111); //

      // Second Pair
        AliDielectronSignalMC FourElePair2_FromPion("FourElePair2_FromPion","FourElePair2_FromPion");
        FourElePair2_FromPion.SetLegPDGs(11,-11);
        FourElePair2_FromPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_FromPion.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_FromPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_FromPion.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_FromPion.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_FromPion.SetMotherPDGs(22, 22); //
        // grand-mother
        FourElePair2_FromPion.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_FromPion.SetGrandMotherPDGs(111, 111); //

  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_FromPion_Dalitz("FourElePair1_FromPion_Dalitz","FourElePair1_FromPion_Dalitz");
        FourElePair1_FromPion_Dalitz.SetLegPDGs(11,-11);
        FourElePair1_FromPion_Dalitz.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_FromPion_Dalitz.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_FromPion_Dalitz.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_FromPion_Dalitz.SetMotherPDGs(111, 111); //
        // check if mother of first pair is grand mother of second pair or vice versa
        FourElePair1_FromPion_Dalitz.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)

      // Second Pair
        AliDielectronSignalMC FourElePair2_FromPion_Dalitz("FourElePair2_FromPion_Dalitz","FourElePair2_FromPion_Dalitz");
        FourElePair2_FromPion_Dalitz.SetLegPDGs(11,-11);
        FourElePair2_FromPion_Dalitz.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_FromPion_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_FromPion_Dalitz.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_FromPion_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_FromPion_Dalitz.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_FromPion_Dalitz.SetMotherPDGs(22, 22); //
        // grand-mother
        FourElePair2_FromPion_Dalitz.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_FromPion_Dalitz.SetGrandMotherPDGs(111, 111); //
        // check if mother of first pair is grand mother of second pair or vice versa
        FourElePair2_FromPion_Dalitz.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)

  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_FromEta("FourElePair1_FromEta","FourElePair1_FromEta");
        FourElePair1_FromEta.SetLegPDGs(11,-11);
        FourElePair1_FromEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_FromEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_FromEta.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_FromEta.SetMotherPDGs(221, 221); //

      // Second Pair
        AliDielectronSignalMC FourElePair2_FromEta("FourElePair2_FromEta","FourElePair2_FromEta");
        FourElePair2_FromEta.SetLegPDGs(11,-11);
        FourElePair2_FromEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_FromEta.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_FromEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_FromEta.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_FromEta.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_FromEta.SetMotherPDGs(22, 22); //
        // grand-mother
        FourElePair2_FromEta.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_FromEta.SetGrandMotherPDGs(221, 221); //

  //________________________________________________________
      // First Pair
        AliDielectronSignalMC FourElePair1_FromEta_Dalitz("FourElePair1_FromEta_Dalitz","FourElePair1_FromEta_Dalitz");
        FourElePair1_FromEta_Dalitz.SetLegPDGs(11,-11);
        FourElePair1_FromEta_Dalitz.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair1_FromEta_Dalitz.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        FourElePair1_FromEta_Dalitz.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair1_FromEta_Dalitz.SetMotherPDGs(221, 221); //
        // check if mother of first pair is grand mother of second pair or vice versa
        FourElePair1_FromEta_Dalitz.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)

      // Second Pair
        AliDielectronSignalMC FourElePair2_FromEta_Dalitz("FourElePair2_FromEta_Dalitz","FourElePair2_FromEta_Dalitz");
        FourElePair2_FromEta_Dalitz.SetLegPDGs(11,-11);
        FourElePair2_FromEta_Dalitz.SetCheckBothChargesLegs(kTRUE,kTRUE);
        FourElePair2_FromEta_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_FromEta_Dalitz.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_FromEta_Dalitz.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        FourElePair2_FromEta_Dalitz.SetMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_FromEta_Dalitz.SetMotherPDGs(22, 22); //
        // grand-mother
        FourElePair2_FromEta_Dalitz.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        FourElePair2_FromEta_Dalitz.SetGrandMotherPDGs(221, 221); //
        // check if mother of first pair is grand mother of second pair or vice versa
        FourElePair2_FromEta_Dalitz.SetCheckMotherGrandmotherDiffPairRelation(kTRUE,kTRUE); // is assuming that both pairs using: SetMothersRelation(AliDielectronSignalMC::kSame)


//________________________________________________________
    task->AddFourPairMCSignal(FourElePair1_FinalState);
    task->AddFourPairMCSignal(FourElePair2_Secondary_from_Photon);

    task->AddFourPairMCSignal(FourElePair1_FinalState_Dalitz);
    task->AddFourPairMCSignal(FourElePair2_Secondary_from_Photon_Dalitz);

    task->AddFourPairMCSignal(FourElePair1_FromPion);
    task->AddFourPairMCSignal(FourElePair2_FromPion);

    task->AddFourPairMCSignal(FourElePair1_FromPion_Dalitz);
    task->AddFourPairMCSignal(FourElePair2_FromPion_Dalitz);

    task->AddFourPairMCSignal(FourElePair1_FromEta);
    task->AddFourPairMCSignal(FourElePair2_FromEta);

    task->AddFourPairMCSignal(FourElePair1_FromEta_Dalitz);
    task->AddFourPairMCSignal(FourElePair2_FromEta_Dalitz);
}
