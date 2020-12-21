//  Macro designed for use with the AliAnalysisTaskPIDBFDptDpt task.
//  Author: Jinjin(Au-Au) Pan & Claude Pruneau, Wayne State University
//
//   PbPb               10:     centralityMethod = 4 (V0),        trigger = kFALSE (AliVEvent::kMB).
//   PbPb_2015_kTRUE    15:     centralityMethod = 4 (V0),        trigger = kTRUE (AliVEvent::kINT7).
//   PbPb_2015_kFALSE   15:     centralityMethod = 4 (V0),        trigger = kTRUE (AliVEvent::kINT7).
//   pPb                13:     centralityMethod = 7 (V0A),       trigger = kTRUE (AliVEvent::kINT7).     
//   pp                 10:     centralityMethod = 3 (nTracks),   trigger = kFALSE (AliVEvent::kMB).
//   pp_V0A_kMB_kTRUE   10:     centralityMethod = 7 (V0A),       trigger = kFALSE (AliVEvent::kMB).
//   pp_V0A_kMB_kFALSE  10:     centralityMethod = 7 (V0A),       trigger = kFALSE (AliVEvent::kMB).
//   pp_V0_kMB_kTRUE    10:     centralityMethod = 4 (V0),        trigger = kFALSE (AliVEvent::kMB).
//   pp_V0_kMB_kFALSE   10:     centralityMethod = 4 (V0),        trigger = kFALSE (AliVEvent::kMB).
//   pp_V0C_kMB_kTRUE   10:     centralityMethod = 8 (V0C),       trigger = kFALSE (AliVEvent::kMB).
//   pp_V0C_kMB_kFALSE  10:     centralityMethod = 8 (V0C),       trigger = kFALSE (AliVEvent::kMB).
//   pp_V0A_kMB_Utils   10:     centralityMethod = 7 (V0A),       trigger = kFALSE (AliVEvent::kMB).
//   pp_V0C_kMB_Utils   10:     centralityMethod = 8 (V0C),       trigger = kFALSE (AliVEvent::kMB).
//   pp_V0_kMB_Utils    10:     centralityMethod = 4 (V0),        trigger = kFALSE (AliVEvent::kMB).
/////////////////////////////////////////////////////////////////////////////////

AliAnalysisTaskPIDBFDptDpt * AddTaskPIDBFDptDpt
(
 TString AnalysisDataType       = "RealData", // "RealData"; "MCAOD" for MC AOD truth; "MCAODreco"
 TString System                 = "PbPb",
 bool    pidparticle            =  1,   // 0: All Charged Particles;       1: PID particles
 bool    Use_PT_Cut             =  1,   // 0: Use_P_Cut ( only TOF lower boundary );       1: Use_PT_Cut
 int    useRapidity             =  1,   // 0: pseudo-rapadity      1: rapidity
 int    CentralityGroup         =  21,  // Diff Cent Groups dealing w/ memory limit & weight file 100M Alien limit
 int    singlesOnly             =  1,   // 0: full correlations    1: singles only
 int    useWeights              =  0,   // 0: no                   1: yes  
 int    chargeSet               =  1,   // 0: ++    1: +-    2: -+    3: --
 double zMax                    =  6.,  // set vertexZ cut   
 double vZwidth                 =  0.5, // zMin, zMax & vZwidth determine _nBins_vertexZ.
 int    trackFilterBit          =  1,   // PbPb10(Global=1;TPConly=128;Hybrid=272); pPb13(Global=?;TPConly=?;Hybrid=768); pp10(Global=1;TPConly=?; Hybrid=?)
 int    nClusterMin             =  70,
 double eta1Max                 =  0.8, // set y1max acturally if useRapidity==1
 double etaBinWidth             =  0.1, // set yBinWidth acturally if useRapidity==1
 double dcaZMax                 =  3.2,
 double dcaXYMax                =  2.4,
 int nCentrality                =  6,
 int particleID                 =  1,   // Pion=0, Kaon=1, Proton=2
 bool Use_CircularCutPID        =  1,   // 0: Not Use_CircularCutPID     1: Use_CircularCutPID TPC+TOF
 double nSigmaCut               =  2.0,
 double nSigmaCut_veto          =  3.0,
 double ElectronVetoCut         =  1.0,
 double ptMin                   =  0.2, // pt range lower limit cut ( also for pt histos )
 double ptTOFlowerMin           =  0.5, // boundary between TPC & TOF region
 double ptCUTupperMax           =  2.0, // pt range upper limit cut
 double ptMax                   =  3.0, // pt range upper limit for histos; NOT pt cut!!!
 double ptWidthBin              =  0.1, // pt bin width in histos
 int nBinsPhi                   =  36,  // 36 is default value
 Bool_t NoResonances            = kTRUE, // only for MCAOD
 Bool_t NoElectron              = kTRUE, // only for MCAOD
 bool   PurePIDinMC             = 0,   // 0: MisID in MCAODreco;       1: No MisID in MCAODreco
 bool   PureNoWeakinMC          = 0,   // 0: No MisID but Secondaries from weak decays in MCAODreco;       1: No MisID and No Secondaries from weak decays in MCAODreco
 bool   PureNoWeakMaterialinMC  = 0,   // 0: No MisID and No Secondaries from weak decays but Secondaries from material in MCAODreco;       1: No MisID and No Secondaries from weak decays and material in MCAODreco
 bool   NoMisIDWeakMaterialInClosure = 0,   // 0: allow MisID and No Secondaries from weak decays and material in MC Closure test       1: No MisID and Secondaries from weak decays and material in MC Closure test
 double SharedFractionPairCut   = 0.2, // check track splitting
 const char* taskname           = "ChPM",
 TString inputHistogramFileName = "alien:///alice/cern.ch/user/j/jipan/TUNE_rHJ_2eCut_8vZ32_G162_4C4_NOwCut_08y16_36phi_02pt2_pi_Pos_S1S2/TUNE_rHJ_2eCut_8vZ32_G162_4C4_NOwCut_08y16_36phi_02pt2_pi_Pos_S1S2.root" )

{
  // Set Default Configuration of this analysis
  // ==========================================
  int    debugLevel             = 0;
  int    rejectPileup           = 1;
  int    rejectPairConversion   = 1;
  int    sameFilter             = 1;
  int    centralityMethod       = 4; 
  Bool_t trigger                = kFALSE;
  Bool_t remove_Tracks_T0       = 1;
  bool    useEventPlane         = 0;   // 0: No      1: Yes
  double  EventPlaneMin         = -3.1415927/6;
  double  EventPlaneMax         =  3.1415927/6;
  bool Use_AliHelperPID         =  0;   // 0: Not Use_AliHelperPID       1: Use_AliHelperPID
  int pidType                   =  2;  // kNSigmaTPC,kNSigmaTOF, kNSigmaTPCTOF // for AliHelperPID
  Bool_t requestTOFPID          =  1;  // for AliHelperPID
  Bool_t isMC                   =  0;  // for AliHelperPID
  double eta2Max                =  eta1Max; // set y2max acturally if useRapidity==1
  double eta1Min                = -eta1Max; // set y1min acturally if useRapidity==1
  double eta2Min                = -eta1Max; // set y2min acturally if useRapidity==1
  double dcaZMin                = -dcaZMax;
  double dcaXYMin               = -dcaXYMax;
  double zMin                   = -zMax;  // |zMin| should = zMax due to the design of the code
	
  if      ( System == "PbPb" )                { centralityMethod = 4; trigger = kFALSE; }
  else if ( System == "PbPb_2015_kTRUE" )     { centralityMethod = 4; trigger = kTRUE;  }
  else if ( System == "PbPb_2015_kFALSE" )    { centralityMethod = 4; trigger = kTRUE;  }
  else if ( System == "pPb" )                 { centralityMethod = 7; trigger = kTRUE;  }
  else if ( System == "pp" )                  { centralityMethod = 3; trigger = kFALSE; }
  else if ( System == "pp_V0A_kMB_kTRUE" )    { centralityMethod = 7; trigger = kFALSE; }
  else if ( System == "pp_V0A_kMB_kFALSE" )   { centralityMethod = 7; trigger = kFALSE; }
  else if ( System == "pp_V0C_kMB_kTRUE" )    { centralityMethod = 8; trigger = kFALSE; }
  else if ( System == "pp_V0C_kMB_kFALSE" )   { centralityMethod = 8; trigger = kFALSE; }
  else if ( System == "pp_V0_kMB_kTRUE" )     { centralityMethod = 4; trigger = kFALSE; }
  else if ( System == "pp_V0_kMB_kFALSE" )    { centralityMethod = 4; trigger = kFALSE; }
  else if ( System == "pp_V0A_kMB_Utils" )    { centralityMethod = 7; trigger = kFALSE; }
  else if ( System == "pp_V0C_kMB_Utils" )    { centralityMethod = 8; trigger = kFALSE; }
  else if ( System == "pp_V0_kMB_Utils" )     { centralityMethod = 4; trigger = kFALSE; }
  else    return 0;

  
  double minCentrality[10];
  double maxCentrality[10];

  if ( CentralityGroup == 1 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 30.;
      minCentrality[3] = 30.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 50.;
      minCentrality[5] = 50.;     maxCentrality[5]  = 60.;
      minCentrality[6] = 60.;     maxCentrality[6]  = 70.;
      minCentrality[7] = 70.;     maxCentrality[7]  = 80.; }
  else if ( CentralityGroup == 2 )
    { minCentrality[0] = 40.;     maxCentrality[0]  = 50.;
      minCentrality[1] = 50.;     maxCentrality[1]  = 60.;
      minCentrality[2] = 60.;     maxCentrality[2]  = 70.;
      minCentrality[3] = 70.;     maxCentrality[3]  = 80.; }
  else if ( CentralityGroup == 3 )
    { minCentrality[0] = 20.;     maxCentrality[0]  = 30.;
      minCentrality[1] = 30.;     maxCentrality[1]  = 40.; }
  else if ( CentralityGroup == 4 )
    { minCentrality[0] = 60.;     maxCentrality[0]  = 70.;
      minCentrality[1] = 70.;     maxCentrality[1]  = 80.; }
  else if ( CentralityGroup == 5 )
    { minCentrality[0] = 10.;     maxCentrality[0]  = 20.; }
  else if ( CentralityGroup == 6 )
    { minCentrality[0] = 30.;     maxCentrality[0]  = 40.; }
  else if ( CentralityGroup == 7 )
    { minCentrality[0] = 50.;     maxCentrality[0]  = 60.; }
  else if ( CentralityGroup == 8 )
    { minCentrality[0] = 70.;     maxCentrality[0]  = 80.; }
  else if ( CentralityGroup == 9 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 20.;
      minCentrality[1] = 20.;     maxCentrality[1]  = 40.;
      minCentrality[2] = 40.;     maxCentrality[2]  = 60.;
      minCentrality[3] = 60.;     maxCentrality[3]  = 80.; }
  else if ( CentralityGroup == 10 )
    { minCentrality[0] = 40.;     maxCentrality[0]  = 60.;
      minCentrality[1] = 60.;     maxCentrality[1]  = 80.; }
  else if ( CentralityGroup == 11 )
    { minCentrality[0] = 20.;     maxCentrality[0]  = 40.; }
  else if ( CentralityGroup == 12 )
    { minCentrality[0] = 60.;     maxCentrality[0]  = 80.; }
  else if ( CentralityGroup == 13 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 100.; }
  else if ( CentralityGroup == 14 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 80.; }
  else if ( CentralityGroup == 15 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 60.;
      minCentrality[7] = 60.;     maxCentrality[7]  = 70.;
      minCentrality[8] = 70.;     maxCentrality[8]  = 80.;  }
  else if ( CentralityGroup == 16 )
    { minCentrality[0] = 20.;     maxCentrality[0]  = 100.; }
  else if ( CentralityGroup == 17 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 100.; }
  else if ( CentralityGroup == 18 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 60.;
      minCentrality[7] = 60.;     maxCentrality[7]  = 70.;
      minCentrality[8] = 70.;     maxCentrality[8]  = 100.;  }
  else if ( CentralityGroup == 19 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 60.;
      minCentrality[7] = 60.;     maxCentrality[7]  = 70.;
      minCentrality[8] = 70.;     maxCentrality[8]  = 90.;  }
  else if ( CentralityGroup == 20 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 30.;     maxCentrality[1]  = 40.;
      minCentrality[2] = 70.;     maxCentrality[2]  = 80.; }
  else if ( CentralityGroup == 21 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 30.;
      minCentrality[3] = 30.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 50.;
      minCentrality[5] = 50.;     maxCentrality[5]  = 80.; }
  else if ( CentralityGroup == 22 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 20.;
      minCentrality[1] = 20.;     maxCentrality[1]  = 40.;
      minCentrality[2] = 40.;     maxCentrality[2]  = 80.; }
  else if ( CentralityGroup == 23 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 80.; }
  else if ( CentralityGroup == 24 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 60.;
      minCentrality[7] = 60.;     maxCentrality[7]  = 80.; }
  else if ( CentralityGroup == 25 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 90.; }
  else if ( CentralityGroup == 26 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 60.;
      minCentrality[7] = 60.;     maxCentrality[7]  = 90.; }
  else if ( CentralityGroup == 27 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 30.;
      minCentrality[3] = 30.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 50.;
      minCentrality[5] = 50.;     maxCentrality[5]  = 90.; }
  else if ( CentralityGroup == 28 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 40.;
      minCentrality[3] = 40.;     maxCentrality[3]  = 60.;
      minCentrality[4] = 60.;     maxCentrality[4]  = 80.; }
  else if ( CentralityGroup == 29 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 30.;     maxCentrality[1]  = 40.;
      minCentrality[2] = 70.;     maxCentrality[2]  = 90.; }
  else if ( CentralityGroup == 30 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 30.;
      minCentrality[3] = 30.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 50.;
      minCentrality[5] = 50.;     maxCentrality[5]  = 70.;
      minCentrality[6] = 70.;     maxCentrality[6]  = 100.; }
  else if ( CentralityGroup == 31 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 60.;
      minCentrality[5] = 60.;     maxCentrality[5]  = 80.;
      minCentrality[6] = 80.;     maxCentrality[6]  = 100.; }
  else if ( CentralityGroup == 32 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 70.;
      minCentrality[7] = 70.;     maxCentrality[7]  = 100.; }
  else if ( CentralityGroup == 33 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 20.;
      minCentrality[1] = 20.;     maxCentrality[1]  = 40.;
      minCentrality[2] = 40.;     maxCentrality[2]  = 90.; }
  else if ( CentralityGroup == 34 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 30.;
      minCentrality[3] = 30.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 60.;
      minCentrality[5] = 60.;     maxCentrality[5]  = 90.; }
  else if ( CentralityGroup == 35 )
    { minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 60.;
      minCentrality[6] = 60.;     maxCentrality[6]  = 90.; }
  else    return 0;
  
  double dedxMin                =  0.0;
  double dedxMax                =  20000.0;
  int    requestedCharge1       =  1; //default
  int    requestedCharge2       = -1; //default
    
  // Get the pointer to the existing analysis manager via the static access method.
  // ==============================================================================
  AliAnalysisManager *analysisManager = AliAnalysisManager::GetAnalysisManager();
    
  if (!analysisManager)
    {
      ::Error("AddTaskPIDBFDptDpt", "No analysis manager to connect to.");
      return NULL;
    }
    
  TString part1Name;
  TString part2Name;
  TString eventName;
  TString prefixName        = "Corr_";
  TString pileupRejecSuffix = "_PileupRejec";
  TString pairRejecSuffix   = "_PairRejec";
  TString calibSuffix       = "_calib";
  TString singlesOnlySuffix = "_SO";
  TString suffix;
    
  TString inputPath         = ".";
  TString outputPath        = ".";
  TString baseName;
  TString listName;
  TString taskName;
  //TString inputHistogramFileName;
  TString outputHistogramFileName;
    
  // Create the task and add subtask.
  // ===========================================================================
  int iTask = 0; // task counter
  AliAnalysisDataContainer *taskInputContainer;
  AliAnalysisDataContainer *taskOutputContainer;
  AliAnalysisTaskPIDBFDptDpt* task;
    
  for (int iCentrality=0; iCentrality < nCentrality; ++iCentrality)
    {
      switch ( chargeSet )
        {
	case 0: part1Name = "P_"; part2Name = "P_"; requestedCharge1 =  1; requestedCharge2 =  1; sameFilter = 1;   break;
	case 1: part1Name = "P_"; part2Name = "M_"; requestedCharge1 =  1; requestedCharge2 = -1; sameFilter = 0;   break;
	case 2: part1Name = "M_"; part2Name = "P_"; requestedCharge1 = -1; requestedCharge2 =  1; sameFilter = 0;   break;
	case 3: part1Name = "M_"; part2Name = "M_"; requestedCharge1 = -1; requestedCharge2 = -1; sameFilter = 1;   break;
        }
        
      part1Name += "eta";
      part1Name += int(1000*eta1Max);
      part1Name += "_";
      part1Name += int(1000*ptMin);
      part1Name += "pt";
      part1Name += int(1000*ptMax);
      part1Name += "_";
      part1Name += int(1000*dcaZMin);
      part1Name += "DCA";
      part1Name += int(1000*dcaZMax);
      part1Name += "_";
        
      part2Name += "eta";
      part2Name += int(1000*eta2Max);
      part2Name += "_";
      part2Name += int(1000*ptMin);
      part2Name += "pt";
      part2Name += int(1000*ptMax);
      part2Name += "_";
      part2Name += int(1000*dcaZMin);
      part2Name += "DCA";
      part2Name += int(1000*dcaZMax);
      part2Name += "_";
        
      eventName =  "";
      eventName += int(10.*minCentrality[iCentrality] );
      eventName += "Vo";
      eventName += int(10.*maxCentrality[iCentrality] );
        
      baseName     =   prefixName;
      baseName     +=  part1Name;
      baseName     +=  part2Name;
      baseName     +=  eventName;
      baseName     +=  "_";
      baseName     +=  particleID;
      listName     =   baseName;
      taskName     =   baseName;        
        
      outputHistogramFileName = baseName;
      //if (singlesOnly) outputHistogramFileName += singlesOnlySuffix;
      outputHistogramFileName += ".root";      
        
      TFile  * inputFile  = 0;
      TList  * histoList  = 0;
      TH3F   * weight_1   = 0;
      TH3F   * weight_2   = 0;
      if ( useWeights )
        {
	  TGrid::Connect("alien:");
	  inputFile = TFile::Open(inputHistogramFileName,"OLD");
	  if (!inputFile)
            {
	      //cout << "Requested file:" << inputHistogramFileName << " was not opened. ABORT." << endl;
	      return 0;
            }
	  TString nameHistoBase = "correction_";
	  TString nameHisto;
	  nameHistoBase += eventName;
	  if (requestedCharge1 == 1)
            {
	      nameHisto = nameHistoBase + "_p";
	      //cout << "Input Histogram named: " << nameHisto << endl;
	      weight_1 = (TH3F *) inputFile->Get(nameHisto);
            }
	  else
            {
	      nameHisto = nameHistoBase + "_m";
	      //cout << "Input Histogram named: " << nameHisto << endl;
	      weight_1 = (TH3F *) inputFile->Get(nameHisto);
            }
	  if (!weight_1)
            {
	      //cout << "Requested histogram 'correction_p/m' was not found. ABORT." << endl;
	      return 0;
            }
            
	  if (!sameFilter)
            {
	      weight_2 = 0;
	      if (requestedCharge2 == 1)
                {
		  nameHisto = nameHistoBase + "_p";
		  //cout << "Input Histogram named: " << nameHisto << endl;
		  weight_2 = (TH3F *) inputFile->Get(nameHisto);
                }
	      else
                {
		  nameHisto = nameHistoBase + "_m";
		  //cout << "Input Histogram named: " << nameHisto << endl;
		  weight_2 = (TH3F *) inputFile->Get(nameHisto);
                }
	      if (!weight_2)
                {
		  //cout << "Requested histogram 'correction_p/m' was not found. ABORT." << endl;
		  return 0;
                }
            }
        }
      
      task = new  AliAnalysisTaskPIDBFDptDpt(taskName);
      //configure my task
      task->SetDebugLevel(          debugLevel      );
      task->SetSameFilter(          sameFilter      );
      task->SetSinglesOnly(         singlesOnly     );
      task->SetPIDparticle(         pidparticle     );
      task->SetUse_pT_cut(          Use_PT_Cut      );
      task->SetIfContaminationInMC(   PurePIDinMC   );
      task->SetIfContaminationWeakInMC( PureNoWeakinMC );
      task->SetIfContaminationWeakMaterialInMC( PureNoWeakMaterialinMC );
      task->SetIfMisIDWeakMaterialInMCClosure( NoMisIDWeakMaterialInClosure );
      task->SetUseWeights(          useWeights      );
      task->SetUseRapidity(         useRapidity     );
      task->SetEventPlane(         useEventPlane     );
      task->SetEPmin(              EventPlaneMin     );
      task->SetEPmax(              EventPlaneMax     );
      task->SetRejectPileup(        rejectPileup    );
      task->SetRejectPairConversion(rejectPairConversion);
      task->SetVertexZMin(          zMin            );
      task->SetVertexZMax(          zMax            );
      task->SetVertexZWidth(        vZwidth         );
      task->SetEtaWidth(        etaBinWidth         );
      task->SetVertexXYMin(         -1.             );
      task->SetVertexXYMax(          1.             );
      task->SetCentralityMethod(    centralityMethod);
      task->SetCentrality(          minCentrality[iCentrality], maxCentrality[iCentrality]);
      task->SetPtMin1(              ptMin           );
      task->SetPtMax1(              ptMax           );
      task->SetPtBinWidth1(         ptWidthBin      );
      task->SetNPhiBins1(           nBinsPhi        );
      task->SetEtaMin1(             eta1Min         ); // SetYMin1 acturally
      task->SetEtaMax1(             eta1Max         ); // SetYMax1 acturally
      task->SetPtMin2(              ptMin           );
      task->SetPtMax2(              ptMax           );
      task->SetPtBinWidth2(         ptWidthBin      );
      task->SetNPhiBins2(           nBinsPhi        );
      task->SetEtaMin2(             eta2Min         ); // SetYMin2 acturally
      task->SetEtaMax2(             eta2Max         ); // SetYMax2 acturally
      task->SetDcaZMin(             dcaZMin         );
      task->SetDcaZMax(             dcaZMax         );
      task->SetDcaXYMin(            dcaXYMin        );
      task->SetDcaXYMax(            dcaXYMax        );
      task->SetDedxMin(             dedxMin         );
      task->SetDedxMax(             dedxMax         );
      task->SetNClusterMin(         nClusterMin     );
      task->SetTrackFilterBit(      trackFilterBit  );
      task->SetRequestedCharge_1(   requestedCharge1);
      task->SetRequestedCharge_2(   requestedCharge2);
      task->SetWeigth_1(            weight_1        );
      task->SetWeigth_2(            weight_2        );
      task->SetParticleSpecies(     particleID      );
      task->SetAnalysisType(        AnalysisDataType);
      task->SetSystemType(          System          );
      task->SetResonancesCut(       NoResonances    );
      task->SetElectronCut(         NoElectron      );
      task->SetSharedFractionPairCut( SharedFractionPairCut );
      task->SetNSigmaCut( nSigmaCut );
      task->SetNSigmaCut_veto( nSigmaCut_veto );
      task->SetPtCutUpperLimit( ptCUTupperMax );
      task->SetPtTOFlowerBoundary( ptTOFlowerMin );
      task->SetElectronNSigmaVetoCut( ElectronVetoCut );
      task->SetfRemoveTracksT0Fill( remove_Tracks_T0 );
      task->SetUse_AliHelperPID(  Use_AliHelperPID  );
      task->SetUse_CircularCutPID( Use_CircularCutPID );

      // assign initial values to AliHelperPID object
      AliHelperPID* helperpid = new AliHelperPID();
      helperpid -> SetNSigmaCut( nSigmaCut );
      helperpid -> SetPIDType( (PIDType_t)pidType );// kNSigmaTPC,kNSigmaTOF, kNSigmaTPCTOF
      helperpid -> SetfRequestTOFPID( requestTOFPID );
      helperpid -> SetfPtTOFPID( ptTOFlowerMin );
      helperpid -> SetisMC( isMC );
      task->SetHelperPID( helperpid );
      
      if(trigger) task -> SelectCollisionCandidates(AliVEvent::kINT7); //pPb, PbPb_2015
      else task -> SelectCollisionCandidates(AliVEvent::kMB); // PbPb & pp
        
      cout << "Creating task output container" << endl;
        
      taskOutputContainer = analysisManager->CreateContainer(listName,
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname));

      cout << "Add task to analysis manager and connect it to input and output containers" << endl;
        
      analysisManager->AddTask(task);
      analysisManager->ConnectInput( task,  0, analysisManager->GetCommonInputContainer());
      analysisManager->ConnectOutput(task,  1, taskOutputContainer );
      //analysisManager->ConnectOutput(task,  0, taskOutputContainer );
        
      iTask++;
    }            
  return task;
}
