//  Macro designed for use with the AliAnalysisTaskR2P2 task.
//  Author: (Baidyanath Sahoo, Basanta Kumar Nandi), IIT Bombay & Claude Pruneau, Wayne State University
//                      Year
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
// pp18_V0_kMB_kFALSE   18:     centralityMethod = 4 (V0),        trigger = kFALSE (AliVEvent::kMB).
/////////////////////////////////////////////////////////////////////////////////

AliAnalysisTaskR2P2 * AddTaskR2P2
(
 TString AnalysisDataType       = "MCAODreco", // "RealData"; "MCAOD" for MC AOD truth; "MCAODreco"
 TString prefixName             = "reco_",
 
 //************ for Stage2_Step1 starts ************************

 int    singlesOnly             =  1,   // 0: full correlations    1: singles only

 int    usePtEff                =  0,   // 0: no                   1: yes  
 char *inputPtEffFileName       = (char*)"alien:///alice/cern.ch/user/s/sbaidyan/efficiency_rmPi0nPhoton_14jan22.root", 
 
 int    useWeights              =  0,   // 0: no                   1: yes  
 char *inputHistogramFileName   = (char*)"alien:///alice/cern.ch/user/s/sbaidyan/data_lhc18i_runAll.root", 
 
 //************ for Stage2_Step1 ends ************************
 
 TString System                 = "pp18_V0_kMB_kFALSE",
 bool    pidparticle            =  0,   // 0: All Charged Particles;       1: PID particles
 bool    Use_PT_Cut             =  1,   // 0: Use_P_Cut ( only TOF lower boundary );       1: Use_PT_Cut
 int    useRapidity             =  0,   // 0: pseudo-rapadity      1: rapidity
 int    CentralityGroup         =  13,  // Diff Cent Groups dealing w/ memory limit & weight file 100M Alien limit
 int    PtGroup                 =  1,  // Diff Pt Groups dealing w/ memory limit & weight file 500M Alien limit
 double ptWidthBin              =  0.1, // pt bin width in histos
 double ptTOFlowerMin           =  0.8, // boundary between TPC & TOF region
 double zMax                    =  8.,  // set vertexZ cut   
 double vZwidth                 =  0.5, // zMin, zMax & vZwidth determine _nBins_vertexZ.
 int    trackFilterBit          =  96,/// 1,   // PbPb10(Global=1;TPConly=128;Hybrid=272); pPb13(Global=?;TPConly=?;Hybrid=768); pp10(Global=1;TPConly=?; Hybrid=?)
 double  chi2ndfMax             =  4.0,
 int    nClusterMin             =  70,//70
 double eta1Max                 =  0.8, // set y1max acturally if useRapidity==1
 double etaBinWidth             =  0.1, // set yBinWidth acturally if useRapidity==1
 double dcaZMax                 =  0.2,/// 3.2,
 double dcaXYMax                =  0.2,/// 2.4,
 int particleID                 =  3,   // Pion=0, Kaon=1, Proton=2, pionKaonProton = 3
 bool Use_CircularCutPID        =  0,///1,   // 0: Not Use_CircularCutPID     1: Use_CircularCutPID TPC+TOF
 double nSigmaCut               =  2.0,
 double nSigmaCut_veto          =  3.0,
 double ElectronVetoCut         =  1.0,
 int nBinsPhi                   =  72,  // 72 is default value
 Bool_t NoResonances            = kFALSE, // only for MCAOD
 Bool_t NoElectron              = kTRUE, // only for MCAOD
 bool   PurePIDinMC             = 1,///0,   // 0: MisID in MCAODreco;       1: No MisID in MCAODreco
 bool   PureNoWeakinMC          = 1,///0,   // 0: No MisID but Secondaries from weak decays in MCAODreco;       1: No MisID and No Secondaries from weak decays in MCAODreco
 bool   PureNoWeakMaterialinMC  = 1,///0,   // 0: No MisID and No Secondaries from weak decays but Secondaries from material in MCAODreco;       1: No MisID and No Secondaries from weak decays and material in MCAODreco
 
 bool   NoMisIDWeakMaterialInClosure = 1,///0,   // 0: allow MisID and No Secondaries from weak decays and material in MC Closure test       1: No MisID and Secondaries from weak decays and material in MC Closure test
 double SharedFractionPairCut   = 0.005 // check track splitting
 
 )
{
  // Set Default Configuration of this analysis
  // ==========================================
  int    debugLevel             = 0;
  int    rejectPileup           = 1;
  int    rejectPairConversion   = 1;
  int    sameFilter             = 1;
  int    centralityMethod       = 4; 
  Bool_t trigger                = kFALSE;
  Bool_t remove_Tracks_T0       = 0;///1;
  bool    useEventPlane         = 0;   // 0: No      1: Yes
  double  EventPlaneMin         = -3.1415927/6;
  double  EventPlaneMax         =  3.1415927/6;
  bool Use_AliHelperPID         =   0;   // 0: Not Use_AliHelperPID       1: Use_AliHelperPID
  //int pidType                   =  0;///2;  // kNSigmaTPC,kNSigmaTOF, kNSigmaTPCTOF // for AliHelperPID
  AliHelperPIDNameSpace::PIDType_t  pidType =  kNSigmaTPC;  // kNSigmaTPC,kNSigmaTOF, kNSigmaTPCTOF // for AliHelperPID
  Bool_t requestTOFPID          =  0;///baidya1;  // for AliHelperPID
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
  else if ( System == "pp18_V0_kMB_kFALSE" )    { centralityMethod = 4; trigger = kTRUE; }
  else if ( System == "pp_V0A_kMB_Utils" )    { centralityMethod = 7; trigger = kFALSE; }
  else if ( System == "pp_V0C_kMB_Utils" )    { centralityMethod = 8; trigger = kFALSE; }
  else if ( System == "pp_V0_kMB_Utils" )     { centralityMethod = 4; trigger = kFALSE; }
  else    return 0;
  
  int    chargeSet[4]; // 0: ++    1: +-    2: -+    3: --
  const char* taskname[4];
  int nChSet;
  
  if(singlesOnly){
    nChSet = 1;
    chargeSet[0] = 1;  taskname[0] = "ChPM";
  }
  else {
    nChSet = 3;
    chargeSet[0] = 1;  taskname[0] = "ChPM";
    chargeSet[1] = 0;  taskname[1] = "ChPP";                                 
    chargeSet[2] = 3;  taskname[2] = "ChMM";                                 
       
  }

  
  double ptMin[15]; // pt range lower limit cut ( also for pt histos )
  double ptCUTupperMax[15]; // pt range upper limit cut
  double ptMax[15]; // pt range upper limit for histos; NOT pt cut!!!
  int    nPt;
  
  if ( PtGroup == 1 )
    {
      nPt = 1;
      ptMin[0] = 0.2; ptCUTupperMax[0]  = 2.;   ptMax[0] = 3.0;
    }
  else  if ( PtGroup == 2 )
    {
      nPt = 2;
      ptMin[0] = 0.2; ptCUTupperMax[0]  = 2.;   ptMax[0] = 3.0;
      ptMin[1] = 2.;  ptCUTupperMax[1]  = 4.;   ptMax[1] = 5.0;
    }
  else  if ( PtGroup == 3 )
    {
      nPt = 3;
      ptMin[0] = 0.2; ptCUTupperMax[0]  = 2.;   ptMax[0] = 3.0;
      ptMin[1] = 2.;  ptCUTupperMax[1]  = 4.;   ptMax[1] = 5.0;
      ptMin[2] = 4.;  ptCUTupperMax[2]  = 10.;  ptMax[2] = 11.0;
    }
  else  if ( PtGroup == 5 )
    {
      nPt = 5;
      ptMin[0] = 0.2; ptCUTupperMax[0]  = 2.;   ptMax[0] = 3.0;
      ptMin[1] = 0.2;  ptCUTupperMax[1]  = 1.;   ptMax[1] = 2.0;
      ptMin[2] = 1.;  ptCUTupperMax[2]  = 2.;  ptMax[2] = 3.0;
      ptMin[3] = 2.;  ptCUTupperMax[3]  = 4.;  ptMax[3] = 5.0;
      ptMin[4] = 4.;  ptCUTupperMax[4]  = 10.;  ptMax[4] = 10.0;
    }
  else   if ( PtGroup == 9 )
    {
      nPt = 9;
      ptMin[0] = 0.2; ptCUTupperMax[0]  = 1.0;   ptMax[0] = 2.0;
      ptMin[1] = 1.0; ptCUTupperMax[1]  = 2.0;   ptMax[1] = 3.0;
      ptMin[2] = 2.0; ptCUTupperMax[2]  = 3.0;   ptMax[2] = 4.0;
      ptMin[3] = 3.0; ptCUTupperMax[3]  = 4.0;   ptMax[3] = 5.0;
      ptMin[4] = 4.0; ptCUTupperMax[4]  = 5.0;   ptMax[4] = 6.0;
      ptMin[5] = 4.0; ptCUTupperMax[5]  = 6.0;   ptMax[5] = 7.0;
      ptMin[6] = 4.0; ptCUTupperMax[6]  = 7.0;   ptMax[6] = 8.0;
      ptMin[7] = 4.0; ptCUTupperMax[7]  = 8.0;   ptMax[7] = 9.0;
      ptMin[8] = 4.0; ptCUTupperMax[8]  = 9.0;   ptMax[8] = 10.0;
      ptMin[9] = 4.0; ptCUTupperMax[9]  = 10.0;   ptMax[9] = 11.0;           
      ptMin[10] = 4.0; ptCUTupperMax[10]  = 100.0;   ptMax[10] = 101.0;           
      
    }

  else   if ( PtGroup == 8 )
    {
      nPt = 8;
      ptMin[0] = 0.2; ptCUTupperMax[0]  = 2.0;   ptMax[0] = 3.0;
      ptMin[1] = 2.; ptCUTupperMax[1]  = 4.0;   ptMax[1] = 5.0;
      ptMin[2] = 4.; ptCUTupperMax[2]  = 10.0;   ptMax[2] = 11.0;
      ptMin[3] = 4.; ptCUTupperMax[3]  = 100.0;   ptMax[3] = 101.0;

      ptMin[4] = 0.2; ptCUTupperMax[4]  = 1.0;   ptMax[4] = 2.0;
      ptMin[5] = 1.0; ptCUTupperMax[5]  = 2.0;   ptMax[5] = 3.0;
      ptMin[6] = 2.0; ptCUTupperMax[6]  = 3.0;   ptMax[6] = 4.0;
      ptMin[7] = 3.0; ptCUTupperMax[7]  = 4.0;   ptMax[7] = 5.0;

      
    }
  
  
  double minCentrality[15];
  double maxCentrality[15];
  int    nCentrality;
  
  
  if ( CentralityGroup == 1 )
    {
      nCentrality = 8;
      minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 30.;
      minCentrality[3] = 30.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 50.;
      minCentrality[5] = 50.;     maxCentrality[5]  = 60.;
      minCentrality[6] = 60.;     maxCentrality[6]  = 70.;
      minCentrality[7] = 70.;     maxCentrality[7]  = 80.; }
  else if ( CentralityGroup == 2 )
    {
      nCentrality = 4;
      minCentrality[0] = 40.;     maxCentrality[0]  = 50.;
      minCentrality[1] = 50.;     maxCentrality[1]  = 60.;
      minCentrality[2] = 60.;     maxCentrality[2]  = 70.;
      minCentrality[3] = 70.;     maxCentrality[3]  = 80.; }
  else if ( CentralityGroup == 3 )
    {
      nCentrality = 2;
      minCentrality[0] = 20.;     maxCentrality[0]  = 30.;
      minCentrality[1] = 30.;     maxCentrality[1]  = 40.; }
  else if ( CentralityGroup == 4 )
    {
      nCentrality = 2;
      minCentrality[0] = 60.;     maxCentrality[0]  = 70.;
      minCentrality[1] = 70.;     maxCentrality[1]  = 80.; }
  else if ( CentralityGroup == 5 )
    {
      nCentrality = 1;
      minCentrality[0] = 10.;     maxCentrality[0]  = 20.; }
  else if ( CentralityGroup == 6 )
    {
      nCentrality = 1;
      minCentrality[0] = 30.;     maxCentrality[0]  = 40.; }
  else if ( CentralityGroup == 7 )
    {
      nCentrality = 1;
      minCentrality[0] = 50.;     maxCentrality[0]  = 60.; }
  else if ( CentralityGroup == 8 )
    {
      nCentrality = 1;
      minCentrality[0] = 70.;     maxCentrality[0]  = 80.; }
  else if ( CentralityGroup == 9 )
    {
      nCentrality = 4;
      minCentrality[0] = 0;       maxCentrality[0]  = 20.;
      minCentrality[1] = 20.;     maxCentrality[1]  = 40.;
      minCentrality[2] = 40.;     maxCentrality[2]  = 60.;
      minCentrality[3] = 60.;     maxCentrality[3]  = 80.; }
  else if ( CentralityGroup == 10 )
    {
      nCentrality = 2;
      minCentrality[0] = 40.;     maxCentrality[0]  = 60.;
      minCentrality[1] = 60.;     maxCentrality[1]  = 80.; }
  else if ( CentralityGroup == 11 )
    {
      nCentrality = 1;
      minCentrality[0] = 20.;     maxCentrality[0]  = 40.; }
  else if ( CentralityGroup == 12 )
    {
      nCentrality = 1;
      minCentrality[0] = 60.;     maxCentrality[0]  = 80.; }
  else if ( CentralityGroup == 13 )
    {
      nCentrality = 1;
      minCentrality[0] = 0;       maxCentrality[0]  = 100.; }
  else if ( CentralityGroup == 14 )
    {
      nCentrality = 1;
      minCentrality[0] = 0;       maxCentrality[0]  = 80.; }
  else if ( CentralityGroup == 15 )
    {
      nCentrality = 9;
      minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 60.;
      minCentrality[7] = 60.;     maxCentrality[7]  = 70.;
      minCentrality[8] = 70.;     maxCentrality[8]  = 80.;  }
  else if ( CentralityGroup == 16 )
    {
      nCentrality = 1;
      minCentrality[0] = 20.;     maxCentrality[0]  = 100.; }
  else if ( CentralityGroup == 17 )
    {
      nCentrality = 1;
      minCentrality[0] = 0;       maxCentrality[0]  = 100.; }
  else if ( CentralityGroup == 18 )
    {
      nCentrality = 9;
      minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 60.;
      minCentrality[7] = 60.;     maxCentrality[7]  = 70.;
      minCentrality[8] = 70.;     maxCentrality[8]  = 100.;  }
  else if ( CentralityGroup == 19 )
    {
      nCentrality = 9;
      minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 60.;
      minCentrality[7] = 60.;     maxCentrality[7]  = 70.;
      minCentrality[8] = 70.;     maxCentrality[8]  = 90.;  }
  else if ( CentralityGroup == 20 )
    {
      nCentrality = 3;
      minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 30.;     maxCentrality[1]  = 40.;
      minCentrality[2] = 70.;     maxCentrality[2]  = 80.; }
  else if ( CentralityGroup == 21 )
    {
      nCentrality = 6;
      minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 30.;
      minCentrality[3] = 30.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 50.;
      minCentrality[5] = 50.;     maxCentrality[5]  = 80.; }
  else if ( CentralityGroup == 22 )
    {
      nCentrality = 3;
      minCentrality[0] = 0;       maxCentrality[0]  = 20.;
      minCentrality[1] = 20.;     maxCentrality[1]  = 40.;
      minCentrality[2] = 40.;     maxCentrality[2]  = 80.; }
  else if ( CentralityGroup == 23 )
    {
      nCentrality = 7;
      minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 80.; }
  else if ( CentralityGroup == 24 )
    {
      nCentrality = 8;
      minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 60.;
      minCentrality[7] = 60.;     maxCentrality[7]  = 80.; }
  else if ( CentralityGroup == 25 )
    {
      nCentrality = 7;
      minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 90.; }
  else if ( CentralityGroup == 26 )
    {
      nCentrality = 8;
      minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 60.;
      minCentrality[7] = 60.;     maxCentrality[7]  = 90.; }
  else if ( CentralityGroup == 27 )
    {
      nCentrality = 6;
      minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 30.;
      minCentrality[3] = 30.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 50.;
      minCentrality[5] = 50.;     maxCentrality[5]  = 90.; }
  else if ( CentralityGroup == 28 )
    {
      nCentrality = 5;
      minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 40.;
      minCentrality[3] = 40.;     maxCentrality[3]  = 60.;
      minCentrality[4] = 60.;     maxCentrality[4]  = 80.; }
  else if ( CentralityGroup == 29 )
    {
      nCentrality = 3;
      minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 30.;     maxCentrality[1]  = 40.;
      minCentrality[2] = 70.;     maxCentrality[2]  = 90.; }
  else if ( CentralityGroup == 30 )
    {
      nCentrality = 7;
      minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 30.;
      minCentrality[3] = 30.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 50.;
      minCentrality[5] = 50.;     maxCentrality[5]  = 70.;
      minCentrality[6] = 70.;     maxCentrality[6]  = 100.; }
  else if ( CentralityGroup == 31 )
    {
      nCentrality = 7;
      minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 60.;
      minCentrality[5] = 60.;     maxCentrality[5]  = 80.;
      minCentrality[6] = 80.;     maxCentrality[6]  = 100.; }
  else if ( CentralityGroup == 32 )
    {
      nCentrality = 8;
      minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 50.;
      minCentrality[6] = 50.;     maxCentrality[6]  = 70.;
      minCentrality[7] = 70.;     maxCentrality[7]  = 100.; }
  else if ( CentralityGroup == 33 )
    {
      nCentrality = 3;
      minCentrality[0] = 0;       maxCentrality[0]  = 20.;
      minCentrality[1] = 20.;     maxCentrality[1]  = 40.;
      minCentrality[2] = 40.;     maxCentrality[2]  = 90.; }
  else if ( CentralityGroup == 34 )
    {
      nCentrality = 6;
      minCentrality[0] = 0;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 30.;
      minCentrality[3] = 30.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 60.;
      minCentrality[5] = 60.;     maxCentrality[5]  = 90.; }
  else if ( CentralityGroup == 35 )
    {
      nCentrality = 7;
      minCentrality[0] = 0;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5.;      maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 30.;
      minCentrality[4] = 30.;     maxCentrality[4]  = 40.;
      minCentrality[5] = 40.;     maxCentrality[5]  = 60.;
      minCentrality[6] = 60.;     maxCentrality[6]  = 90.; }
  else if ( CentralityGroup == 36 )
    {
      nCentrality = 5;
      minCentrality[0] = 0;       maxCentrality[0]  = 20.;
      minCentrality[1] = 20.;      maxCentrality[1]  = 40.;
      minCentrality[2] = 40.;     maxCentrality[2]  = 100.;
      minCentrality[3] = 100.;     maxCentrality[3]  = 200.;
      minCentrality[4] = 200.;     maxCentrality[4]  = 400.;
    }
  //for pp13 TeV multiplicity classes
  else if ( CentralityGroup == 37 )
    {
      nCentrality = 6;
      minCentrality[0] = 2;       maxCentrality[0]  = 10.;
      minCentrality[1] = 10.;      maxCentrality[1]  = 20.;
      minCentrality[2] = 20.;     maxCentrality[2]  = 40.;
      minCentrality[3] = 40.;     maxCentrality[3]  = 100.;
      minCentrality[4] = 100.;     maxCentrality[4]  = 200.;
      minCentrality[5] = 200.;     maxCentrality[5]  = 400.;
    }
  else if ( CentralityGroup == 38 )
    {
      nCentrality = 1;
      minCentrality[0] = 2;       maxCentrality[0]  = 400.;
    }
  //pp 13TeV MB multiplicity classes
  else if ( CentralityGroup == 39 )
    {
      nCentrality = 10;
      minCentrality[0] = 0;       maxCentrality[0]  = 100.;
      minCentrality[1] = 0;       maxCentrality[1]  = 5.;
      minCentrality[2] = 5.;      maxCentrality[2]  = 10.;
      minCentrality[3] = 10.;     maxCentrality[3]  = 20.;
      minCentrality[4] = 20.;     maxCentrality[4]  = 30.;
      minCentrality[5] = 30.;     maxCentrality[5]  = 40.;
      minCentrality[6] = 40.;     maxCentrality[6]  = 50.;
      minCentrality[7] = 50.;     maxCentrality[7]  = 60.;
      minCentrality[8] = 60.;     maxCentrality[8]  = 70.;
      minCentrality[9] = 70.;     maxCentrality[9]  = 100.;  }

  //for pp13 TeV multiplicity classes
  else if ( CentralityGroup == 40 )
    {
      nCentrality = 9;
      minCentrality[0] = 2;       maxCentrality[0]  = 5.;
      minCentrality[1] = 5;       maxCentrality[1]  = 10.;
      minCentrality[2] = 10.;     maxCentrality[2]  = 20.;
      minCentrality[3] = 20.;     maxCentrality[3]  = 40.;
      minCentrality[4] = 40.;     maxCentrality[4]  = 100.;
      minCentrality[5] = 100.;    maxCentrality[5]  = 200.;
      minCentrality[6] = 200.;    maxCentrality[6]  = 400.;
      minCentrality[7] = 2.;      maxCentrality[7]  = 400.;
      minCentrality[8] = 5.;      maxCentrality[8]  = 400.;
    }
  else if ( CentralityGroup == 41 )
    {
      nCentrality = 11;
      minCentrality[0] = 0;       maxCentrality[0]  = 100.;
      minCentrality[1] = 0;       maxCentrality[1]  = 1.;
      minCentrality[2] = 1.;      maxCentrality[2]  = 5.;
      minCentrality[3] = 5.;      maxCentrality[3]  = 10.;
      minCentrality[4] = 10.;     maxCentrality[4]  = 15.;
      minCentrality[5] = 15.;     maxCentrality[5]  = 20.;
      minCentrality[6] = 20.;     maxCentrality[6]  = 30.;
      minCentrality[7] = 30.;     maxCentrality[7]  = 40.;
      minCentrality[8] = 40.;     maxCentrality[8]  = 50.;
      minCentrality[9] = 50.;     maxCentrality[9]  = 70.;
      minCentrality[10] = 70.;     maxCentrality[10]  = 100.;
    }
  else if ( CentralityGroup == 42 )
    {
      nCentrality = 6;
      minCentrality[0] = 0;       maxCentrality[0]  = 100.;
      minCentrality[1] = 0;       maxCentrality[1]  = 1.;
      minCentrality[2] = 1.;      maxCentrality[2]  = 5.;
      minCentrality[3] = 5.;      maxCentrality[3]  = 10.;
      minCentrality[4] = 10.;     maxCentrality[4]  = 15.;
      minCentrality[5] = 15.;     maxCentrality[5]  = 20.;

    }

  else if ( CentralityGroup == 43 )
    {
      nCentrality = 5;

      minCentrality[0] = 20.;     maxCentrality[0]  = 30.;
      minCentrality[1] = 30.;     maxCentrality[1]  = 40.;
      minCentrality[2] = 40.;     maxCentrality[2]  = 50.;
      minCentrality[3] = 50.;     maxCentrality[3]  = 70.;
      minCentrality[4] = 70.;     maxCentrality[4]  = 100.;
    }
  else if ( CentralityGroup == 44 )
    {
      nCentrality = 3;

      minCentrality[0] = 0.;     maxCentrality[0]  = 100.;
      minCentrality[1] = 0.;     maxCentrality[1]  = 90.;
      minCentrality[2] = 0.;     maxCentrality[2]  = 80.;

    }
  else if ( CentralityGroup == 45 )
    {
      nCentrality = 8;

      minCentrality[0] = 0.;     maxCentrality[0]  = 100.;
      minCentrality[1] = 0.;     maxCentrality[1]  = 5.;
      minCentrality[2] = 5.;     maxCentrality[2]  = 10.;
      minCentrality[3] = 10.;     maxCentrality[3]  = 20.;
      minCentrality[4] = 80.;     maxCentrality[4]  = 85.;
      minCentrality[5] = 85.;     maxCentrality[5]  = 90.;
      minCentrality[6] = 90.;     maxCentrality[6]  = 95.;
      minCentrality[7] = 95.;     maxCentrality[7]  = 100.;
    }
  else if ( CentralityGroup == 46 )
    {
      nCentrality = 10;

      minCentrality[0] = 0.;     maxCentrality[0]  = 100.;
      minCentrality[1] = 0.;     maxCentrality[1]  = 90.;
      minCentrality[2] = 0.;     maxCentrality[2]  = 80.;
      minCentrality[3] = 0.;     maxCentrality[3]  = 5.;
      minCentrality[4] = 5.;     maxCentrality[4]  = 10.;
      minCentrality[5] = 10.;     maxCentrality[5]  = 20.;
      minCentrality[6] = 80.;     maxCentrality[6]  = 85.;
      minCentrality[7] = 85.;     maxCentrality[7]  = 90.;
      minCentrality[8] = 90.;     maxCentrality[8]  = 95.;
      minCentrality[9] = 95.;     maxCentrality[9]  = 100.;
    }

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
      printf("AddTaskR2P2 No analysis manager to connect to.");
      return NULL;
    }

   
 
  TString part1Name;
  TString part2Name;
  TString partName;
  TString eventName;

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
  AliAnalysisTaskR2P2* task;
  
  for (int iChSet=0; iChSet < nChSet; ++iChSet)
    {
      for (int iCentrality=0; iCentrality < nCentrality; ++iCentrality)
	{
	  for (int iPt = 0; iPt < nPt; ++iPt)
	    {
	      switch ( chargeSet[iChSet] )
		{
		case 0: part1Name = "P1_"; part2Name = "P2"; requestedCharge1 =  1; requestedCharge2 =  1; sameFilter = 1;   break;
		case 1: part1Name = "P1_"; part2Name = "M2"; requestedCharge1 =  1; requestedCharge2 = -1; sameFilter = 0;   break;
		case 2: part1Name = "M1_"; part2Name = "P2"; requestedCharge1 = -1; requestedCharge2 =  1; sameFilter = 0;   break;
		case 3: part1Name = "M1_"; part2Name = "M2"; requestedCharge1 = -1; requestedCharge2 = -1; sameFilter = 1;   break;
		}
	      
	      

	      partName = "pt";
	      partName += int(10*ptMin[iPt]);
	      partName += "to";
	      partName += int(10*ptCUTupperMax[iPt]);

	      part1Name += "pt";
	      part1Name += int(10*ptMin[iPt]);
	      part1Name += "to";
	      part1Name += int(10*ptCUTupperMax[iPt]);

	      eventName =  "_";
	      eventName += int(minCentrality[iCentrality] );
	      eventName += "Vo";
	      eventName += int(maxCentrality[iCentrality] );
	      

	      baseName     =   prefixName;	      
	      baseName     +=  part1Name;
	      baseName     +=  "_";
	      baseName     +=  part2Name;
	      baseName     +=  eventName;
	      baseName     +=  "_";
	      baseName     +=  particleID;
	      listName     =   baseName;
	      taskName     =   baseName;        
	      
	      outputHistogramFileName = baseName;
	      outputHistogramFileName += ".root";      





	      
	      //------pt-efficiency starts ---------------
	      TFile  * inputFilePtEff  = 0;
	      TH1F   * hPtEff_1   = 0;
	      TH1F   * hPtEff_2   = 0;
	      if ( usePtEff )
		{
		  TGrid::Connect("alien:");
		  inputFilePtEff = TFile::Open(inputPtEffFileName,"OLD");
		  if (inputFilePtEff)
		    {
		      cout << "\n\n\n\n ====================================================" << endl;
		      cout << " Requested file:" << inputPtEffFileName << " was opened." << endl;
		      cout << "\n\n\n\n ====================================================" << endl;
		    }
		  else if (!inputFilePtEff)
		    {
		      cout << "\n\n\n\n ====================================================" << endl;
		      cout << " Requested file:" << inputPtEffFileName << " was not opened. ABORT." << endl;
		      cout << "\n\n\n\n ====================================================" << endl;
		      return 0;
		    }
		  TString nameHistoBasePtEff = "hEff_";
		  TString nameHistoPtEff;
		  nameHistoBasePtEff += partName;
		  nameHistoBasePtEff += eventName;
		  cout <<"\n\n\n nameHistoBasePtEff: "<<nameHistoBasePtEff<<endl;
		  if (requestedCharge1 == 1)
		    {
		      nameHistoPtEff = nameHistoBasePtEff + "_p";
		      cout << "\n\n\n Input Histogram named: " << nameHistoPtEff << endl;
		      hPtEff_1 = (TH1F *) inputFilePtEff->Get(nameHistoPtEff);
		    }
		  else
		    {
		      nameHistoPtEff = nameHistoBasePtEff + "_m";
		      cout << "\n\n\n Input Histogram named: " << nameHistoPtEff << endl;
		      hPtEff_1 = (TH1F *) inputFilePtEff->Get(nameHistoPtEff);
		    }
		  if (!hPtEff_1)
		    {
		      cout << "\n\n\n Requested histogram 'ptEff_p/m' was not found. ABORT." << endl;
		      return 0;
		    }
		  
		  if (!sameFilter)
		    {
		      hPtEff_2 = 0;
		      if (requestedCharge2 == 1)
			{
			  nameHistoPtEff = nameHistoBasePtEff + "_p";
			  cout << "\n\n\n Input Histogram named: " << nameHistoPtEff << endl;
			  hPtEff_2 = (TH1F *) inputFilePtEff->Get(nameHistoPtEff);
			}
		      else
			{
			  nameHistoPtEff = nameHistoBasePtEff + "_m";
			  cout << "\n\n\n Input Histogram named: " << nameHistoPtEff << endl;
			  hPtEff_2 = (TH1F *) inputFilePtEff->Get(nameHistoPtEff);
			}
		      if (!hPtEff_2)
			{
			  cout << "\n \n \n Requested histogram 'ptEff_p/m' was not found. ABORT." << endl;
			  return 0;
			}
		    }
		}


	      //------pt-efficiency ends ---------------


	      

	      
	      TFile  * inputFile  = 0;
	      TList  * histoList  = 0;
	      TH3F   * weight_1   = 0;
	      TH3F   * weight_2   = 0;
	      if ( useWeights )
		{
		  TGrid::Connect("alien:");
		  inputFile = TFile::Open(inputHistogramFileName,"OLD");
		  if (inputFile)
		    {
		      cout << "\n\n\n\n ====================================================" << endl;
		      cout << " Requested file:" << inputHistogramFileName << " was opened." << endl;
		      cout << "\n\n\n\n ====================================================" << endl;
		    }
		  else if (!inputFile)
		    {
		      cout << "\n\n\n\n ====================================================" << endl;
		      cout << " Requested file:" << inputHistogramFileName << " was not opened. ABORT." << endl;
		      cout << "\n\n\n\n ====================================================" << endl;
		      return 0;
		    }
		  TString nameHistoBase = "correction_";
		  TString nameHisto;
		  nameHistoBase += partName;
		  nameHistoBase += eventName;
		  cout <<"\n\n\n nameHistoBase: "<<nameHistoBase<<endl;
		  if (requestedCharge1 == 1)
		    {
		      nameHisto = nameHistoBase + "_p";
		      //    cout << "\n\n\n Input Histogram named: " << nameHisto << endl;
		      weight_1 = (TH3F *) inputFile->Get(nameHisto);
		    }
		  else
		    {
		      nameHisto = nameHistoBase + "_m";
		      // cout << "\n\n\n Input Histogram named: " << nameHisto << endl;
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
	      //       cout << "baidyaTest before task: "  << endl;
	      task = new  AliAnalysisTaskR2P2(taskName);
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
	      task->SetUsePtEff(          usePtEff      );
	      task->SetUseWeights(          useWeights      );
	      task->SetUseRapidity(         useRapidity     );
	      task->SetEventPlane(         useEventPlane     );
	      task->SetEPmin(              EventPlaneMin     );
	      task->SetEPmax(              EventPlaneMax     );
	      task->SetRejectPileup(        rejectPileup    );
	      task->SetRejectPairConversion(rejectPairConversion);
	      task->SetVertexZMin(          zMin            );
	      task->SetChi2PerNDF(          chi2ndfMax      );        
	      task->SetVertexZMax(          zMax            );
	      task->SetVertexZWidth(        vZwidth         );
	      task->SetEtaWidth(        etaBinWidth         );
	      task->SetVertexXYMin(         -1.             );
	      task->SetVertexXYMax(          1.             );
	      task->SetCentralityMethod(    centralityMethod);
	      task->SetCentrality(          minCentrality[iCentrality], maxCentrality[iCentrality]);
	      task->SetPtMin1(              ptMin[iPt]      );
	      task->SetPtMax1(              ptMax[iPt]      );
	      task->SetPtBinWidth1(         ptWidthBin      );
	      task->SetNPhiBins1(           nBinsPhi        );
	      task->SetEtaMin1(             eta1Min         ); // SetYMin1 acturally
	      task->SetEtaMax1(             eta1Max         ); // SetYMax1 acturally
	      task->SetPtMin2(              ptMin[iPt]      );
	      task->SetPtMax2(              ptMax[iPt]      );
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
	      task->SetPtEff_1(            hPtEff_1        );
	      task->SetPtEff_2(            hPtEff_2        );
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
	      task->SetPtCutUpperLimit( ptCUTupperMax[iPt] );
	      task->SetPtTOFlowerBoundary( ptTOFlowerMin );
	      task->SetElectronNSigmaVetoCut( ElectronVetoCut );
	      task->SetfRemoveTracksT0Fill( remove_Tracks_T0 );
	      task->SetUse_AliHelperPID(  Use_AliHelperPID  );
	      task->SetUse_CircularCutPID( Use_CircularCutPID );
	      
	      // assign initial values to AliHelperPID object
	      AliHelperPID* helperpid = new AliHelperPID();
	      helperpid -> SetNSigmaCut( nSigmaCut );
	      //helperpid -> SetPIDType( (PIDType_t)pidType );// kNSigmaTPC,kNSigmaTOF, kNSigmaTPCTOF
	      helperpid -> SetPIDType( pidType );
	      helperpid -> SetfRequestTOFPID( requestTOFPID );
	      helperpid -> SetfPtTOFPID( ptTOFlowerMin );
	      helperpid -> SetisMC( isMC );
	      task->SetHelperPID( helperpid );
	      
	      if(trigger) task -> SelectCollisionCandidates(AliVEvent::kINT7); //pPb, PbPb_2015
	      else task -> SelectCollisionCandidates(AliVEvent::kMB); // PbPb & pp
	      
	      cout << "Creating task output container" << endl;
	      
	      taskOutputContainer =
		analysisManager->CreateContainer(listName, TList::Class(),
						 AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname[iChSet]));
	      
	      cout << "Add task to analysis manager and connect it to input and output containers" << endl;
	      
	      analysisManager->AddTask(task);
	      analysisManager->ConnectInput( task,  0, analysisManager->GetCommonInputContainer());
	      analysisManager->ConnectOutput(task,  1, taskOutputContainer );
	      //analysisManager->ConnectOutput(task,  0, taskOutputContainer );
	      
	      iTask++;
	    }//iPt loop
	} //iCent loop
    } //iChSet loop
  
  return task;
  
}
