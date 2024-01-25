//  Macro designed for use with the AliAnalysisTaskR2P2spheroDiffpT task.
//  Author: Subhadeep Roy, IIT Bombay & Claude Pruneau, Wayne State University
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
//   pp18_V0_kMB_kFALSE   18:     centralityMethod = 4 (V0),        trigger = kFALSE (AliVEvent::kMB).
/////////////////////////////////////////////////////////////////////////////////


AliAnalysisTaskR2P2spheroDiffpT * AddTaskR2P2spheroDiffpT
(

  TString AnalysisDataType       = "RealData", // "RealData"; "MCAOD" for MC AOD truth; "MCAODreco"
  TString prefixName             = "data_",  
  int    singlesOnly             =  0,   // 0: full correlations    1: singles only // (0 for stage 2 & 1 for stage 1)
  int    usePtEff                =  1,   // 0: no                   1: yes
  TString inputPtEffFileName     = "alien:///alice/cern.ch/user/s/subhadee/inv_efficiency_correction_factor.root",
    
  int    chargeSet               =  1,   // 0: ++    1: +-    2: -+    3: --
  const char* taskname           = "ChPM",// ChPM, ChPP, ChMM
  //************ for Stage2_Step1 ends ************************
  
  TString System                 = "pp18_V0_kMB_kFALSE", //for pp system
  int    CentralityGroup         =  13,  // Diff Cent Groups dealing w/ memory limit & weight file 100M Alien limit
  int    SpherocityGroup         =  4,   //spherocity groups
  int    PtGroup                 =  2,  // Diff Pt Groups dealing w/ memory limit & weight file 500M Alien limit

  double zMax                    =  8.,  // set vertexZ cut   
  double vZwidth                 =  0.5, // zMin, zMax & vZwidth determine _nBins_vertexZ.
  int    trackFilterBit          =  96,   // PbPb10(Global=1;TPConly=128;Hybrid=272); pPb13(Global=?;TPConly=?;Hybrid=768); pp10(Global=1;TPConly=?; Hybrid=?)
  int    nClusterMin             =  70,
  double eta1Max                 =  0.8, // set y1max acturally if useRapidity==1
  double etaBinWidth             =  0.1, // set yBinWidth acturally if useRapidity==1
  double dcaZMax                 =  0.2,//3.2,
  double dcaXYMax                =  0.2,//2.4,
  double ElectronVetoCut         =  1.0,
  int nBinsPhi                   =  72,  // 36 is default value
  //double ptMin                   =  0.2, // pt range lower limit cut ( also for pt histos )
  //double ptCUTupperMax           =  2.0, // pt range upper limit cut
  //double ptMax                   =  3.0, // pt range upper limit for histos; NOT pt cut!!!
  Bool_t NoResonances            = kFALSE, // only for MCAOD
  Bool_t NoElectron              = kTRUE, // only for MCAOD
  bool   PurePIDinMC             = 1,   // 0: MisID in MCAODreco;       1: No MisID in MCAODreco
  bool   PureNoWeakinMC          = 1,   // 0: No MisID but Secondaries from weak decays in MCAODreco;       1: No MisID and No Secondaries from weak decays in MCAODreco
  bool   PureNoWeakMaterialinMC  = 1,   // 0: No MisID and No Secondaries from weak decays but Secondaries from material in MCAODreco;       1: No MisID and No Secondaries from weak decays and material in MCAODreco
  bool   NoMisIDWeakMaterialInClosure = 1,   // 0: allow MisID and No Secondaries from weak decays and material in MC Closure test       1: No MisID and Secondaries from weak decays and material in MC Closure test
  double SharedFractionPairCut   = 0.005, // check track splitting
  int    pileUpEventPP13             =  1,   // 0: with Event Pileup    1: without Event Pileup
  int    pileUpTrackPP13             =  1   // 0: with track Pileup    1: without track Pileup
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
  bool Use_AliHelperPID         =  0;   // 0: Not Use_AliHelperPID       1: Use_AliHelperPID
  int pidType                   = 0;///2;  // kNSigmaTPC,kNSigmaTOF, kNSigmaTPCTOF // for AliHelperPID
  
  Bool_t requestTOFPID          = 0;///baidya1;  // for AliHelperPID
  
  Bool_t isMC                   =  0;  // for AliHelperPID
  double eta2Max                =  eta1Max; // set y2max acturally if useRapidity==1
  double eta1Min                = -eta1Max; // set y1min acturally if useRapidity==1
  double eta2Min                = -eta1Max; // set y2min acturally if useRapidity==1
  double dcaZMin                = -dcaZMax;
  double dcaXYMin               = -dcaXYMax;
  double zMin                   = -zMax;  // |zMin| should = zMax due to the design of the code

  bool    pidparticle            =  0;   // 0: All Charged Particles;       1: PID particles
  bool    Use_PT_Cut             =  1;   // 0: Use_P_Cut ( only TOF lower boundary );       1: Use_PT_Cut
  int    useRapidity             =  0;   // 0: pseudo-rapadity      1: rapidity
  
  double ptTOFlowerMin           =  0.5; // boundary between TPC & TOF region
  double ptWidthBin              =  0.1; // pt bin width in histos

  int particleID                 =  3;   // Pion=0, Kaon=1, Proton=2, pionKaonProton = 3
  bool Use_CircularCutPID        =  0;   // 0: Not Use_CircularCutPID     1: Use_CircularCutPID TPC+TOF
  double nSigmaCut               =  2.0;
  double nSigmaCut_veto          =  3.0;

  
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
  else if ( System == "pp18_V0_kMB_kFALSE" )    { centralityMethod = 4; trigger = kTRUE; }
  else    return 0;


  //different pT groups

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
    ptMin[0] = 0.2; ptCUTupperMax[0]  = 1.;   ptMax[0] = 2.0;
    ptMin[1] = 1.;  ptCUTupperMax[1]  = 2.;   ptMax[1] = 3.0;
  }

  //Different spherocity class values


  double minSherocity[15];
  double maxSherocity[15];
  int    nSherocity;

  //For 0-100% multiplicity class with Nch>=5

  if ( SpherocityGroup ==1) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=1;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.985; //0- 100% sphero class
  }

  else if ( SpherocityGroup ==2) //corresponds to the spherocity for  Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.435; //0- 20%
    minSherocity[1] = 0.755;       maxSherocity[1]  =0.985;  //80-100%
  }

  else if ( SpherocityGroup ==3) //corresponds to the spherocity for  Gen/Truth level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.455; //0- 20%
    minSherocity[1] = 0.765;       maxSherocity[1]  =0.985;  //80-100%
  }

  //For 0-10% multiplicity class with Nch>=5

  else if ( SpherocityGroup ==4) //corresponds to the spherocity for  Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.595; //0- 20%
    minSherocity[1] = 0.825;       maxSherocity[1]  =0.985;  //80-100%
  }

  else if ( SpherocityGroup ==5) //corresponds to the spherocity for  Gen/Truth level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.635; //0- 20%
    minSherocity[1] = 0.845;       maxSherocity[1]  =0.985;  //80-100%
  }


  //For 0-10% multiplicity class with Nch>=10

  else if ( SpherocityGroup ==6) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.595; //0- 20%
    minSherocity[1] = 0.825;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  else if ( SpherocityGroup ==7) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.635; //0- 20%
    minSherocity[1] = 0.845;       maxSherocity[1]  =0.985;  //80-100%
    
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
    minCentrality[7] = 70.;     maxCentrality[7]  = 80.; 
  }

  else if ( CentralityGroup == 2 )
  {
    nCentrality = 4;
    minCentrality[0] = 40.;     maxCentrality[0]  = 50.;
    minCentrality[1] = 50.;     maxCentrality[1]  = 60.;
    minCentrality[2] = 60.;     maxCentrality[2]  = 70.;
    minCentrality[3] = 70.;     maxCentrality[3]  = 80.; 
  }

  else if ( CentralityGroup == 3 )
  {
    nCentrality = 2;
    minCentrality[0] = 20.;     maxCentrality[0]  = 30.;
    minCentrality[1] = 30.;     maxCentrality[1]  = 40.; 
  }

  else if ( CentralityGroup == 4 )
  {
    nCentrality = 2;
    minCentrality[0] = 60.;     maxCentrality[0]  = 70.;
    minCentrality[1] = 70.;     maxCentrality[1]  = 80.; 
  }

  else if ( CentralityGroup == 5 )
  {
    nCentrality = 1;
    minCentrality[0] = 10.;     maxCentrality[0]  = 20.;  
  }
  else if ( CentralityGroup == 6 )
  {
    nCentrality = 1;
    minCentrality[0] = 30.;     maxCentrality[0]  = 40.; 
  }
  else if ( CentralityGroup == 7 )
  {
    nCentrality = 1;
    minCentrality[0] = 50.;     maxCentrality[0]  = 60.; 
  }
  else if ( CentralityGroup == 8 )
  {
    nCentrality = 1;
    minCentrality[0] = 70.;     maxCentrality[0]  = 80.; 
  }
  else if ( CentralityGroup == 9 )
  {
    nCentrality = 4;
    minCentrality[0] = 0;       maxCentrality[0]  = 20.;
    minCentrality[1] = 20.;     maxCentrality[1]  = 40.;
    minCentrality[2] = 40.;     maxCentrality[2]  = 60.;
    minCentrality[3] = 60.;     maxCentrality[3]  = 80.; 
  }
  else if ( CentralityGroup == 10 )
  {
    nCentrality = 2;
    minCentrality[0] = 40.;     maxCentrality[0]  = 60.;
    minCentrality[1] = 60.;     maxCentrality[1]  = 80.; 
  }
  else if ( CentralityGroup == 11 )
  {
    nCentrality = 1;
    minCentrality[0] = 20.;     maxCentrality[0]  = 40.; 
  }
  else if ( CentralityGroup == 12 )
  {
    nCentrality = 1;
    minCentrality[0] = 60.;     maxCentrality[0]  = 80.; 
  }
  else if ( CentralityGroup == 13 )
  {
    nCentrality = 1;
    minCentrality[0] = 0;       maxCentrality[0]  = 10.; 
  }
  else if ( CentralityGroup == 14 )
  {
    nCentrality = 1;
    minCentrality[0] = 0;       maxCentrality[0]  = 80.; 
  }
  
  else if ( CentralityGroup == 101 )
  {
    nCentrality = 1;
    minCentrality[0] = 0;       maxCentrality[0]  = 1.;
  }
  else if ( CentralityGroup == 102 )
  {
    nCentrality = 1;
    minCentrality[0] = 1.;      maxCentrality[0]  = 5.;
  }
  else if ( CentralityGroup == 103 )
  {
    nCentrality = 1;
    minCentrality[0] = 5.;      maxCentrality[0]  = 10.;
  }
  else if ( CentralityGroup == 104 )
  {
    nCentrality = 1;
    minCentrality[0] = 10.;	 maxCentrality[0]  = 15.;
  }
  else if ( CentralityGroup == 105 )
  {
    nCentrality = 1;
    minCentrality[0] = 15.;	 maxCentrality[0]  = 20.;
  }
  else if ( CentralityGroup == 106 )
  {
    nCentrality = 1;
    minCentrality[0] = 20.;	 maxCentrality[0]  = 30.;
  }
  else if ( CentralityGroup == 107 )
  {
    nCentrality = 1;
    minCentrality[0] = 30.;	 maxCentrality[0]  = 40.;
  }
  else if ( CentralityGroup == 108 )
  {
    nCentrality = 1;
    minCentrality[0] = 40.;	 maxCentrality[0]  = 50.;
  }
  else if ( CentralityGroup == 109 )
  {
    nCentrality = 1;
    minCentrality[0] = 50.;	 maxCentrality[0]  = 70.;
  }
  else if ( CentralityGroup == 110 )
  {
    nCentrality = 1;
    minCentrality[0] = 70.;	 maxCentrality[0]  = 100.;
  }
  else if ( CentralityGroup == 111 )
  {
    nCentrality = 2;
    minCentrality[0] = 0;       maxCentrality[0]  = 1.;
    minCentrality[1] = 1.;      maxCentrality[1]  = 5.;
  }
  else if ( CentralityGroup == 112 )
  {
    nCentrality = 2;
    minCentrality[0] = 5.;      maxCentrality[0]  = 10.;
    minCentrality[1] = 10.;	  maxCentrality[1]  = 15.;
  }
  else if ( CentralityGroup == 113 )
  {
    nCentrality = 2;
    minCentrality[0] = 15.;	 maxCentrality[0]  = 20.;
    minCentrality[1] = 20.;	 maxCentrality[1]  = 30.;
  }
  else if ( CentralityGroup == 114 )
  {
    nCentrality = 2;
    minCentrality[0] = 30.;	 maxCentrality[0]  = 40.;
    minCentrality[1] = 40.;	 maxCentrality[1]  = 50.;
  }
  else if ( CentralityGroup == 115 )
  {
    nCentrality = 2;
    minCentrality[0] = 50.;	 maxCentrality[0]  = 70.;
    minCentrality[1] = 70.;	 maxCentrality[1]  = 100.;
  }
  else if ( CentralityGroup == 121 ) //victor
  {
    nCentrality = 1;
    minCentrality[0] = 0;       maxCentrality[0]  = 5.;
  }
  else if ( CentralityGroup == 122 ) //victor
  {
    nCentrality = 2;
    minCentrality[0] = 5.;	  maxCentrality[0]  = 10.;
    minCentrality[1] = 10.;	  maxCentrality[1]  = 20.;
  }
  else if ( CentralityGroup == 123 ) //victor
  {
    nCentrality = 2;
    minCentrality[0] = 20.;	  maxCentrality[0]  = 30.;
    minCentrality[1] = 30.;	  maxCentrality[1]  = 40.;
  }
  else if ( CentralityGroup == 124 ) //victor
  {
    nCentrality = 2;
    minCentrality[0] = 40.;	  maxCentrality[0]  = 50.;
    minCentrality[1] = 50.;	  maxCentrality[1]  = 60.;
  }
  else if ( CentralityGroup == 125 ) //victor
  {
    nCentrality = 2;
    minCentrality[0] = 60.;	  maxCentrality[0]  = 70.;
    minCentrality[1] = 70.;	  maxCentrality[1]  = 80.;
  }
  else if ( CentralityGroup == 131 ) //victor
  {
    nCentrality = 1;
    minCentrality[0] = 0;       maxCentrality[0]  = 5.;
  }
  else if ( CentralityGroup == 132 ) //victor
  {
    nCentrality = 1;
    minCentrality[0] = 5.;	  maxCentrality[0]  = 10.;
  }
  else if ( CentralityGroup == 133 ) //victor
  {
    nCentrality = 1;
    minCentrality[0] = 10.;	  maxCentrality[0]  = 20.;
  }
  else if ( CentralityGroup == 134 ) //victor
  {
    nCentrality = 1;
    minCentrality[0] = 20.;	  maxCentrality[0]  = 30.;
  }
  else if ( CentralityGroup == 135 ) //victor
  {
    nCentrality = 1;
    minCentrality[0] = 30.;	  maxCentrality[0]  = 40.;
  }
  else if ( CentralityGroup == 136 ) //victor
  {
    nCentrality = 1;
    minCentrality[0] = 40.;	  maxCentrality[0]  = 50.;
  }
  else if ( CentralityGroup == 137 ) //victor
  {
    nCentrality = 1;
    minCentrality[0] = 50.;	  maxCentrality[0]  = 60.;
  }
  else if ( CentralityGroup == 138 ) //victor
  {
    nCentrality = 1;
    minCentrality[0] = 60.;	  maxCentrality[0]  = 70.;
  }
  else if ( CentralityGroup == 139 ) //victor
  {
    nCentrality = 1;
    minCentrality[0] = 70.;	  maxCentrality[0]  = 80.;
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
    ::Error("AddTaskR2P2spheroDiffpT", "No analysis manager to connect to.");
    return NULL;
  }
    
  TString part1Name;
  TString part2Name;
  TString partName;
  TString eventName;
  TString ptName;
  TString eventName1;

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
  AliAnalysisTaskR2P2spheroDiffpT* task;


  for (int iCentrality=0; iCentrality < nCentrality; ++iCentrality)
  {
    for (int isphero = 0; isphero < nSherocity; ++isphero)
    {
      eventName  = "So";
      eventName += int(1000*minSherocity[isphero] );
      eventName += "to";
      eventName += int(1000*maxSherocity[isphero] );
      eventName += "_";

      for (int iPt = 0; iPt < nPt; ++iPt)
      {
        switch ( chargeSet )
        {
          case 0: part1Name = "P"; part2Name = "P"; requestedCharge1 =  1; requestedCharge2 =  1; sameFilter = 1;   break;
          case 1: part1Name = "P"; part2Name = "M"; requestedCharge1 =  1; requestedCharge2 = -1; sameFilter = 0;   break;
          case 2: part1Name = "M"; part2Name = "P"; requestedCharge1 = -1; requestedCharge2 =  1; sameFilter = 0;   break;
          case 3: part1Name = "M"; part2Name = "M"; requestedCharge1 = -1; requestedCharge2 = -1; sameFilter = 1;   break;
        }

        partName = "pt2to20";

        part1Name += part2Name;
        part1Name += "_";

      
        ptName  = "pt";
        ptName += int(10*ptMin[iPt]);
        ptName += "to";
        ptName += int(10*ptCUTupperMax[iPt]);
        ptName += "_";

        baseName     = prefixName + part1Name + eventName + ptName + particleID;
        
        listName     =   baseName;
        taskName     =   baseName;
        outputHistogramFileName = baseName;
        outputHistogramFileName += ".root";

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

          TString nameHistoPtEff;

          TString nameHistoBasePtEff = "hEff_";
          nameHistoBasePtEff += partName;
          nameHistoBasePtEff += "_0Vo100";

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

          } //sameFilter
        }//Use pT efficiency

        task = new  AliAnalysisTaskR2P2spheroDiffpT(taskName);
        //configure my task
        task->SetDebugLevel(          debugLevel      );
        task->SetSameFilter(          sameFilter      );
        task->SetSinglesOnly(         singlesOnly     );
        task->SetPileUpEventPP13(     pileUpEventPP13 );
        task->SetPileUpTrackPP13(     pileUpTrackPP13 );
        task->SetPIDparticle(         pidparticle     );
        task->SetUse_pT_cut(          Use_PT_Cut      );
        task->SetIfContaminationInMC(   PurePIDinMC   );
        task->SetIfContaminationWeakInMC( PureNoWeakinMC );
        task->SetIfContaminationWeakMaterialInMC( PureNoWeakMaterialinMC );
        task->SetIfMisIDWeakMaterialInMCClosure( NoMisIDWeakMaterialInClosure );
        task->SetUsePtEff(          usePtEff      );
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
        task->SetTSpherocity(          minSherocity[isphero], maxSherocity[isphero]);
        task->SetPtMin1(              ptMin[iPt]           );
        task->SetPtMax1(              ptMax[iPt]         );
        task->SetPtBinWidth1(         ptWidthBin      );
        task->SetNPhiBins1(           nBinsPhi        );
        task->SetEtaMin1(             eta1Min         ); // SetYMin1 acturally
        task->SetEtaMax1(             eta1Max         ); // SetYMax1 acturally
        task->SetPtMin2(              ptMin[iPt]           );
        task->SetPtMax2(              ptMax[iPt]           );
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
        helperpid -> SetPIDType( (AliHelperPIDNameSpace::PIDType_t)pidType );// kNSigmaTPC,kNSigmaTOF, kNSigmaTPCTOF
        helperpid -> SetfRequestTOFPID( requestTOFPID );
        helperpid -> SetfPtTOFPID( ptTOFlowerMin );
        helperpid -> SetisMC( isMC );
        task->SetHelperPID( helperpid );

        if(trigger) task -> SelectCollisionCandidates(AliVEvent::kINT7); //pPb, PbPb_2015
        else task -> SelectCollisionCandidates(AliVEvent::kMB); // PbPb & pp

        cout << "Creating task output container" << endl;

        taskOutputContainer = analysisManager->CreateContainer(listName, TList::Class(), AliAnalysisManager::kOutputContainer,
							     Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname));

        cout << "Add task to analysis manager and connect it to input and output containers" << endl;


        analysisManager->AddTask(task);
        analysisManager->ConnectInput( task,  0, analysisManager->GetCommonInputContainer());
        analysisManager->ConnectOutput(task,  1, taskOutputContainer );


        iTask++;

      }//pt loop
    }//spherocity loop
  }  //centrality loop         
  return task;
}
