//  Macro designed for use with the AliAnalysisTaskR2P2wSphero task.
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
// pp18_V0_kMB_kFALSE   18:     centralityMethod = 4 (V0),        trigger = kFALSE (AliVEvent::kMB).
/////////////////////////////////////////////////////////////////////////////////


AliAnalysisTaskR2P2wSphero * AddTaskR2P2wSphero
(
 TString AnalysisDataType       = "RealData", // "RealData"; "MCAOD" for MC AOD truth; "MCAODreco"
 TString prefixName             = "data_",
  //************ for Stage2_Step1 starts ************************
 
 int    singlesOnly             =  0,   // 0: full correlations    1: singles only
 
 int    usePtEff                =  1,   // 0: no                   1: yes
 TString inputPtEffFileName     = "alien:///alice/cern.ch/user/s/subhadee/inv_efficiency_correction_factor.root",
  
 //int    chargeSet               =  1,   // 0: ++    1: +-    2: -+    3: --
 //const char* taskname           = "ChPM",// ChPM, ChPP, ChMM
 //************ for Stage2_Step1 ends ************************
 
 TString System                 = "pp18_V0_kMB_kFALSE",
 int    CentralityGroup         =  2,  // Diff Cent Groups dealing w/ memory limit & weight file 100M Alien limit
 int    minNch                  =  5, //Min number of charge particles for spherocity calculation
 int    SpherocityGroup         =  13,   //spherocity groups
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
 double ptMin                   =  0.2, // pt range lower limit cut ( also for pt histos )
 double ptCUTupperMax           =  2.0, // pt range upper limit cut
 double ptMax                   =  3.0, // pt range upper limit for histos; NOT pt cut!!!
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

  int    chargeSet[4]; // 0: ++    1: +-    2: -+    3: --
  const char* taskname[4];
  int nChSet;
  
  if(singlesOnly)
  {
    nChSet = 1;
    chargeSet[0] = 1;  taskname[0] = "ChPM";
  }
  else 
  {
    nChSet = 3;
    chargeSet[0] = 1;  taskname[0] = "ChPM";
    chargeSet[1] = 0;  taskname[1] = "ChPP";                                 
    chargeSet[2] = 3;  taskname[2] = "ChMM";                                 
       
  }

  double minSherocity[15];
  double maxSherocity[15];
  int    nSherocity;

//////////////////////////////////////////**************** spherocity values calculated with Nch>=5 & 0-100% mult class**************///////////////////////

 if ( SpherocityGroup ==1) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.435; //0- 20%
    minSherocity[1] = 0.755;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==2) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.455; //0- 20%
    minSherocity[1] = 0.765;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==3) //corresponds to the spherocity for  Data/Reco or Gen level
  {
    nSherocity=1;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.985; //0- 100%
  }

  ///////////////////////////////////******************** spherocity values calculated with Nch>=10 & 0-100% mult class**************************//////////////////////

 if ( SpherocityGroup ==4) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.525; //0- 20%
    minSherocity[1] = 0.795;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==5) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.545; //0- 20%
    minSherocity[1] = 0.805;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==6) //corresponds to the spherocity for  Data/Reco or Gen level
  {
    nSherocity=1;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.985; //0- 100%
  }

  //////////////////////////////////////////////////********* spherocity values calculated with Nch>=10 & 0-10% mult class**************//////////////////

 if ( SpherocityGroup ==7) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.595; //0- 20%
    minSherocity[1] = 0.825;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==8) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.635; //0- 20%
    minSherocity[1] = 0.845;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==9) //corresponds to the spherocity for  Data/Reco or Gen level
  {
    nSherocity=1;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.985; //0- 100%
  }

//////////////////////////////////////////////////********* spherocity values calculated with Nch>=10 & 10-30% mult class**************//////////////////

 if ( SpherocityGroup ==21) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.545; //0- 20%
    minSherocity[1] = 0.795;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==22) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.575; //0- 20%
    minSherocity[1] = 0.815;       maxSherocity[1]  =0.985;  //80-100%
    
  }

//////////////////////////////////////////////////********* spherocity values calculated with Nch>=10 & 30-60% mult class**************//////////////////

 if ( SpherocityGroup ==23) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.495; //0- 20%
    minSherocity[1] = 0.765;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==24) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.515; //0- 20%
    minSherocity[1] = 0.785;       maxSherocity[1]  =0.985;  //80-100%
    
  }

//////////////////////////////////////////////////********* spherocity values calculated with Nch>=10 & 60-100% mult class**************//////////////////

 if ( SpherocityGroup ==25) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.455; //0- 20%
    minSherocity[1] = 0.745;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==26) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.465; //0- 20%
    minSherocity[1] = 0.755;       maxSherocity[1]  =0.985;  //80-100%
    
  }
////////////////////////////////////**************** spherocity values calculated with Nch>=10 & 0-5% mult class **********************//////////////////////////////

 if ( SpherocityGroup ==10) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.615; //0- 20%
    minSherocity[1] = 0.825;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==11) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.655; //0- 20%
    minSherocity[1] = 0.845;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==12) //corresponds to the spherocity for  Data/Reco or Gen level
  {
    nSherocity=1;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.985; //0- 100%

  }


//**************************//////spherocity values calculated with Nch>=5 & for 4 different Multiplicity classes/////////***************************************//

///////////////////////// Mult 0-10% //////////////////////////

 if ( SpherocityGroup ==13) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.595; //0- 20%
    minSherocity[1] = 0.825;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  if ( SpherocityGroup ==14) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.635; //0- 20%
    minSherocity[1] = 0.845;       maxSherocity[1]  =0.985;  //80-100%
    
  }

  
  ///////////////////////// Mult 10-30% //////////////////////////

 if ( SpherocityGroup ==15) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.515; //0- 20%
    minSherocity[1] = 0.785;       maxSherocity[1]  =0.975;  //80-100%
    
  }

  if ( SpherocityGroup ==16) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.565; //0- 20%
    minSherocity[1] = 0.805;       maxSherocity[1]  =0.975;  //80-100%
    
  }

  ///////////////////////// Mult 30-60% //////////////////////////

 if ( SpherocityGroup ==17) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.425; //0- 20%
    minSherocity[1] = 0.735;       maxSherocity[1]  =0.975;  //80-100%
    
  }

  if ( SpherocityGroup ==18) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.455; //0- 20%
    minSherocity[1] = 0.755;       maxSherocity[1]  =0.985;  //80-100%
    
  }

///////////////////////// Mult 60-100% //////////////////////////

 if ( SpherocityGroup ==19) //corresponds to the spherocity for Data/Reco level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.345; //0- 20%
    minSherocity[1] = 0.685;       maxSherocity[1]  =0.965;  //80-100%
    
  }

  if ( SpherocityGroup ==20) //corresponds to the spherocity for Gen level
  {
    nSherocity=2;
    minSherocity[0] = 0.005;       maxSherocity[0]  =0.365; //0- 20%
    minSherocity[1] = 0.695;       maxSherocity[1]  =0.975;  //80-100%
    
  }



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  double minCentrality[15];
  double maxCentrality[15];
  int    nCentrality;

  if ( CentralityGroup == 1 ) 
  {
    nCentrality = 1;
    minCentrality[0] = 0.;	  maxCentrality[0]  = 100.;
  }
  
  else if ( CentralityGroup == 2 ) 
  {
    nCentrality = 1;
    minCentrality[0] = 0.;	  maxCentrality[0]  = 10.;
  }
  
  else if ( CentralityGroup == 3 )
  {
    nCentrality = 4;
    minCentrality[0] = 0;       maxCentrality[0]  = 10.;
    minCentrality[1] = 10.;     maxCentrality[1]  = 30.;
    minCentrality[2] = 30.;     maxCentrality[2]  = 60.;
    minCentrality[3] = 60.;     maxCentrality[3]  = 100.;
  }
  else if ( CentralityGroup == 4 ) 
  {
    nCentrality = 1;
    minCentrality[0] = 0.;	  maxCentrality[0]  = 10.;
  }

  else if ( CentralityGroup == 5 ) 
  {
    nCentrality = 1;
    minCentrality[0] = 10.;	  maxCentrality[0]  = 30.;
  }

  else if ( CentralityGroup == 6 ) 
  {
    nCentrality = 1;
    minCentrality[0] = 30.;	  maxCentrality[0]  = 60.;
  }

  else if ( CentralityGroup == 7 ) 
  {
    nCentrality = 1;
    minCentrality[0] = 60.;	  maxCentrality[0]  = 100.;
  }

  else if ( CentralityGroup == 8 )
  {
    nCentrality = 5;
    minCentrality[0] = 0;       maxCentrality[0]  = 10.;
    minCentrality[1] = 10.;     maxCentrality[1]  = 20.;
    minCentrality[2] = 20.;     maxCentrality[2]  = 40.;
    minCentrality[3] = 40.;     maxCentrality[3]  = 70.;
    minCentrality[4] = 70.;     maxCentrality[4]  = 100.;
  }

  else if ( CentralityGroup == 9 )
  {
    nCentrality = 1;
    minCentrality[0] = 0;       maxCentrality[0]  = 10.;
  }

  else if ( CentralityGroup == 10 )
  {
    nCentrality = 1;
    minCentrality[0] = 10;       maxCentrality[0]  = 20.;
  }

  else if ( CentralityGroup == 11 )
  {
    nCentrality = 1;
    minCentrality[0] = 20;       maxCentrality[0]  = 40.;
  }

  else if ( CentralityGroup == 12 )
  {
    nCentrality = 1;
    minCentrality[0] = 40;       maxCentrality[0]  = 70.;
  }

  else if ( CentralityGroup == 13 )
  {
    nCentrality = 1;
    minCentrality[0] = 70;       maxCentrality[0]  = 100.;
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
      ::Error("AddTaskR2P2wSphero", "No analysis manager to connect to.");
      return NULL;
    }
    
  TString part1Name;
  TString part2Name;
  TString partName;
  TString eventName;
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
  AliAnalysisTaskR2P2wSphero* task;

  for (int iChSet=0; iChSet < nChSet; ++iChSet)
  {
    
  for (int iCentrality=0; iCentrality < nCentrality; ++iCentrality)
    {
            
      for (int isphero = 0; isphero < nSherocity; ++isphero)
      {
      switch ( chargeSet[iChSet] )
        {
	case 0: part1Name = "P"; part2Name = "P"; requestedCharge1 =  1; requestedCharge2 =  1; sameFilter = 1;   break;
	case 1: part1Name = "P"; part2Name = "M"; requestedCharge1 =  1; requestedCharge2 = -1; sameFilter = 0;   break;
	case 2: part1Name = "M"; part2Name = "P"; requestedCharge1 = -1; requestedCharge2 =  1; sameFilter = 0;   break;
	case 3: part1Name = "M"; part2Name = "M"; requestedCharge1 = -1; requestedCharge2 = -1; sameFilter = 1;   break;
        }
        
      partName = "pt";
      partName += int(10*ptMin);
      partName += "to";
      partName += int(10*ptCUTupperMax);

      part1Name += part2Name;
      part1Name += "_";
      part1Name += "pt";
      part1Name += int(10*ptMin);
      part1Name += "to";
      part1Name += int(10*ptCUTupperMax);

      eventName =  "_";
      eventName += int(1000*minSherocity[isphero] );
      eventName += "So";
      eventName += int(1000*maxSherocity[isphero] );


      eventName1 =  "_";
      eventName1 += int(minCentrality[iCentrality] );
      eventName1 += "Vo";
      eventName1 += int(maxCentrality[iCentrality] );
            
      baseName     =   prefixName;      
      baseName     +=  part1Name;
      //baseName     +=  "_";
      //   baseName     +=  part2Name;
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
		    }
		}
	      //------pt-efficiency ends ---------------


     task = new  AliAnalysisTaskR2P2wSphero(taskName);
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
      task->SetNchMin    (          minNch          );
      task->SetVertexZWidth(        vZwidth         );
      task->SetEtaWidth(        etaBinWidth         );
      task->SetVertexXYMin(         -1.             );
      task->SetVertexXYMax(          1.             );
      task->SetCentralityMethod(    centralityMethod);
      task->SetCentrality(          minCentrality[iCentrality], maxCentrality[iCentrality]);
      task->SetTSpherocity(          minSherocity[isphero], maxSherocity[isphero]);
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
      task->SetPtCutUpperLimit( ptCUTupperMax );
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
        
      taskOutputContainer = analysisManager->CreateContainer(listName,
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname[iChSet]));

      cout << "Add task to analysis manager and connect it to input and output containers" << endl;


      analysisManager->AddTask(task);
      analysisManager->ConnectInput( task,  0, analysisManager->GetCommonInputContainer());
      analysisManager->ConnectOutput(task,  1, taskOutputContainer );
      //analysisManager->ConnectOutput(task,  0, taskOutputContainer );
        
      iTask++;
      }//spherocity loop
    }  //centrality      
  }//ichset    
  return task;
}
