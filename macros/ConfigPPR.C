void Config()
{
  //
  // Set Random Number seed
  gRandom->SetSeed(12345);
  
  // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
  new AliGeant3("C++ Interface to Geant3");

  //=======================================================================
  //  Create the output file
   
  TFile *rootfile = new TFile("galice.root","recreate");
  rootfile->SetCompressionLevel(2);
  TGeant3 *geant3 = (TGeant3*)gMC;

  //
  // Set External decayer
  AliDecayer* decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);

  //
  //
  //=======================================================================
  // ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
  geant3->SetTRIG(1); //Number of events to be processed 
  geant3->SetSWIT(4,10);
  geant3->SetDEBU(0,0,1);
  //geant3->SetSWIT(2,2);
  geant3->SetDCAY(1);
  geant3->SetPAIR(1);
  geant3->SetCOMP(1);
  geant3->SetPHOT(1);
  geant3->SetPFIS(0);
  geant3->SetDRAY(0);
  geant3->SetANNI(1);
  geant3->SetBREM(1);
  geant3->SetMUNU(1);
  geant3->SetCKOV(1);
  geant3->SetHADR(1); //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
  geant3->SetLOSS(2);
  geant3->SetMULS(1);
  geant3->SetRAYL(1);
  geant3->SetAUTO(1); //Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
  geant3->SetABAN(0); //Restore 3.16 behaviour for abandoned tracks
  geant3->SetOPTI(2); //Select optimisation level for GEANT geometry searches (0,1,2)
  geant3->SetERAN(5.e-7);

  Float_t cut    = 1.e-3; // 1MeV cut by default
  Float_t tofmax = 1.e10;
  //             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
  geant3->SetCUTS(cut,cut, cut, cut, cut, cut,  cut,  cut, cut,  cut, tofmax);
  //
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV
  AliGenHIJINGpara *gener = new AliGenHIJINGpara(84210);
  gener->SetMomentumRange(0,999);
  gener->SetPhiRange(0,360);
  //  gener->SetThetaRange(0.28,179.72);
  gener->SetThetaRange(0.,180.);
  gener->SetOrigin(0,0,0);        //vertex position
  gener->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
  gener->Init();
  // 
  // Activate this line if you want the vertex smearing to happen
  // track by track
  //
  //gener->SetVertexSmear(perTrack); 

  gAlice->SetField(-999,2);    //Specify maximum magnetic field in Tesla (neg. ==> default field)

  Int_t iABSO=1;
  Int_t iCASTOR=1;
  Int_t iDIPO=1;
  Int_t iFMD=1;
  Int_t iFRAME=1;
  Int_t iHALL=1;
  Int_t iITS=1;
  Int_t iMAG=1;
  Int_t iMUON=1;
  Int_t iPHOS=1;
  Int_t iPIPE=1;
  Int_t iPMD=1;
  Int_t iRICH=1;
  Int_t iSHIL=1;
  Int_t iSTART=1;
  Int_t iTOF=1;
  Int_t iTPC=1;
  Int_t iTRD=1;
  Int_t iZDC=1;

  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY","Alice envelop");

  if(iMAG) {
    //=================== MAG parameters ============================
    // --- Start with Magnet since detector layouts may be depending ---
    // --- on the selected Magnet dimensions ---
    AliMAG *MAG  = new AliMAG("MAG","Magnet");
  }


  if(iABSO) {
    //=================== ABSO parameters ============================
    AliABSO *ABSO  = new AliABSOv0("ABSO","Muon Absorber");
  }

  if(iDIPO) {
    //=================== DIPO parameters ============================

    AliDIPO *DIPO  = new AliDIPOv2("DIPO","Dipole version 2");
  }

  if(iHALL) {
    //=================== HALL parameters ============================

    AliHALL *HALL  = new AliHALL("HALL","Alice Hall");
  }


  if(iFRAME) {
    //=================== FRAME parameters ============================

    AliFRAME *FRAME  = new AliFRAMEv1("FRAME","Space Frame");

  }

  if(iSHIL) {
    //=================== SHIL parameters ============================

    AliSHIL *SHIL  = new AliSHILv0("SHIL","Shielding");
  }


  if(iPIPE) {
    //=================== PIPE parameters ============================

    AliPIPE *PIPE  = new AliPIPEv0("PIPE","Beam Pipe");
  }


  if(iITS) {
    //=================== ITS parameters ============================
    //
    // As the innermost detector in ALICE, the Inner Tracking System "impacts" on
    // almost all other detectors. This involves the fact that the ITS geometry 
    // still has several options to be followed in parallel in order to determine 
    // the best set-up which minimizes the induced background. All the geometries
    // available to date are described in the following. Read carefully the comments 
    // and use the default version (the only one uncommented) unless you are making
    // comparisons and you know what you are doing. In this case just uncomment the
    // ITS geometry you want to use and run Aliroot. 
    //
    // Detailed geometries:
    // ====================
    //
    //
    //AliITS *ITS  = new AliITSv3("ITS","Old ITS detailed version as of the ALICE TP");
    //
    //AliITS *ITS  = new AliITSv5("ITS","Current ITS detailed version used for the ITS TDR");
    //
    //AliITS *ITS  = new AliITSv5symm("ITS","Updated ITS TDR detailed version with symmetric services");
    //
    AliITS *ITS  = new AliITSv5asymm("ITS","Updates ITS TDR detailed version with asymmetric services");
    //
    //
    // Coarse geometries (warning: no hits are produced with these coarse geometries and they unuseful for reconstruction !):
    // ======================================================================================================================
    //
    //
    //AliITS *ITS  = new AliITSv1("ITS","Old ITS coarse version as of the ALICE TP");
    //
    //AliITS *ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS coarse version with asymmetric services");
    //
    //AliITS *ITS  = new AliITSvPPRcoarsesymm("ITS","New ITS coarse version with symmetric services");
    //
    //
    // Geant3 <-> EUCLID conversion
    // ============================
    //
    // SetEUCLID is a flag to output (=1) or not to output (=0) both geometry and 
    // media to two ASCII files (called by default ITSgeometry.euc and 
    // ITSgeometry.tme) in a format understandable to the CAD system EUCLID. 
    // The default (=0) means that you dont want to use this facility.
    //
    ITS->SetEUCLID(0);
  }



  if(iTPC) {
    //============================ TPC parameters ================================
    // --- This allows the user to specify sectors for the SLOW (TPC geometry 2)
    // --- Simulator. SecAL (SecAU) <0 means that ALL lower (upper)
    // --- sectors are specified, any value other than that requires at least one 
    // --- sector (lower or upper)to be specified!
    // --- Reminder: sectors 1-24 are lower sectors (1-12 -> z>0, 13-24 -> z<0)
    // ---           sectors 25-72 are the upper ones (25-48 -> z>0, 49-72 -> z<0)
    // --- SecLows - number of lower sectors specified (up to 6)
    // --- SecUps - number of upper sectors specified (up to 12)
    // --- Sens - sensitive strips for the Slow Simulator !!!
    // --- This does NOT work if all S or L-sectors are specified, i.e.
    // --- if SecAL or SecAU < 0
    //
    //
    //-----------------------------------------------------------------------------

    //  gROOT->LoadMacro("SetTPCParam.C");
    //  AliTPCParam *param = SetTPCParam();
    AliTPC *TPC  = new AliTPCv2("TPC","Default");
    // All sectors included 
    TPC->SetSecAL(-1);
    TPC->SetSecAU(-1);

  }


  if(iTOF) {
    //=================== TOF parameters ============================
    AliTOF *TOF  = new AliTOFv2("TOF","normal TOF");
  }

  if(iRICH) {
    //=================== RICH parameters ===========================
    AliRICH *RICH  = new AliRICHv1("RICH","normal RICH");    

  }


  if(iZDC) {
    //=================== ZDC parameters ============================

    AliZDC *ZDC  = new AliZDCv1("ZDC","normal ZDC");
  }

  if(iCASTOR) {
    //=================== CASTOR parameters ============================

    AliCASTOR *CASTOR  = new AliCASTORv1("CASTOR","normal CASTOR");
  }

  if(iTRD) {
    //=================== TRD parameters ============================
  
    AliTRD *TRD  = new AliTRDv1("TRD","TRD slow simulator");
  
    // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 90% Xe + 10% CO2)
    TRD->SetGasMix(1);
  
    // With hole in front of PHOS
    TRD->SetPHOShole();
    // With hole in front of RICH
    TRD->SetRICHhole();
    // Switch on TR
    AliTRDsim *TRDsim = TRD->CreateTR();
  }

  if(iFMD) {
    //=================== FMD parameters ============================

    AliFMD *FMD  = new AliFMDv0("FMD","normal FMD");
  }

  if(iMUON) {
    //=================== MUON parameters ===========================
    AliMUON *MUON  = new AliMUONv0("MUON","normal MUON");
  }
 
  //=================== PHOS parameters ===========================

  if(iPHOS) {
    AliPHOS *PHOS  = new AliPHOSv1("PHOS","GPS2");
  }


  if(iPMD) {
    //=================== PMD parameters ============================

    AliPMD *PMD  = new AliPMDv1("PMD","normal PMD");
    PMD->SetPAR(1., 1., 0.8, 0.02);
    PMD->SetIN(6., 18., -580., 27., 27.);
    PMD->SetGEO(0.0, 0.2, 4.);
    PMD->SetPadSize(0.8, 1.0, 1.0, 1.5);

  }

  if(iSTART) {
    //=================== START parameters ============================
    AliSTART *START  = new AliSTARTv1("START","START Detector");
  }

         
}
