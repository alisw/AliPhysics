static Int_t    eventsPerRun = 100;
void Config()
{
  cout << "==> Config.C..." << endl;
  
  // Set Random Number seed
  // gRandom->SetSeed(12345);
  
  
  // libraries required by fluka21

  if (!gSystem->Getenv("WITH_ROOT")) {
     cout << "=== RUNNING TFluka with FLUGG ===\n";
     char * gvmc = gSystem->ExpandPathName("$(G4VMC)/examples/macro/g4libs.C");
     gROOT->LoadMacro(gvmc);
     g4libs();

     cout << "\t* Loading Flugg..." << endl;  
     gSystem->Load("libFlugg");    
  } else {
     cout << "=== RUNNING TFluka with TGeo ===\n";
     gSystem->Load("libGeom");
  }     
  cout << "\t* Loading TFluka..." << endl;  
  gSystem->Load("libTFluka");    
    
  cout << "\t* Instantiating TFluka..." << endl;
  new  TFluka("C++ Interface to Fluka", 3/*verbositylevel*/);
  
  cout << "\t* Recreating galice.root if needed..." << endl;
  
  if (!gSystem->Getenv("CONFIG_FILE"))
  {
      cout<<"Config.C: Creating Run Loader ..."<<endl;
      AliRunLoader* rl = AliRunLoader::Open("galice.root",AliConfig::fgkDefaultEventFolderName,
					    "recreate");
      if (rl == 0x0)
      {
	  gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
	  return;
      }
      rl->SetCompressionLevel(2);
      rl->SetNumberOfEventsPerFile(3);
      gAlice->SetRunLoader(rl);
  }

  
  TFluka *fluka = (TFluka *) gMC;
  fluka->SetCoreInputFileName("corealice.inp");
  fluka->SetInputFileName("alice.inp");
  //
  // Set External decayer
  TVirtualMCDecayer *decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
  //
  //
  // Physics process control
  gMC ->SetProcess("DCAY",1);
  gMC ->SetProcess("PAIR",0);
  gMC ->SetProcess("COMP",0);
  gMC ->SetProcess("PHOT",0);
  gMC ->SetProcess("PFIS",0);
  gMC ->SetProcess("DRAY",0);
  gMC ->SetProcess("ANNI",0);
  gMC ->SetProcess("BREM",0);
  gMC ->SetProcess("MUNU",1);
  gMC ->SetProcess("HADR",1); //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
  gMC ->SetProcess("LOSS",2);
  gMC ->SetProcess("MULS",1);
  //xx gMC ->SetProcess("RAYL",1);

  // Energy cuts
  // (in development)
  Float_t cut    = 1.e-1; // 100 MeV cut by default
  Float_t tofmax = 1.e10; 

  gMC ->SetCut("CUTGAM",cut);
  gMC ->SetCut("CUTELE",cut);
  gMC ->SetCut("CUTNEU",cut);
  gMC ->SetCut("CUTHAD",cut);
  gMC ->SetCut("CUTMUO",cut);
  gMC ->SetCut("BCUTE",cut);
  gMC ->SetCut("BCUTM",cut);
  gMC ->SetCut("DCUTE",cut);
  gMC ->SetCut("DCUTM",cut);
  gMC ->SetCut("PPCUTM",cut);
  gMC ->SetCut("TOFMAX",tofmax);

  //
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV
  if (gSystem->Getenv("CONFIG_NPARTICLES"))
      int     nParticles = atoi(gSystem->Getenv("CONFIG_NPARTICLES"));
  else
      int     nParticles = 10;
  
  cout << "\t* Creating and configuring generator for " << nParticles 
       << " particles..." << endl;
  
  AliGenHIJINGpara *gener = new AliGenHIJINGpara(nParticles);
  
  gener->SetMomentumRange(0, 999);
  gener->SetPhiRange(0, 360);
  // Set pseudorapidity range from -8 to 8.
  Float_t thmin = EtaToTheta( 0.1);   // theta min. <---> eta max
  Float_t thmax = EtaToTheta(-0.1);  // theta max. <---> eta min 
  gener->SetThetaRange(thmin,thmax);
  gener->SetOrigin(0, 0, 0);  //vertex position
  gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position
  gener->Init();
  // 
  // Activate this line if you want the vertex smearing to happen
  // track by track
  //

  gAlice->SetField(-999, 2);  //Specify maximum magnetic field in Tesla (neg. ==> default field)
  gAlice->SetDebug(10);
  
  Int_t   iABSO  = 0; //1
  Int_t   iCRT   = 0; //Not good ?
  Int_t   iDIPO  = 0; //1
  Int_t   iFMD   = 0; //1
  Int_t   iFRAME = 0; //1
  Int_t   iHALL  = 0; //1
  Int_t   iITS   = 0; //1
  Int_t   iMAG   = 0; //1
  Int_t   iMUON  = 0; //1. Not good (newFlagLttc=10000 is outside array bounds)
  Int_t   iPHOS  = 0; //1
  Int_t   iPIPE  = 0; //1
  Int_t   iPMD   = 0; //Not good (too many regions)
  Int_t   iRICH  = 1; //1. Not good (no tracking with FRAME)
  Int_t   iSHIL  = 0; //1. Not good (no tracking) (it works alone)
  Int_t   iSTART = 0; //1. Not good (no tracking) (it works alone)
  Int_t   iTOF   = 0; //1. Not good (no tracking) (newFlagLttc=10000 is outside array bounds if alone)
  Int_t   iTPC   = 1;
  Int_t   iTRD   = 0; //1. Not good (no tracking) (Crash alone with FRAME)
  Int_t   iZDC   = 0; //1. Needs SHIL and others
  Int_t   iEMCAL = 0; //Not good (Crash)
  Int_t   iVZERO = 0;
 
  cout << "\t* Creating the detectors ..." << endl;
  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY", "Alice envelop");
  cout << "\t\t+ BODY..." << endl;
  
  
  if (iMAG)
    {
      //=================== MAG parameters ============================
      // --- Start with Magnet since detector layouts may be depending ---
      // --- on the selected Magnet dimensions ---
      cout << "\t\t+ Magnet..." << endl;
      AliMAG *MAG = new AliMAG("MAG", "Magnet");
    }
  
  
  if (iABSO)
    {
      //=================== ABSO parameters ============================
      cout << "\t\t+ ABSO..." << endl;
      AliABSO *ABSO = new AliABSOv0("ABSO", "Muon Absorber");
    }
  
  if (iDIPO)
    {
      //=================== DIPO parameters ============================
      cout << "\t\t+ DIPO..." << endl;     
      AliDIPO *DIPO = new AliDIPOv2("DIPO", "Dipole version 2");
    }
  
  if (iHALL)
    {
      //=================== HALL parameters ============================
      cout << "\t\t+ HALL..." << endl;
      AliHALL *HALL = new AliHALL("HALL", "Alice Hall");
    }
  
  
  if (iFRAME)
    {
      //=================== FRAME parameters ============================
      
      cout << "\t\t+ FRAME..." << endl;
      AliFRAME *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
      
    }
  
  if (iSHIL)
    {
      //=================== SHIL parameters ============================
      
      cout << "\t\t+ SHIL..." << endl;
      AliSHIL *SHIL = new AliSHILv2("SHIL", "Shielding");
    }
  
  
  if (iPIPE)
    {
      //=================== PIPE parameters ============================
      
      cout << "\t\t+ PIPE..." << endl;
      AliPIPE *PIPE = new AliPIPEv0("PIPE", "Beam Pipe");
    }
  
  if(iITS) {
      cout << "\t\t+ ITS..." << endl;
    
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
    //
    //
    //AliITS *ITS  = new AliITSv5symm("ITS","Updated ITS TDR detailed version with symmetric services");
    //
    //AliITS *ITS  = new AliITSv5asymm("ITS","Updates ITS TDR detailed version with asymmetric services");
    //
    AliITSvPPRasymm *ITS  = new AliITSvPPRasymm("ITS","New ITS PPR detailed version with asymmetric services");
    ITS->SetMinorVersion(2);					 // don't touch this parameter if you're not an ITS developer
    ITS->SetReadDet(kFALSE);					 // don't touch this parameter if you're not an ITS developer
    //    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
    ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
    ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
    ITS->SetThicknessChip1(200.);  // chip thickness on layer 1 must be in the range [150,300]
    ITS->SetThicknessChip2(200.);  // chip thickness on layer 2 must be in the range [150,300]
    ITS->SetRails(1);	     // 1 --> rails in ; 0 --> rails out
    ITS->SetCoolingFluid(1);   // 1 --> water ; 0 --> freon
    //
    //AliITSvPPRsymm *ITS  = new AliITSvPPRsymm("ITS","New ITS PPR detailed version with symmetric services");
    //ITS->SetMinorVersion(2);                                       // don't touch this parameter if you're not an ITS developer
    //ITS->SetReadDet(kFALSE);                                       // don't touch this parameter if you're not an ITS developer
    //ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRsymm2.det"); // don't touch this parameter if you're not an ITS developer
    //ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
    //ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
    //ITS->SetThicknessChip1(200.);  // chip thickness on layer 1 must be in the range [150,300]
    //ITS->SetThicknessChip2(200.);  // chip thickness on layer 2 must be in the range [150,300]
    //ITS->SetRails(1);              // 1 --> rails in ; 0 --> rails out
    //ITS->SetCoolingFluid(1);       // 1 --> water ; 0 --> freon
    //
    //
    // Coarse geometries (warning: no hits are produced with these coarse geometries and they unuseful 
    // for reconstruction !):
    //                                                     
    //
    //AliITSvPPRcoarseasymm *ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS PPR coarse version with asymmetric services");
    //ITS->SetRails(1);                // 1 --> rails in ; 0 --> rails out
    //ITS->SetSupportMaterial(0);      // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    //
    //AliITS *ITS  = new AliITSvPPRcoarsesymm("ITS","New ITS PPR coarse version with symmetric services");
    //ITS->SetRails(1);                // 1 --> rails in ; 0 --> rails out
    //ITS->SetSupportMaterial(0);      // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    //                      
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
  
  
  if (iTPC)
    {
      cout << "\t\t+ TPC..." << endl;
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
      //
      //  gROOT->LoadMacro("SetTPCParam.C");
      //  AliTPCParam *param = SetTPCParam();
      AliTPC *TPC = new AliTPCv1("TPC", "Default");
      // All sectors included 
      TPC->SetSecAL(-1);
      TPC->SetSecAU(-1);
    }
  
  if (iTOF)
    {
      cout << "\t\t+ TOF..." << endl;
      //=================== TOF parameters ============================
      AliTOF *TOF = new AliTOFv2("TOF", "normal TOF");
    }
  
  if (iRICH)
    {
      cout << "\t\t+ RICH..." << endl;
      //=================== RICH parameters ===========================
        AliRICH *RICH = new AliRICHv1("RICH", "normal RICH");
	
    }
  
  
  if (iZDC)
    {
      cout << "\t\t+ ZDC..." << endl;
      //=================== ZDC parameters ============================
      
      AliZDC *ZDC = new AliZDCv2("ZDC", "normal ZDC");
    }
  
  if (iCRT)
    {
      cout << "\t\t+ CRT..." << endl;
      //=================== CRT parameters ============================
      
      AliCRT *CRT = new AliCRTv0("CRT", "normal CRT");
    }
  
  if (iTRD)
    {
      cout << "\t\t+ TRD..." << endl;
      //=================== TRD parameters ============================
      
      AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
      
      // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 90% Xe + 10% CO2)
      TRD->SetGasMix(1);
      
      // With hole in front of PHOS
      TRD->SetPHOShole();
      // With hole in front of RICH
      TRD->SetRICHhole();
      // Switch on TR
      AliTRDsim *TRDsim = TRD->CreateTR();
    }
  
  if (iFMD)
    {
      cout << "\t\t+ FMD..." << endl;
      //=================== FMD parameters ============================
      
      AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
      FMD->SetRingsSi1(256);
      FMD->SetRingsSi2(64);
      FMD->SetSectorsSi1(20);
      FMD->SetSectorsSi2(24);
    }
  
  if (iMUON)
    {
      cout << "\t\t+ MUON..." << endl;
      //=================== MUON parameters ===========================
      
      AliMUON *MUON = new AliMUONv1("MUON", "default");
    }
  //=================== PHOS parameters ===========================
  
  if (iPHOS)
    {
      cout << "\t\t+ PHOS..." << endl;
      AliPHOS *PHOS = new AliPHOSv1("PHOS", "GPS2");
    }
  
  
  if (iPMD)
    {
      cout << "\t\t+ PMD..." << endl;
      //=================== PMD parameters ============================
      
      AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");
      
      PMD->SetPAR(1., 1., 0.8, 0.02);
      PMD->SetIN(6., 18., -580., 27., 27.);
      PMD->SetGEO(0.0, 0.2, 4.);
      PMD->SetPadSize(0.8, 1.0, 1.0, 1.5);
      
    }
  
  if (iEMCAL && !iRICH)
    {
      cout << "\t\t+ EMCAL (no RICH)..." << endl;
      //=================== EMCAL parameters ============================
      AliEMCAL *EMCAL = new AliEMCALv1("EMCAL", "EMCALArch1a");
    }
  
  if (iSTART)
    {
      cout << "\t\t+ START..." << endl;
      //=================== START parameters ============================
      AliSTART *START = new AliSTARTv1("START", "START Detector");
    }
  if (iVZERO)
    {
      cout << "\t\t+ VZERO..." << endl;
      //=================== CRT parameters ============================
      AliVZERO *VZERO = new AliVZEROv2("VZERO", "normal VZERO");
    }
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
