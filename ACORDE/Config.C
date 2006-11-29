static Int_t    eventsPerRun = 100;
enum PprGeo_t 
{
    kHoles, kNoHoles
};
static PprGeo_t geo = kHoles;

void Config()
{
    // 7-DEC-2000 09:00
    // Switch on Transition Radiation simulation. 6/12/00 18:00
    // iZDC=1  7/12/00 09:00
    // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
    // Theta range given through pseudorapidity limits 22/6/2001




    // Set Random Number seed
    // gRandom->SetSeed(12345);
    TDatime dt;
    UInt_t curtime=dt.Get();
    UInt_t procid=gSystem->GetPid();
    UInt_t seed=curtime-procid;
    gRandom->SetSeed(seed);
    cerr<<"Seed for random number generation= "<<seed<<endl; 
    cout<<"Seed for random number generation= "<<seed<<endl;


    new     AliGeant3("C++ Interface to Geant3");

    if (!gSystem->Getenv("CONFIG_FILE"))
    {
        TFile  *rootfile = new TFile("galice.root", "recreate");

        rootfile->SetCompressionLevel(2);
    }

    TGeant3 *geant3 = (TGeant3 *) gMC;

    //
    // Set External decayer
    AliDecayer *decayer = new AliDecayerPythia();

    decayer->SetForceDecay(kAll);
    decayer->Init();
    gMC->SetExternalDecayer(decayer);
    //
    //
    //=======================================================================
    // ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
    geant3->SetTRIG(1);         //Number of events to be processed 
    geant3->SetSWIT(4, 10);
    geant3->SetDEBU(0, 0, 1);
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
    geant3->SetHADR(4);         //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
    geant3->SetLOSS(2);
    geant3->SetMULS(1);
    geant3->SetRAYL(1);
    geant3->SetAUTO(1);         //Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
    geant3->SetABAN(0);         //Restore 3.16 behaviour for abandoned tracks
    geant3->SetOPTI(2);         //Select optimisation level for GEANT geometry searches (0,1,2)
    geant3->SetERAN(5.e-7);

    Float_t cut = 1.e-3;        // 1MeV cut by default
    Float_t tofmax = 1.e10;

    //             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
    geant3->SetCUTS(cut, cut, cut, cut, cut, cut, cut, cut, cut, cut,
                    tofmax);
    //
    //=======================================================================
    // ************* STEERING parameters FOR ALICE SIMULATION **************
    // --- Specify event type to be tracked through the ALICE setup
    // --- All positions are in cm, angles in degrees, and P and E in GeV
    if (gSystem->Getenv("CONFIG_NPARTICLES"))
    {
        int     nParticles = atoi(gSystem->Getenv("CONFIG_NPARTICLES"));
    } else
    {
        int     nParticles = 50;
    }

    AliGenACORDE *gener = new AliGenACORDE(); // One generation per event
    //TF1* zenithFunc = new TF1("zenith","1/cos(x)", 0. 65.);
    //gener->SetZenithalDistributionFunction(zenithFunc);
    //Float_t tend = 65.*TMath::Pi()/180.;
    //TF1* momFunc = new TF1("momFunc", "x", 10., 1000.);
    //gener->SetMomentumDistrubutionFunction(momFunc);
    gener->SetZenithalAngleRange(0., 65.);
    gener->SetAzimuthalAngleRange(0., 359.);
    gener->SetMomentumResolution(1.); // GeVs
    gener->SetSigma(0,0,0); // Sigma for smearing
    gener->SetMomentumRange(10., 1000.); // in GeVs.

    gener->SetPart(kMuonMinus);
    // Two choices:
    //           kMuonMinus
    //           kMuonPlus

    //gener->SetMode(kSingleMuons);// Single muons
    //gener->SetMode(kMuonBundle);
    gener->SetMode(kMuonFlux);
    // three choises
    //          kSingleMuons
    //          kMuonBundle
    //          kMuonFlux  <- Not available, yet

    // If you choose kMuonBundle, you should set how many particles do you
    // want and in wich region should the muons be generated, for instance
    // if you want to generate n^2 muons, in a square of side 300 cm
    //gener->SetRange(4, 300., 4, 300.);

    //gener->SetVertexSmear(kPerEvent); // Set the smearing per Event

    gener->Init(); // Intialize the generator.

    // Set the tracking limits to include the Molasse
    gAlice->TrackingLimits( 5000, 5000);

    // 
    // Activate this line if you want the vertex smearing to happen
    // track by track
    //
    //gener->SetVertexSmear(perTrack); 
    // Field (L3 0.4 T)
    AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
    rootfile->cd();
    gAlice->SetField(field);    


    Int_t   iABSO  =  0;
    Int_t   iDIPO  =  0;
    Int_t   iFMD   =  0;
    Int_t   iFRAME =  0;
    Int_t   iHALL  =  0;
    Int_t   iITS   =  0;
    Int_t   iMAG   =  0;
    Int_t   iMUON  =  0;
    Int_t   iPHOS  =  0;
    Int_t   iPIPE  =  0;
    Int_t   iPMD   =  0;
    Int_t   iHMPID  =  0;
    Int_t   iSHIL  =  0;
    Int_t   iT0 =  0;
    Int_t   iTOF   =  0;
    Int_t   iTPC   =  0;
    Int_t   iTRD   =  0;
    Int_t   iZDC   =  0;
    Int_t   iEMCAL =  0;
    Int_t   iACORDE   =  1;
    Int_t   iVZERO =  0;

    //=================== Alice BODY parameters =============================
    AliBODY *BODY = new AliBODY("BODY", "Alice envelop");

    if (iMAG)
    {
        //=================== MAG parameters ============================
        // --- Start with Magnet since detector layouts may be depending ---
        // --- on the selected Magnet dimensions ---
        AliMAG *MAG = new AliMAG("MAG", "Magnet");
    }


    if (iABSO)
    {
        //=================== ABSO parameters ============================
        AliABSO *ABSO = new AliABSOv0("ABSO", "Muon Absorber");
    }

    if (iDIPO)
    {
        //=================== DIPO parameters ============================

        AliDIPO *DIPO = new AliDIPOv2("DIPO", "Dipole version 2");
    }

    if (iHALL)
    {
        //=================== HALL parameters ============================

        AliHALL *HALL = new AliHALL("HALL", "Alice Hall");
    }


    if (iFRAME)
    {
        //=================== FRAME parameters ============================

        AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
	if (geo == kHoles) {
	    FRAME->SetHoles(1);
	} else {
	    FRAME->SetHoles(0);
	}
    }

    if (iSHIL)
    {
        //=================== SHIL parameters ============================

        AliSHIL *SHIL = new AliSHILv2("SHIL", "Shielding Version 2");
    }


    if (iPIPE)
    {
        //=================== PIPE parameters ============================

        AliPIPE *PIPE = new AliPIPEv0("PIPE", "Beam Pipe");
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
	ITS->SetRails(0);	     // 1 --> rails in ; 0 --> rails out
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
    //ITS->SetRails(0);              // 1 --> rails in ; 0 --> rails out
    //ITS->SetCoolingFluid(1);       // 1 --> water ; 0 --> freon
    //
    //
    // Coarse geometries (warning: no hits are produced with these coarse geometries and they unuseful 
    // for reconstruction !):
    //                                                     
    //
    //AliITSvPPRcoarseasymm *ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS PPR coarse version with asymmetric services");
    //ITS->SetRails(0);                // 1 --> rails in ; 0 --> rails out
    //ITS->SetSupportMaterial(0);      // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    //
    //AliITS *ITS  = new AliITSvPPRcoarsesymm("ITS","New ITS PPR coarse version with symmetric services");
    //ITS->SetRails(0);                // 1 --> rails in ; 0 --> rails out
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
        AliTPC *TPC = new AliTPCv2("TPC", "Default");

        // All sectors included 
        TPC->SetSecAL(-1);
        TPC->SetSecAU(-1);

    }


    if (iTOF) {
	if (geo == kHoles) {
        //=================== TOF parameters ============================
	    AliTOF *TOF = new AliTOFv2FHoles("TOF", "TOF with Holes");
	} else {
	    AliTOF *TOF = new AliTOFv4T0("TOF", "normal TOF");
	}
    }


    if (iHMPID)
    {
        //=================== HMPID parameters ===========================
        AliHMPID *HMPID = new AliHMPIDv3("HMPID", "normal HMPID");

    }


    if (iZDC)
    {
        //=================== ZDC parameters ============================

        AliZDC *ZDC = new AliZDCv2("ZDC", "normal ZDC");
    }

    if (iTRD)
    {
        //=================== TRD parameters ============================

        AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");

        // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 90% Xe + 10% CO2)
        TRD->SetGasMix(1);
	if (geo == kHoles) {
	    // With hole in front of PHOS
	    TRD->SetPHOShole();
	    // With hole in front of HMPID
	    TRD->SetHMPIDhole();
	}
	    // Switch on TR
	    AliTRDsim *TRDsim = TRD->CreateTR();
    }

    if (iFMD)
    {
        //=================== FMD parameters ============================
	AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
        FMD->SetRingsSi1(256);
        FMD->SetRingsSi2(128);
        FMD->SetSectorsSi1(20);
        FMD->SetSectorsSi2(40);      
   }

    if (iMUON)
    {
        //=================== MUON parameters ===========================

        AliMUON *MUON = new AliMUONv1("MUON", "default");
    }
    //=================== PHOS parameters ===========================

    if (iPHOS)
    {
        AliPHOS *PHOS = new AliPHOSv1("PHOS", "IHEP");
    }


    if (iPMD)
    {
        //=================== PMD parameters ============================
        AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");
    }

    if (iT0)
    {
        //=================== T0 parameters ============================
        AliT0 *T0 = new AliT0v1("T0", "T0 Detector");
    }

    if (iEMCAL)
    {
        //=================== EMCAL parameters ============================
        AliEMCAL *EMCAL = new AliEMCALv1("EMCAL", "EMCALArch1a");
    }

     if (iACORDE)
    {
        //=================== ACORDE parameters ============================
        //AliACORDE *ACORDE = new AliACORDEv0("ACORDE", "normal ACORDE");
	AliACORDE* ACORDE = new AliACORDEv1("ACORDE", "Special version Only with the ACORDE");
    }

     if (iVZERO)
    {
        //=================== ACORDE parameters ============================
        AliVZERO *VZERO = new AliVZEROv2("VZERO", "normal VZERO");
    }

}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
