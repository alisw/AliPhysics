enum PprRun_t 
{
    test50,
    kParam_8000,   kParam_4000,  kParam_2000, 
    kHijing_cent1, kHijing_cent2, 
    kHijing_per1,  kHijing_per2, kHijing_per3, kHijing_per4,  kHijing_per5,
    kHijing_jj25,  kHijing_jj50, kHijing_jj75, kHijing_jj100, kHijing_jj200, 
    kHijing_gj25,  kHijing_gj50, kHijing_gj75, kHijing_gj100, kHijing_gj200
};

enum PprGeo_t 
{
    kHoles, kNoHoles
};

enum PprRad_t
{
    kGluonRadiation, kNoGluonRadiation
};

enum PprMag_t
{
    k2kG, k4kG, k5kG
};


// This part for configuration    
static PprRun_t run = test50;
static PprGeo_t geo = kHoles;
static PprRad_t rad = kGluonRadiation;
static PprMag_t mag = k4kG;

// Comment line 
static TString  comment;



void Config()
{

    // 7-DEC-2000 09:00
    // Switch on Transition adiation simulation. 6/12/00 18:00
    // iZDC=1  7/12/00 09:00
    // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
    // Theta range given through pseudorapidity limits 22/6/2001

    // Set Random Number seed
    // gRandom->SetSeed(12345);

    new AliGeant3("C++ Interface to Geant3");

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
    geant3->SetHADR(1);         //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
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


// Generator Configuration
    gAlice->SetDebug(1);
    AliGenerator* gener = GeneratorFactory(run);
    gener->SetOrigin(0, 0, 0);    // vertex position
    gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position
    gener->SetCutVertexZ(1.);     // Truncate at 1 sigma
    gener->SetVertexSmear(kPerEvent); 
    gener->SetTrackingFlag(1);
    gener->Init();
    
    if (mag == k2kG) {
	comment = comment.Append(" | L3 field 0.2 T");
    } else if (mag == k4kG) {
	comment = comment.Append(" | L3 field 0.4 T");
    } else if (mag == k5kG) {
	comment = comment.Append(" | L3 field 0.5 T");
    }
    
    
    if (rad == kGluonRadiation)
    {
	comment = comment.Append(" | Gluon Radiation On");
	
    } else {
	comment = comment.Append(" | Gluon Radiation Off");
    }

    if (geo == kHoles)
    {
	comment = comment.Append(" | Holes for PHOS/RICH");
	
    } else {
	comment = comment.Append(" | No holes for PHOS/RICH");
    }

    printf("\n \n Comment: %s \n \n", (char*) comment);
    
    
// Field (L3 0.4 T)
    AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., mag);
    rootfile->cd();
    gAlice->SetField(field);    
    
//
    Int_t   iABSO   = 1;
    Int_t   iDIPO   = 1;
    Int_t   iFMD    = 1;
    Int_t   iFRAME  = 1;
    Int_t   iHALL   = 1;
    Int_t   iITS    = 1;
    Int_t   iMAG    = 1;
    Int_t   iMUON   = 1;
    Int_t   iPHOS   = 1;
    Int_t   iPIPE   = 1;
    Int_t   iPMD    = 1;
    Int_t   iRICH   = 1;
    Int_t   iSHIL   = 1;
    Int_t   iSTART  = 1;
    Int_t   iTOF    = 1;
    Int_t   iTPC    = 1;
    Int_t   iTRD    = 1;
    Int_t   iZDC    = 1;
    Int_t   iEMCAL  = 1;
    Int_t   iVZERO  = 1;
    Int_t   iCRT    = 0;

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


    if (iRICH)
    {
        //=================== RICH parameters ===========================
        AliRICH *RICH = new AliRICHv3("RICH", "normal RICH");

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
	    // With hole in front of RICH
	    TRD->SetRICHhole();
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

    if (iSTART)
    {
        //=================== START parameters ============================
        AliSTART *START = new AliSTARTv1("START", "START Detector");
    }

    if (iEMCAL)
    {
        //=================== EMCAL parameters ============================
        AliEMCAL *EMCAL = new AliEMCALv1("EMCAL", "EMCALArch1a");
    }

     if (iCRT)
    {
        //=================== CRT parameters ============================
        AliCRT *CRT = new AliCRTv0("CRT", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== CRT parameters ============================
        AliVZERO *VZERO = new AliVZEROv2("VZERO", "normal VZERO");
    }
 
             
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}



AliGenerator* GeneratorFactory(PprRun_t run) {
    Int_t isw = 3;
    if (rad == kNoGluonRadiation) isw = 0;
    

    switch (run) {
    case test50:
	comment = comment.Append(":HIJINGparam test 50 particles");
	AliGenHIJINGpara *gener = new AliGenHIJINGpara(50);
	gener->SetMomentumRange(0, 999999.);
	gener->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
	Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
	gener->SetThetaRange(thmin,thmax);
	break;
    case kParam_8000:
	comment = comment.Append(":HIJINGparam N=8000");
	AliGenHIJINGpara *gener = new AliGenHIJINGpara(86030);
	gener->SetMomentumRange(0, 999999.);
	gener->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
	Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
	gener->SetThetaRange(thmin,thmax);
	break;
    case kParam_4000:
	comment = comment.Append("HIJINGparam N=4000");
	AliGenHIJINGpara *gener = new AliGenHIJINGpara(43015);
	gener->SetMomentumRange(0, 999999.);
	gener->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
	Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
	gener->SetThetaRange(thmin,thmax);
	break;
    case kParam_2000:
	comment = comment.Append("HIJINGparam N=2000");
	AliGenHIJINGpara *gener = new AliGenHIJINGpara(21507);
	gener->SetMomentumRange(0, 999999.);
	gener->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
	Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
	gener->SetThetaRange(thmin,thmax);
	break;
//
//  Hijing Central
//
    case kHijing_cent1:
	comment = comment.Append("HIJING cent1");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	break;
    case kHijing_cent2:
	comment = comment.Append("HIJING cent2");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 2.);
	break;
//
// Hijing Peripheral 
//
    case kHijing_per1:
	comment = comment.Append("HIJING per1");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(5., 8.6);
	break;
    case kHijing_per2:
	comment = comment.Append("HIJING per2");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(8.6, 11.2);
	break;
    case kHijing_per3:
	comment = comment.Append("HIJING per3");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(11.2, 13.2);
	break;
    case kHijing_per4:
	comment = comment.Append("HIJING per4");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(13.2, 15.);
	break;
    case kHijing_per5:
	comment = comment.Append("HIJING per5");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(15., 100.);
	break;
//
//  Jet-Jet
//
    case kHijing_jj25:
	comment = comment.Append("HIJING Jet 25 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(1);
	gener->SetPtJet(25.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.3,0.3);
	gener->SetJetPhiRange(75., 165.);   
	break;

    case kHijing_jj50:
	comment = comment.Append("HIJING Jet 50 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(1);
	gener->SetPtJet(50.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.3,0.3);
	gener->SetJetPhiRange(75., 165.);   
	break;

    case kHijing_jj75:
	comment = comment.Append("HIJING Jet 75 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(1);
	gener->SetPtJet(75.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.3,0.3);
	gener->SetJetPhiRange(75., 165.);   
	break;

    case kHijing_jj100:
	comment = comment.Append("HIJING Jet 100 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(1);
	gener->SetPtJet(100.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.3,0.3);
	gener->SetJetPhiRange(75., 165.);   
	break;

    case kHijing_jj200:
	comment = comment.Append("HIJING Jet 200 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(1);
	gener->SetPtJet(200.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.3,0.3);
	gener->SetJetPhiRange(75., 165.);   
	break;
//
// Gamma-Jet
//
    case kHijing_gj25:
	comment = comment.Append("HIJING Gamma 25 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(2);
	gener->SetPtJet(25.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.12, 0.12);
        gener->SetJetPhiRange(220., 320.);
	break;

    case kHijing_gj50:
	comment = comment.Append("HIJING Gamma 50 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(2);
	gener->SetPtJet(50.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.12, 0.12);
        gener->SetJetPhiRange(220., 320.);
	break;

    case kHijing_gj75:
	comment = comment.Append("HIJING Gamma 75 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(2);
	gener->SetPtJet(75.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.12, 0.12);
        gener->SetJetPhiRange(220., 320.);
	break;

    case kHijing_gj100:
	comment = comment.Append("HIJING Gamma 100 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(2);
	gener->SetPtJet(100.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.12, 0.12);
        gener->SetJetPhiRange(220., 320.);
	break;

    case kHijing_gj200:
	comment = comment.Append("HIJING Gamma 200 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(2);
	gener->SetPtJet(200.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.12, 0.12);
        gener->SetJetPhiRange(220., 320.);
	break;
    }
    return gener;
}

AliGenHijing* HijingStandard()
{
    AliGenHijing *gener = new AliGenHijing(-1);
// centre of mass energy 
    gener->SetEnergyCMS(5500.);
// reference frame
    gener->SetReferenceFrame("CMS");
// projectile
     gener->SetProjectile("A", 208, 82);
     gener->SetTarget    ("A", 208, 82);
// tell hijing to keep the full parent child chain
     gener->KeepFullEvent();
// enable jet quenching
     gener->SetJetQuenching(4);
// enable shadowing
     gener->SetShadowing(1);
// neutral pion and heavy particle decays switched off
     gener->SetDecaysOff(1);
// Don't track spectators
     gener->SetSpectators(0);
// kinematic selection
     gener->SetSelectAll(0);
     return gener;
}



