enum gentype_t {hijing, hijingParam, gun, box, pythia,                          
                param1, param2, param3, param4,                  
		cocktail, fluka, halo, ntuple, scan, doublescan};               

//gentype_t gentype=hijingParam;           
//gentype_t gentype=hijing;
  gentype_t gentype=
  //box;
  cocktail;
  //hijingParam;
  //hijing;

void Config()
{
    // 7-DEC-2000 09:00
    // Switch on Transition Radiation simulation. 6/12/00 18:00
    // iZDC=1  7/12/00 09:00
    // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00

    // Set Random Number seed
    // gRandom->SetSeed(12345);

    new     AliGeant3("C++ Interface to Geant3");

    if (!gSystem->Getenv("CONFIG_FILE"))
    {
        TFile  *rootfile = new TFile("galice.root", "recreate");
	//TFile  *rootfile = new TFile("galice.root", "recreate");
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
    if (gSystem->Getenv("CONFIG_NPARTICLES"))
    {
        int     nParticles = atoi(gSystem->Getenv("CONFIG_NPARTICLES"));
    } else
    {
        int     nParticles = 50;
    }
//-----------------------------------------------------------
switch(gentype)                                                                
{                                                                              
   case hijing:      
//-----------------------------------------HIJING------------
AliGenHijing *gener = new AliGenHijing(-1); 
    gener->SetEnergyCMS(5500);                                                 
    gener->SetReferenceFrame("CMS     ");                                      
    gener->SetProjectile("A       ", 208, 82);                                 
    gener->SetTarget    ("A       ", 208, 82);                                 
    gener->SetImpactParameterRange(0, 3.);                                     
    gener->SetEvaluate(1);                                                     
    gener->KeepFullEvent();                                                    
    gener->SetJetQuenching(1);                                                 
    gener->SetShadowing(1);   
    gener->SetDecaysOff(1);                                                      
    gener->SetTrigger(0);                                                      
    gener->SetSelectAll(0);                                                    
    gener->SetMomentumRange( 0.05 , 1.3 );                                            
    gener->SetPhiRange( 0.0 , 1.0 );                                               
//    gener->SetPhiRange(-180.,180.) 
    gener->SetThetaRange(45.0,135.);                                              
    gener->SetFlavor(0);   // 0 - all, 4 - charm and beaty                                                       
    gener->SetOrigin(0., 0.0 ,0);                                              
    gener->SetSigma(0,0,5.3);                                                  
    gener->SetVertexSmear(kPerEvent);                                          
//    gener->SetTrackingFlag(0);                                                 
    gener->Init();
    break;
//-----------------------------------------HijingPARA----------
    case hijingParam:
    
    AliGenHIJINGpara *gener = new AliGenHIJINGpara(20);

    gener->SetMomentumRange(0.1, .7);
    gener->SetThetaRange(45.,135.);
    gener->SetPhiRange(0, 360);
    //  gener->SetThetaRange(0.28,179.72);
    //    gener->SetThetaRange(45., 135.);
    gener->SetOrigin(0, 0, 0);  //vertex position
    gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position
    gener->Init();
    break;
//------------------------------------------ Box --------------
 case box:  
   AliGenBox *gener = new AliGenBox(4000);  //ntracks=100
     gener->SetMomentumRange(0.2,1.0);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(45, 135. );
     gener->SetOrigin(0,0,0);   
     //vertex position
     gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(321);             // K+
     //gener->SetPart(211);           // Pi+
     break;

    //------------------------------------------------------------
 case cocktail:
   AliGenCocktail *gener = new AliGenCocktail();
   int Nmax=50;
   AliGenBox *gener1 = new AliGenBox(2*Nmax);  //ntracks=100
     gener1->SetMomentumRange(.050,1.00);
     //gener1->SetMomentumRange(.050,0.80);
     gener1->SetPhiRange(0,360);
     gener1->SetThetaRange(45, 135. );
     gener1->SetOrigin(0,0,0);   
     gener1->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener1->SetPart(321);             // K+

   AliGenBox *gener2 = new AliGenBox(2*Nmax);  //ntracks=100
     gener2->SetMomentumRange(.050,1.00);
     //gener1->SetMomentumRange(.050,0.80);
     gener2->SetPhiRange(0,360);
     gener2->SetThetaRange(45, 135. );
     gener2->SetOrigin(0,0,0);   
     gener2->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener2->SetPart(2212);             // Protons

   AliGenBox *gener3 = new AliGenBox(Nmax);  //ntracks=100
     gener3->SetMomentumRange(.050,1.00);
     //gener1->SetMomentumRange(.050,0.80);
     gener3->SetPhiRange(0,360);
     gener3->SetThetaRange(45, 135. );
     gener3->SetOrigin(0,0,0);   
     gener3->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener3->SetPart(11);             // e+


   AliGenBox *gener4 = new AliGenBox(Nmax);  //ntracks=100
     gener4->SetMomentumRange(.050,1.00);
     //gener1->SetMomentumRange(.050,0.80);
     gener4->SetPhiRange(0,360);
     gener4->SetThetaRange(45, 135. );
     gener4->SetOrigin(0,0,0);   
     gener4->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener4->SetPart(211);             // pi+

     gener->SetMomentumRange(.050,1.00);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(45., 135. );
     gener->SetOrigin(0,0,0);   
     gener->SetSigma(0,0,0); 
     gener->AddGenerator(gener1,"box1",0.5);
     gener->AddGenerator(gener2,"box2",0.5);
     gener->AddGenerator(gener3,"box3",0.5);
     gener->AddGenerator(gener4,"box4",0.5);
     gener->Init();
   break;

    default:
    break;
}
    // 
    // Activate this line if you want the vertex smearing to happen
    // track by track
    //
    //gener->SetVertexSmear(perTrack); 

    gAlice->SetField(-999, 2);  //Specify maximum magnetic field in Tesla (neg. ==> default field)

    Int_t   iABSO = 0;
    Int_t   iCASTOR = 0;
    Int_t   iDIPO = 0;
    Int_t   iFMD = 0;
    Int_t   iFRAME = 0;
    Int_t   iHALL = 0;
    Int_t   iITS = 1;
    Int_t   iMAG = 0;
    Int_t   iMUON = 0;
    Int_t   iPHOS = 0;
    Int_t   iPIPE = 1;
    Int_t   iPMD = 0;
    Int_t   iRICH = 0;
    Int_t   iSHIL = 0;
    Int_t   iSTART = 0;
    Int_t   iTOF = 0;
    Int_t   iTPC = 1;
    Int_t   iTRD = 0;
    Int_t   iZDC = 0;

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

        AliFRAME *FRAME = new AliFRAMEv2("FRAME", "Space Frame");

    }

    if (iSHIL)
    {
        //=================== SHIL parameters ============================

        AliSHIL *SHIL = new AliSHILv0("SHIL", "Shielding");
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
    ITS->SetMinorVersion(2);                                         // don't touch this parameter if you're not an ITS developer
    ITS->SetReadDet(kFALSE);                                         // don't touch this parameter if you're not an ITS developer
    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
    ITS->SetThicknessDet1(300.);   // detector thickness on layer 1 must be in the range [100,300]
    ITS->SetThicknessDet2(300.);   // detector thickness on layer 2 must be in the range [100,300]
    ITS->SetThicknessChip1(300.);  // chip thickness on layer 1 must be in the range [150,300]
    ITS->SetThicknessChip2(300.);  // chip thickness on layer 2 must be in the range [150,300]
    ITS->SetRails(1);	           // 1 --> rails in ; 0 --> rails out
    ITS->SetCoolingFluid(1);       // 1 --> water ; 0 --> freon
    //
    //AliITSvPPRsymm *ITS  = new AliITSvPPRsymm("ITS","New ITS PPR detailed version with symmetric services");
    //ITS->SetMinorVersion(2);                                       // don't touch this parameter if you're not an ITS developer
    //ITS->SetReadDet(kFALSE);                                       // don't touch this parameter if you're not an ITS developer
    //ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRsymm2.det"); // don't touch this parameter if you're not an ITS developer
    //ITS->SetThicknessDet1(300.);   // detector thickness on layer 1 must be in the range [100,300]
    //ITS->SetThicknessDet2(300.);   // detector thickness on layer 2 must be in the range [100,300]
    //ITS->SetThicknessChip1(300.);  // chip thickness on layer 1 must be in the range [150,300]
    //ITS->SetThicknessChip2(300.);  // chip thickness on layer 2 must be in the range [150,300]
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

    if (iTOF)
    {
        //=================== TOF parameters ============================
        AliTOF *TOF = new AliTOFv2("TOF", "normal TOF");
    }

    if (iRICH)
    {
        //=================== RICH parameters ===========================
        AliRICH *RICH = new AliRICHv1("RICH", "normal RICH");

    }


    if (iZDC)
    {
        //=================== ZDC parameters ============================

        AliZDC *ZDC = new AliZDCv1("ZDC", "normal ZDC");
    }

    if (iCASTOR)
    {
        //=================== CASTOR parameters ============================

        AliCASTOR *CASTOR = new AliCASTORv1("CASTOR", "normal CASTOR");
    }

    if (iTRD)
    {
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
        //=================== FMD parameters ============================

        AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
    }

    if (iMUON)
    {
        //=================== MUON parameters ===========================

        AliMUON *MUON = new AliMUONv1("MUON", "default");
    }
    //=================== PHOS parameters ===========================

    if (iPHOS)
    {
        AliPHOS *PHOS = new AliPHOSv1("PHOS", "GPS2");
    }


    if (iPMD)
    {
        //=================== PMD parameters ============================

        AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");

        PMD->SetPAR(1., 1., 0.8, 0.02);
        PMD->SetIN(6., 18., -580., 27., 27.);
        PMD->SetGEO(0.0, 0.2, 4.);
        PMD->SetPadSize(0.8, 1.0, 1.0, 1.5);

    }

    if (iSTART)
    {
        //=================== START parameters ============================
        AliSTART *START = new AliSTARTv1("START", "START Detector");
    }


}
