static Int_t    eventsPerRun = 100;
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
static PprGeo_t sgeo = kHoles;
static PprRad_t srad = kGluonRadiation;
static PprMag_t smag = k5kG;
                                                                                
// Comment line
static TString  comment;
                                                                                
// Functions
Float_t EtaToTheta(Float_t arg);


void Config()
{
  cout << "==> Config.C..." << endl;
  
  // Set Random Number seed
  gRandom->SetSeed(12345);
  cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<endl;

  
  
  // libraries required by fluka21

  Bool_t isFluka = kTRUE;
  if (isFluka) {
    gSystem->Load("libGeom");
    cout << "\t* Loading TFluka..." << endl;  
    gSystem->Load("libTFluka");    
    
    cout << "\t* Instantiating TFluka..." << endl;
    new  TFluka("C++ Interface to Fluka", 0/*verbositylevel*/);
  }
  else {
    cout << "\t* Loading Geant3..." << endl;  
    gSystem->Load("libgeant321");
    
    cout << "\t* Instantiating Geant3TGeo..." << endl;
    new     TGeant3TGeo("C++ Interface to Geant3");
  }
  
  AliRunLoader* rl=0x0;
                                                                                
  cout<<"Config.C: Creating Run Loader ..."<<endl;
  rl = AliRunLoader::Open("galice.root",
			    AliConfig::GetDefaultEventFolderName(),
			    "recreate");
  if (rl == 0x0)
    {
      gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
      return;
    }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(3);
  gAlice->SetRunLoader(rl);
                                                                                
  //
  // Set External decayer
  AliDecayer *decayer = new AliDecayerPythia();
                                                                               
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);

  //
  //
  //
  // Physics process control

  gMC->SetProcess("DCAY",1);
  gMC->SetProcess("PAIR",1);
  gMC->SetProcess("COMP",1);
  gMC->SetProcess("PHOT",1);
  gMC->SetProcess("PFIS",0);
  gMC->SetProcess("DRAY",1);
  gMC->SetProcess("ANNI",1);
  gMC->SetProcess("BREM",1);
  gMC->SetProcess("MUNU",1);
  gMC->SetProcess("CKOV",1); 
  gMC->SetProcess("HADR",1);
  gMC->SetProcess("LOSS",2);
  gMC->SetProcess("MULS",1);
  gMC->SetProcess("RAYL",1);
                                                                                
  Float_t cut = 1.e-3;        // 1MeV cut by default
  Float_t tofmax = 1.e10;
                                                                                
  gMC->SetCut("CUTGAM", cut);
  gMC->SetCut("CUTELE", cut);
  gMC->SetCut("CUTNEU", cut);
  gMC->SetCut("CUTHAD", cut);
  gMC->SetCut("CUTMUO", cut);
  gMC->SetCut("BCUTE",  cut);
  gMC->SetCut("BCUTM",  cut);
  gMC->SetCut("DCUTE",  cut);
  gMC->SetCut("DCUTM",  cut);
  gMC->SetCut("PPCUTM", cut);
  gMC->SetCut("TOFMAX", tofmax);

  //
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV
  if (gSystem->Getenv("CONFIG_NPARTICLES"))
      int     nParticles = atoi(gSystem->Getenv("CONFIG_NPARTICLES"));
  else
      int     nParticles = 1000;
  
  cout << "\t* Creating and configuring generator for " << nParticles 
       << " particles..." << endl;
  
  AliGenHIJINGpara *gener = new AliGenHIJINGpara(nParticles);
  
  gener->SetMomentumRange(0., 999);
  gener->SetPhiRange(0, 360);
  // Set pseudorapidity range from -3 to 3.
  Float_t thmin = EtaToTheta( 3.);   // theta min. <---> eta max
  Float_t thmax = EtaToTheta(-3.);   // theta max. <---> eta min 
  gener->SetThetaRange(thmin,thmax);
  gener->SetOrigin(0, 0, 0);  //vertex position
  gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position
  gener->Init();
  // 
  // Activate this line if you want the vertex smearing to happen
  // track by track
  //


     if (smag == k2kG) {
        comment = comment.Append(" | L3 field 0.2 T");
    } else if (smag == k4kG) {
        comment = comment.Append(" | L3 field 0.4 T");
    } else if (smag == k5kG) {
        comment = comment.Append(" | L3 field 0.5 T");
    }
                                                                                
                                                                                
    if (srad == kGluonRadiation)
    {
        comment = comment.Append(" | Gluon Radiation On");
                                                                                
    } else {
        comment = comment.Append(" | Gluon Radiation Off");
    }
                                                                                
    if (sgeo == kHoles)
    {
        comment = comment.Append(" | Holes for PHOS/HMPID");
                                                                                
    } else {
        comment = comment.Append(" | No holes for PHOS/HMPID");
    }
                                                                                
    printf("\n \n Comment: %s \n \n", comment.Data());
                                                                                
                                                                                
// Field (L3 0.4 T)
    AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., smag);
    field->SetL3ConstField(0); //Using const. field in the barrel
    rl->CdGAFile();
    gAlice->SetField(field);
 
  Int_t   iABSO  = 1; 
  Int_t   iCRT   = 1; 
  Int_t   iDIPO  = 1; 
  Int_t   iFMD   = 1; 
  Int_t   iFRAME = 1; 
  Int_t   iHALL  = 1; 
  Int_t   iITS   = 1; 
  Int_t   iMAG   = 1; 
  Int_t   iMUON  = 1; 
  Int_t   iPHOS  = 1; 
  Int_t   iPIPE  = 1; 
  Int_t   iPMD   = 1; 
  Int_t   iHMPID  = 1; 
  Int_t   iSHIL  = 1; 
  Int_t   iSTART = 1; 
  Int_t   iTOF   = 1; 
  Int_t   iTPC   = 1;
  Int_t   iTRD   = 1; 
  Int_t   iZDC   = 1; 
  Int_t   iEMCAL = 0; 
  Int_t   iVZERO = 1;
 
  cout << "\t* Creating the detectors ..." << endl;
  //=================== Alice BODY parameters =============================
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
	AliITSvPPRasymmFMD *ITS  = new AliITSvPPRasymmFMD("ITS","New ITS PPR detailed version with asymmetric services");
	ITS->SetMinorVersion(2);  // don't touch this parameter if you're not an ITS developer
	ITS->SetReadDet(kFALSE);	  // don't touch this parameter if you're not an ITS developer
    //    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
	ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
	ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
	ITS->SetThicknessChip1(150.);  // chip thickness on layer 1 must be in the range [150,300]
	ITS->SetThicknessChip2(150.);  // chip thickness on layer 2 must be in the range [150,300]
	ITS->SetRails(0);	     // 1 --> rails in ; 0 --> rails out
	ITS->SetCoolingFluid(1);   // 1 --> water ; 0 --> freon

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
      //============================ TPC parameters =====================
        AliTPC *TPC = new AliTPCv2("TPC", "Default");
    }


    if (iTOF) {
        //=================== TOF parameters ============================
	AliTOF *TOF = new AliTOFv5T0("TOF", "normal TOF");
	// Partial geometry: modules at 2,3,4,6,7,11,12,14,15,16
	// starting at 6h in positive direction
	//	Int_t TOFSectors[18]={-1,-1,0,0,0,-1,0,0,-1,-1,-1,0,0,-1,0,0,0,0};
	// Partial geometry: modules at 1,2,6,7,9,10,11,12,15,16,17
	// (ALICE numbering convention)
       	Int_t TOFSectors[18]={-1,0,0,-1,-1,-1,0,0,-1,0,0,0,0,-1,-1,0,0,0};
	TOF->SetTOFSectors(TOFSectors);
    }


    if (iHMPID)
    {
        //=================== HMPID parameters ===========================
        AliHMPID *HMPID = new AliHMPIDv1("HMPID", "normal HMPID");

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
        AliTRDgeometry *geoTRD = TRD->GetGeometry();
	// Partial geometry: modules at 2,3,4,6,11,12,14,15
	// starting at 6h in positive direction
	geoTRD->SetSMstatus(0,0);
        geoTRD->SetSMstatus(1,0);
        geoTRD->SetSMstatus(5,0);
        geoTRD->SetSMstatus(7,0);
        geoTRD->SetSMstatus(8,0);
        geoTRD->SetSMstatus(9,0);
        geoTRD->SetSMstatus(10,0);
        geoTRD->SetSMstatus(13,0);
        geoTRD->SetSMstatus(16,0);
        geoTRD->SetSMstatus(17,0);
    }

    if (iFMD)
    {
        //=================== FMD parameters ============================
	AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
   }

    if (iMUON)
    {
        //=================== MUON parameters ===========================
        // New MUONv1 version (geometry defined via builders)
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
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "SHISH_77_TRD1_2X2_FINAL_110DEG");
    }

     if (iCRT)
    {
        //=================== CRT parameters ============================
        AliCRT *CRT = new AliCRTv0("CRT", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== CRT parameters ============================
        AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
    }                                                                                
}
                                                                                
Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

