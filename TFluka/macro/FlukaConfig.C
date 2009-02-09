static Int_t    eventsPerRun = 100;
enum PprGeo_t
{
    kHoles, kNoHoles
};
                                                                                
enum PprRad_t
{
    kGluonRadiation, kNoGluonRadiation
};
                                                                                
// This part for configuration
static PprGeo_t sgeo = kHoles;
static PprRad_t srad = kGluonRadiation;
static AliMagF::BMap_t smag = AliMagF::k5kG;
                                                                                
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

  gSystem->Load("liblhapdf.so");      // Parton density functions
  gSystem->Load("libEGPythia6.so");   // TGenerator interface
  gSystem->Load("libpythia6.so");     // Pythia
  gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
  
  

  Bool_t isFluka = kTRUE;
  if (isFluka) {
    gSystem->Load("libGeom");
    cout << "\t* Loading TFluka..." << endl;  
    gSystem->Load("libfluka");    
    
    cout << "\t* Instantiating TFluka..." << endl;
    new  TFluka("C++ Interface to Fluka", 0/*verbositylevel*/);
  }
  else {
    cout << "\t* Loading Geant3..." << endl;  
    gSystem->Load("libgeant321");
    
    cout << "\t* Instantiating Geant3TGeo..." << endl;
    new     TGeant3TGeo("C++ Interface to Geant3");
  }
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()){
      AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
      AliCDBManager::Instance()->SetRun(0);
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


    if (smag == AliMagF::k2kG) {
        comment = comment.Append(" | L3 field 0.2 T");
    } else if (smag == AliMagF::k5kG) {
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
    AliMagF *field = new AliMagF("Maps","Maps", 2, 1., 1., 10., smag);
    field->SetL3ConstField(0); //Using const. field in the barrel
    TGeoGlobalMagField::Instance()->SetField(field);
    rl->CdGAFile();
 
  Int_t   iABSO    = 0; 
  Int_t   iACORDE  = 0; 
  Int_t   iDIPO    = 0; 
  Int_t   iFMD     = 0; 
  Int_t   iFRAME   = 0; 
  Int_t   iHALL    = 0; 
  Int_t   iITS     = 0; 
  Int_t   iMAG     = 0; 
  Int_t   iMUON    = 0; 
  Int_t   iPHOS    = 0; 
  Int_t   iPIPE    = 0; 
  Int_t   iPMD     = 0; 
  Int_t   iHMPID   = 0; 
  Int_t   iSHIL    = 0; 
  Int_t   iT0      = 0; 
  Int_t   iTOF     = 0; 
  Int_t   iTPC     = 0;
  Int_t   iTRD     = 0; 
  Int_t   iZDC     = 0; 
  Int_t   iEMCAL   = 0; 
  Int_t   iVZERO   = 0;
 
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
        AliABSO *ABSO = new AliABSOv3("ABSO", "Muon Absorber");
    }

    if (iDIPO)
    {
        //=================== DIPO parameters ============================

        AliDIPO *DIPO = new AliDIPOv3("DIPO", "Dipole version 2");
    }

    if (iHALL)
    {
        //=================== HALL parameters ============================

        AliHALL *HALL = new AliHALLv3("HALL", "Alice Hall");
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

        AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
    }
 
    if (iITS)
    {
        //=================== ITS parameters ============================

	AliITS *ITS  = new AliITSv11Hybrid("ITS","ITS v11Hybrid");
    }


    if (iTPC)
    {
      //============================ TPC parameters =====================
        AliTPC *TPC = new AliTPCv2("TPC", "Default");
    }


    if (iTOF) {
        //=================== TOF parameters ============================
	AliTOF *TOF = new AliTOFv6T0("TOF", "normal TOF");
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
        //Set simulation parameters different from the default ones.
        AliPHOSSimParam* simEmc = AliPHOSSimParam::GetInstance() ;
  
        // APD noise of warm (+20C) PHOS:
        // a2 = a1*(Y1/Y2)*(M1/M2), where a1 = 0.012 is APD noise at -25C,
        // Y1 = 4.3 photo-electrons/MeV, Y2 = 1.7 p.e/MeV - light yields at -25C and +20C,
        // M1 = 50, M2 = 50 - APD gain factors chosen for t1 = -25C and t2 = +20C,
        // Y = MeanLightYield*APDEfficiency.

        Float_t apdNoise = 0.012*2.5; 
        simEmc->SetAPDNoise(apdNoise);

        //Raw Light Yield at +20C
        simEmc->SetMeanLightYield(18800);

        //ADC channel width at +18C.
        simEmc->SetADCchannelW(0.0125);
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
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "SHISH_77_TRD1_2X2_FINAL_110DEG");
    }

     if (iACORDE)
    {
        //=================== ACORDE parameters ============================
        AliACORDE *ACORDE = new AliACORDEv1("ACORDE", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== ACORDE parameters ============================
        AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
    }                                                                                
}
                                                                                
Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

