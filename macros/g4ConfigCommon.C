// $Id$
//
// AliRoot Configuration for running aliroot with Monte Carlo.
// Called from g4Config.C
//
// By I. Hrivnacova, IPN Orsay

enum PprGeo_t
  {
    kHoles, kNoHoles
  };
static PprGeo_t geo = kHoles;

Float_t EtaToTheta(Float_t arg);
void    LoadPythia();
static Int_t    eventsPerRun = 50;

void ConfigCommon(Bool_t setRootGeometry = kTRUE)
{
  cout << "Running ConfigCommon.C ... " << endl;

  // Load Pythia libraries
  LoadPythia();

  // Create the output file
  AliRunLoader* rl = 0;
  rl = AliRunLoader::Open("galice.root",
			   AliConfig::GetDefaultEventFolderName(),
			   "recreate");
  if (!rl) {
    gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(3);
  gAlice->SetRunLoader(rl);
  
  // Set Root geometry file
  //
  if ( setRootGeometry ) {
    gAlice->SetRootGeometry();
    gAlice->SetGeometryFromFile("geometry.root");
  }

  // Set the trigger configuration
  gAlice->SetTriggerDescriptor("Pb-Pb");
  cout<<"Trigger configuration is set to  Pb-Pb"<<endl;

    // Set Random Number seed
    gRandom->SetSeed(123456); // Set 0 to use the currecnt time
    AliLog::Message(AliLog::kInfo, Form("Seed for random number generation = %d",gRandom->GetSeed()), "Config.C", "Config.C", "Config()","Config.C", __LINE__);
    int     nParticles = 100;
    if (gSystem->Getenv("CONFIG_NPARTICLES"))
    {
        nParticles = atoi(gSystem->Getenv("CONFIG_NPARTICLES"));
    }

    AliGenCocktail *gener = new AliGenCocktail();
    gener->SetPhiRange(0, 360);
    // Set pseudorapidity range from -8 to 8.
    Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
    Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
    gener->SetThetaRange(thmin,thmax);
    gener->SetOrigin(0, 0, 0);  //vertex position
    gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position

    AliGenHIJINGpara *hijingparam = new AliGenHIJINGpara(nParticles);
    hijingparam->SetMomentumRange(0.2, 999);
    gener->AddGenerator(hijingparam,"HIJING PARAM",1);
    gener->Init();

    // 
    // Activate this line if you want the vertex smearing to happen
    // track by track
    //
    //gener->SetVertexSmear(perTrack); 

  // ============================= 
  // Magnetic field
  // ============================= 

    // Field (L3 0.4 T)
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG));

    Int_t   iABSO  =  1;
    Int_t   iDIPO  =  1;
    Int_t   iFMD   =  1;
    Int_t   iFRAME =  1;
    Int_t   iHALL  =  1;
    Int_t   iITS   =  1;
    Int_t   iMAG   =  1;
    Int_t   iMUON  =  1;
    Int_t   iPHOS  =  1;
    Int_t   iPIPE  =  1;
    Int_t   iPMD   =  1;
    Int_t   iHMPID =  1;
    Int_t   iSHIL  =  1;
    Int_t   iT0    =  1;
    Int_t   iTOF   =  1;
    Int_t   iTPC   =  1;
    Int_t   iTRD   =  1;
    Int_t   iZDC   =  1;
    Int_t   iEMCAL =  1;
    Int_t   iACORDE = 0;
    Int_t   iVZERO =  1;
/*
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
    Int_t   iHMPID =  0;
    Int_t   iSHIL  =  0;
    Int_t   iT0    =  0;
    Int_t   iTOF   =  0;
    Int_t   iTPC   =  1;
    Int_t   iTRD   =  0;
    Int_t   iZDC   =  0;
    Int_t   iEMCAL =  0;
    Int_t   iACORDE = 0;
    Int_t   iVZERO =  0;
*/   
// Exluded for problems

    rl->CdGAFile();
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

        AliDIPO *DIPO = new AliDIPOv3("DIPO", "Dipole version 3");
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
        if (geo == kHoles) {
	  FRAME->SetHoles(1);
	} else {
	  FRAME->SetHoles(0);
	}
    }

    if (iSHIL)
    {
        //=================== SHIL parameters ============================

        AliSHIL *SHIL = new AliSHILv3("SHIL", "Shielding Version 3");
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
        //============================ TPC parameters ===================
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

        AliZDC *ZDC = new AliZDCv3("ZDC", "normal ZDC");
    }

    if (iTRD)
    {
        //=================== TRD parameters ============================

        AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
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

    if (iT0)
    {
        //=================== T0 parameters ============================
        AliT0 *T0 = new AliT0v1("T0", "T0 Detector");
    }

    if (iEMCAL)
    {
        //=================== EMCAL parameters ============================
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETE");
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

     AliLog::Message(AliLog::kInfo, "End of Config", "Config.C", "Config.C", "Config()"," Config.C", __LINE__);

     cout << "Running ConfigCommon.C finished ... " << endl;

}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

void LoadPythia()
{
    // Load Pythia related libraries
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
}
