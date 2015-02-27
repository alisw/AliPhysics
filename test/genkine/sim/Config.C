// Configuration of simulation

enum PprRad_t
{
    kGluonRadiation, kNoGluonRadiation
};

enum PprTrigConf_t
{
    kDefaultPPTrig, kDefaultPbPbTrig
};

const char * pprTrigConfName[] = {
    "p-p","Pb-Pb"
};


// This part for configuration    

static PprRad_t srad = kGluonRadiation;
static AliMagF::BMap_t smag = AliMagF::k5kG;
static Int_t    sseed = 12345; //Set 0 to use the current time

static PprTrigConf_t strig = kDefaultPPTrig; // default PbPb trigger configuration
// Comment line 
static TString  comment;

// Functions
Float_t EtaToTheta(Float_t arg);
AliGenerator* GeneratorFactory();

void Config()
{
    // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
    // Theta range given through pseudorapidity limits 22/6/2001

    // Set Random Number seed
    gRandom->SetSeed(sseed);
    cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<endl; 
  


   // libraries required by geant321 and Pythia6: loaded in sim.C

    new     TGeant3TGeo("C++ Interface to Geant3");

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
    rl->SetNumberOfEventsPerFile(100);
    gAlice->SetRunLoader(rl);

  // Set the trigger configuration
    AliSimulation::Instance()->SetTriggerConfig(pprTrigConfName[strig]);
    cout<<"Trigger configuration is set to  "<<pprTrigConfName[strig]<<endl;


    //
    // Set External decayer
    AliDecayer *decayer = new AliDecayerPythia();
    decayer->SetForceDecay(kAll);
    decayer->Init();

    //forbid some decays
    AliPythia * py= AliPythia::Instance();
    py->SetMDME(737,1,0); //forbid D*+->D+ + pi0
    py->SetMDME(738,1,0);//forbid D*+->D+ + gamma

    for(Int_t d=747; d<=762; d++){ 
      py->SetMDME(d,1,0);
    }

    for(Int_t d=764; d<=807; d++){ 
      py->SetMDME(d,1,0);
    }

    TVirtualMC * vmc = TVirtualMC::GetMC();
   
    vmc->SetExternalDecayer(decayer);
    //
    //
    //=======================================================================
    //
    //=======================================================================
    // ************* STEERING parameters FOR ALICE SIMULATION **************
    // --- Specify event type to be tracked through the ALICE setup
    // --- All positions are in cm, angles in degrees, and P and E in GeV

    vmc->SetProcess("DCAY",1);
    vmc->SetProcess("PAIR",1);
    vmc->SetProcess("COMP",1);
    vmc->SetProcess("PHOT",1);
    vmc->SetProcess("PFIS",0);
    vmc->SetProcess("DRAY",0);
    vmc->SetProcess("ANNI",1);
    vmc->SetProcess("BREM",1);
    vmc->SetProcess("MUNU",1);
    vmc->SetProcess("CKOV",1);
    vmc->SetProcess("HADR",1);
    vmc->SetProcess("LOSS",2);
    vmc->SetProcess("MULS",1);
    vmc->SetProcess("RAYL",1);

    Float_t cut = 1.e-3;        // 1MeV cut by default
    Float_t tofmax = 1.e10;

    vmc->SetCut("CUTGAM", cut);
    vmc->SetCut("CUTELE", cut);
    vmc->SetCut("CUTNEU", cut);
    vmc->SetCut("CUTHAD", cut);
    vmc->SetCut("CUTMUO", cut);
    vmc->SetCut("BCUTE",  cut); 
    vmc->SetCut("BCUTM",  cut); 
    vmc->SetCut("DCUTE",  cut); 
    vmc->SetCut("DCUTM",  cut); 
    vmc->SetCut("PPCUTM", cut);
    vmc->SetCut("TOFMAX", tofmax); 

    // Debug and log level
    //    AliLog::SetGlobalDebugLevel(0);
    //    AliLog::SetGlobalLogLevel(AliLog::kError);

    // Generator Configuration
    AliGenerator* gener = GeneratorFactory();
    gener->SetOrigin(0, 0, 0);    // vertex position
    gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position
    gener->SetCutVertexZ(1.);     // Truncate at 1 sigma
    gener->SetVertexSmear(kPerEvent); 
    gener->SetTrackingFlag(1);
    gener->Init();
    
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


    printf("\n \n Comment: %s \n \n", comment.Data());
    
    
// Field (L3 0.4 T)
    AliMagF* field = new AliMagF("Maps","Maps",-1., -1., smag);
    TGeoGlobalMagField::Instance()->SetField(field);

    rl->CdGAFile();
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
    Int_t   iHMPID   = 1;
    Int_t   iSHIL   = 1;
    Int_t   iT0  = 1;
    Int_t   iTOF    = 1;
    Int_t   iTPC    = 1;
    Int_t   iTRD    = 1;
    Int_t   iZDC    = 1;
    Int_t   iEMCAL  = 1;
    Int_t   iVZERO  = 1;
    Int_t   iACORDE    = 0;

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
	FRAME->SetHoles(1);
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

	AliITS *ITS  = new AliITSv11("ITS","ITS v11");
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

        AliZDC *ZDC = new AliZDCv4("ZDC", "normal ZDC");
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
        AliPHOS *PHOS = new AliPHOSv1("PHOS", "Run1");
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
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETEV1");
    }

     if (iACORDE)
    {
        //=================== ACORDE parameters ============================
        AliACORDE *ACORDE = new AliACORDEv1("ACORDE", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== VZERO parameters ============================
        AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
    }
 
             
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}



AliGenerator* GeneratorFactory() {

  AliGenExtFile *gener = new AliGenExtFile(-1);
  AliGenReaderTreeK * reader = new AliGenReaderTreeK();

  reader->SetFileName("galice.root");
  reader->AddDir("../gen");
  gener->SetReader(reader);
     
  return gener;

 
}

