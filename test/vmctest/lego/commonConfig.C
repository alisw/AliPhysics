// $Id$
//
// AliRoot Configuration for running aliroot with Monte Carlo.
// commonConfig() includes the common setting for all MCs
// which has to be called before MC is instantiated.
// Called from MC specific configs (g3Config.C, g4Config.C).
//
// Extracted from G3 specific Config.C 
// by I. Hrivnacova, IPN Orsay

enum PprTrigConf_t {
    kDefaultPPTrig, kDefaultPbPbTrig
};

const char * pprTrigConfName[] = {
    "p-p","Pb-Pb"
};

enum ConfigVersion_t {
    kConfigV0,  // default configuration  
    kConfigV1   // configuration for LHC 2010 production
};

// Options 
static AliMagF::BMap_t smag = AliMagF::k5kG;
static PprTrigConf_t strig = kDefaultPPTrig; // default PP trigger configuration
static TString comment;

// Functions
void  LoadPythia();

void commonConfig(const TString& det, 
                  ConfigVersion_t configVersion = kConfigV0)
{
  cout << "Running commonConfig.C ... " << endl;

    // Set Random Number seed
  gRandom->SetSeed(123456); // Set 0 to use the currecnt time
  AliLog::Message(AliLog::kInfo, Form("Seed for random number generation = %d",gRandom->GetSeed()), "Config.C", "Config.C", "Config()","Config.C", __LINE__);


  //=======================================================================
  // Load Pythia libraries
  //=======================================================================

  LoadPythia();

  //=======================================================================
  // ALICE steering object (AliRunLoader)
  //=======================================================================

  AliRunLoader* rl 
    = AliRunLoader::Open("galice.root",
			  AliConfig::GetDefaultEventFolderName(),
			  "recreate");
  if ( ! rl ) {
    gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(3);
  gAlice->SetRunLoader(rl);
  
  //======================================================================
  // Trigger configuration
  //=======================================================================

  AliSimulation::Instance()->SetTriggerConfig(pprTrigConfName[strig]);
  cout << "Trigger configuration is set to  " << pprTrigConfName[strig] << endl;

  // ============================= 
  // Magnetic field
  // ============================= 

  // Field (L3 0.5 T)
  AliMagF* field = new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG);
  TGeoGlobalMagField::Instance()->SetField(field);

  printf("\n \n Comment: %s \n \n", comment.Data());

  // ============================= 
  // Modules
  // ============================= 

  rl->CdGAFile();

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
  Int_t   iTPC   =  0;
  Int_t   iTRD   =  0;
  Int_t   iZDC   =  0;
  Int_t   iEMCAL =  0;
  Int_t   iACORDE = 0;
  Int_t   iVZERO =  0;

  if ( det == "ABSO" )  iABSO  =  1;
  if ( det == "DIPO" )  iDIPO  =  1;
  if ( det == "FMD" )   iFMD   =  1;
  if ( det == "FRAME" )  iFRAME =  1;
  if ( det == "HALL" )  iHALL  =  1;
  if ( det == "ITS" )   iITS   =  1;
  if ( det == "MAG" )   iMAG   =  1;
  if ( det == "MUON" )  iMUON  =  1;
  if ( det == "PHOS" )  iPHOS  =  1;
  if ( det == "PIPE" )  iPIPE  =  1;
  if ( det == "PMD" )   iPMD   =  1;
  if ( det == "HMPID" )  iHMPID =  1;
  if ( det == "SHIL" )  iSHIL  =  1;
  if ( det == "T0" )    iT0    =  1;
  if ( det == "TOF" )   iTOF   =  1;
  if ( det == "TPC" )   iTPC   =  1;
  if ( det == "TRD" )   iTRD   =  1;
  if ( det == "ZDC" )   iZDC   =  1;
  if ( det == "EMCAL" )  iEMCAL =  1;
  if ( det == "ACORDE" ) iACORDE = 1;
  if ( det == "VZERO" )  iVZERO =  1;

  if ( det == "ALL" ) {
    iABSO  =  1;
    iDIPO  =  1;
    iFMD   =  1;
    iFRAME =  1;
    iHALL  =  1;
    iITS   =  1;
    iMAG   =  1;
    iMUON  =  1;
    iPHOS  =  1;
    iPIPE  =  1;
    iPMD   =  1;
    iHMPID =  1;
    iSHIL  =  1;
    iT0    =  1;
    iTOF   =  1;
    iTPC   =  1;
    iTRD   =  1;
    iZDC   =  1;
    iEMCAL =  1;
    iACORDE = 1;
    iVZERO =  1;
  }

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
      if ( configVersion == kConfigV1 ) {
        AliTRDgeometry *geoTRD = TRD->GetGeometry();
        // Partial geometry: modules at 0,1,7,8,9,16,17
        // starting at 3h in positive direction
        geoTRD->SetSMstatus(2,0);
        geoTRD->SetSMstatus(3,0);
        geoTRD->SetSMstatus(4,0);
        geoTRD->SetSMstatus(5,0);
        geoTRD->SetSMstatus(6,0);
        geoTRD->SetSMstatus(11,0);
        geoTRD->SetSMstatus(12,0);
        geoTRD->SetSMstatus(13,0);
        geoTRD->SetSMstatus(14,0);
        geoTRD->SetSMstatus(15,0);
        geoTRD->SetSMstatus(16,0);
      }
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
     if ( configVersion == kConfigV0 ) 
       AliPHOS *PHOS = new AliPHOSv1("PHOS", "IHEP");
     else if ( configVersion == kConfigV1 )  
       AliPHOS *PHOS = new AliPHOSv1("PHOS", "noCPV_Modules123"); 
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
    if ( configVersion == kConfigV0 ) 
      AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETEV1");
    else if ( configVersion == kConfigV1 )  
      AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_FIRSTYEARV1");
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

  AliLog::Message(AliLog::kInfo, "End of Config", "Config.C", "Config.C", "Config()"," Config.C", __LINE__);

  cout << "Running commonConfig.C finished ... " << endl;
}

void LoadPythia()
{
  // Load Pythia related libraries
  gSystem->Load("liblhapdf.so");      // Parton density functions
  gSystem->Load("libEGPythia6.so");   // TGenerator interface
  gSystem->Load("libpythia6.so");     // Pythia
  gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
}
