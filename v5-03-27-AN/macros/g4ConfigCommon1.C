// $Id: g4ConfigCommon.C 30849 2009-02-01 11:42:22Z fca $
//
// AliRoot Configuration for running aliroot with Monte Carlo.
// ConfigCommon1() includes the common setting for all MCs
// which has to be called before MC is instantiated.
// Called from g4Config.C
//
// By I. Hrivnacova, IPN Orsay

// Options 
static AliMagF::BMap_t smag = AliMagF::k5kG;
static TString comment;

// Functions
void  LoadPythia();

void ConfigCommon1(Bool_t setRootGeometry = kTRUE)
{
  cout << "Running ConfigCommon1.C ... " << endl;

  //=======================================================================
  // Load Pythia libraries
  //=======================================================================

  LoadPythia();

  //=======================================================================
  // ALICE steering object (AliRunLoader)
  //=======================================================================

  // Set Root geometry file
  if ( setRootGeometry ) {
    gAlice->SetRootGeometry();
    gAlice->SetGeometryFromFile("geometry.root");
  }

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
  
  //=======================================================================
  // Trigger configuration
  //=======================================================================

  AliSimulation::Instance()->SetTriggerConfig("Pb-Pb");
  cout<<"Trigger configuration is set to  Pb-Pb"<<endl;

  // ============================= 
  // Magnetic field
  // ============================= 

  // Field (L3 0.4 T)
  if (smag == AliMagF::k2kG) {
      comment = comment.Append(" | L3 field 0.2 T");
  } 
  else if (smag == AliMagF::k5kG) {
      comment = comment.Append(" | L3 field 0.5 T");
  }
  // OK
  AliMagF* field = new AliMagF("Maps","Maps", -1., -1., smag);
  TGeoGlobalMagField::Instance()->SetField(field);

  printf("\n \n Comment: %s \n \n", comment.Data());

  // ============================= 
  // Modules
  // ============================= 

  rl->CdGAFile();

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

      AliITS *ITS  = new AliITSv11("ITS","ITS v11");
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

  cout << "Running ConfigCommon1.C finished ... " << endl;
}

void LoadPythia()
{
  // Load Pythia related libraries
  gSystem->Load("liblhapdf.so");      // Parton density functions
  gSystem->Load("libEGPythia6.so");   // TGenerator interface
  gSystem->Load("libpythia6.so");     // Pythia
  gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
}
