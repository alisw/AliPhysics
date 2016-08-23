// The common part contains the description of geometry,
// magnetic field and generator. The initialization of
// the random generator seed is also here.
// The specific transport packages is created in the main
// configuration file. The specific decayer is also set in
// the main configuration file since it is needed by the transport.

// Triggers
enum PprTrigConf_t { kDefaultPPTrig, kDefaultPbPbTrig };
const char * pprTrigConfName[] = { "p-p","Pb-Pb" };

// Defaults
Int_t    sseed = 12345; //Set 0 to use the current time
AliMagF::BMap_t smag = AliMagF::k5kG; // Magnetic field
PprTrigConf_t strig = kDefaultPPTrig; // default PP trigger configuration

// Conversion function that might be used in the generator to set eta ranges
Float_t EtaToTheta(Float_t arg) {
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

void commonConfig(){
  // The logical groups are separated using comment lines
  
  //=======================================================================
  // Set Random Number seed
  gRandom->SetSeed(sseed);
  cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<endl; 

  //=======================================================================
  // Run loader
  AliRunLoader* rl=0x0;

  cout<<"Config.C: Creating Run Loader ..."<<endl;
  rl = AliRunLoader::Open("galice.root",
			  AliConfig::GetDefaultEventFolderName(),
			  "recreate");
  if (rl == 0x0) {
    gAlice->Fatal("commonConfig.C","Can not instatiate the Run Loader");
    return;
  }
  
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(100);
  gAlice->SetRunLoader(rl);
  rl->CdGAFile();

  //=======================================================================
  // Field (L3 0.5 T)
  AliMagF* field = new AliMagF("Maps","Maps",-1., -1., smag);
  TGeoGlobalMagField::Instance()->SetField(field);
  if (smag == AliMagF::k2kG) {
    cout << "L3 field 0.2 T" << endl;
  } else if (smag == AliMagF::k5kG) {
    cout << "L3 field 0.5 T" << endl;
  }

  // Set the trigger configuration
  AliSimulation::Instance()->SetTriggerConfig(pprTrigConfName[strig]);
  cout<<"Trigger configuration is set to  "<<pprTrigConfName[strig]<<endl;
  
  //=======================================================================
  // Debug and log level
  //    AliLog::SetGlobalDebugLevel(0);
  //    AliLog::SetGlobalLogLevel(AliLog::kError);
  
  //=======================================================================
  // Generator Configuration
  AliGenExtFile *gener = new AliGenExtFile(-1);
  AliGenReaderTreeK * reader = new AliGenReaderTreeK();
  
  reader->SetFileName("galice.root");
  reader->AddDir("../gen");
  gener->SetReader(reader);
  
  gener->SetOrigin(0, 0, 0);    // vertex position
  gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetCutVertexZ(1.);     // Truncate at 1 sigma
  gener->SetVertexSmear(kPerEvent); 
  gener->SetTrackingFlag(1);
  gener->Init();
  
  //=======================================================================
  // Switches for the detector simulation
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
  Int_t   iHMPID  = 1;
  Int_t   iSHIL   = 1;
  Int_t   iT0     = 1;
  Int_t   iTOF    = 1;
  Int_t   iTPC    = 1;
  Int_t   iTRD    = 1;
  Int_t   iZDC    = 1;
  Int_t   iEMCAL  = 1;
  Int_t   iVZERO  = 1;
  Int_t   iACORDE = 0;
  
  // Alice BODY parameters
  AliBODY *BODY = new AliBODY("BODY", "Alice envelop");
  
  if (iMAG) {
    // MAG parameters 
    // Start with Magnet since detector layouts may be depending
    // on the selected Magnet dimensions
    AliMAG *MAG = new AliMAG("MAG", "Magnet");
  }
  
  if (iABSO) {
    // ABSO parameters 
    AliABSO *ABSO = new AliABSOv3("ABSO", "Muon Absorber");
  }
  
  if (iDIPO) {
    // DIPO parameters 
    AliDIPO *DIPO = new AliDIPOv3("DIPO", "Dipole version 3");
  }
  
  if (iHALL) {
    // HALL parameters 
    AliHALL *HALL = new AliHALLv3("HALL", "Alice Hall");
  }
  
  if (iFRAME) {
    // FRAME parameters 
    AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
    FRAME->SetHoles(1);
  }
  
  if (iSHIL) {
    // SHIL parameters 
    AliSHIL *SHIL = new AliSHILv3("SHIL", "Shielding Version 3");
  }
  
  if (iPIPE) {
    // PIPE parameters 
    AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
  }
  
  if (iITS) {
    // ITS parameters 
    AliITS *ITS  = new AliITSv11("ITS","ITS v11");
  }
  
  if (iTPC) {
    // TPC parameters
    AliTPC *TPC = new AliTPCv2("TPC", "Default");
  }
  
  if (iTOF) {
    // TOF parameters 
    AliTOF *TOF = new AliTOFv6T0("TOF", "normal TOF");
  }
    
  if (iHMPID) {
    // HMPID parameters
    AliHMPID *HMPID = new AliHMPIDv3("HMPID", "normal HMPID");
  }
  
  if (iZDC) {
    // ZDC parameters 
    AliZDC *ZDC = new AliZDCv4("ZDC", "normal ZDC");
  }
  
  if (iTRD) {
    // TRD parameters 
    AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
  }
  
  if (iFMD) {
    // FMD parameters 
    AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
  }
  
  if (iMUON) {
    // MUON parameters
    // New MUONv1 version (geometry defined via builders)
    AliMUON *MUON = new AliMUONv1("MUON", "default");
  }
  
  if (iPHOS) {
    // PHOS parameters
    AliPHOS *PHOS = new AliPHOSv1("PHOS", "Run1");
  }
  
  if (iPMD) {
    // PMD parameters 
    AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");
  }
  
  if (iT0) {
    // T0 parameters 
    AliT0 *T0 = new AliT0v1("T0", "T0 Detector");
  }
  
  if (iEMCAL) {
    // EMCAL parameters 
    AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETEV1");
  }
  
  if (iACORDE) {
    // ACORDE parameters 
    AliACORDE *ACORDE = new AliACORDEv1("ACORDE", "normal ACORDE");
  }
  
  if (iVZERO) {
    // VZERO parameters 
    AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
  }
}
