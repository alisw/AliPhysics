// Function converting pseudorapidity
// interval to polar angle interval. It is used to set 
// the limits in the generator
Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

// Set Random Number seed using the current time
TDatime dat;
static UInt_t sseed = dat.Get();

void Config()
{
  gRandom->SetSeed(sseed);
  cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<endl; 
    
  // Load GEANT 3 library. It has to be in LD_LIBRARY_PATH
  gSystem->Load("libgeant321");

  // Instantiation of the particle transport package. gMC is set internaly
  new TGeant3TGeo("C++ Interface to Geant3");

  // Create run loader and set some properties
  AliRunLoader* rl =  AliRunLoader::Open("galice.root",
					 AliConfig::GetDefaultEventFolderName(),
					 "recreate");
  if (!rl) Fatal("Config.C","Can not instatiate the Run Loader");
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(3);

  // Register the run loader in gAlice
  gAlice->SetRunLoader(rl);

  // Set external decayer
  TVirtualMCDecayer *decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll); // kAll means no specific decay is forced
  decayer->Init();

  // Register the external decayer in the transport package
  gMC->SetExternalDecayer(decayer);

  // STEERING parameters FOR ALICE SIMULATION
  // Specify event type to be transported through the ALICE setup
  // All positions are in cm, angles in degrees, and P and E in GeV
  // For the details see the GEANT 3 manual

  // Switch on/off the physics processes (global)
  // Please consult the file data/galice.cuts for detector
  // specific settings, i.e. DRAY
  gMC->SetProcess("DCAY",1); // Particle decay
  gMC->SetProcess("PAIR",1); // Pair production
  gMC->SetProcess("COMP",1); // Compton scattering
  gMC->SetProcess("PHOT",1); // Photo effect
  gMC->SetProcess("PFIS",0); // Photo fission
  gMC->SetProcess("DRAY",0); // Delta rays
  gMC->SetProcess("ANNI",1); // Positron annihilation
  gMC->SetProcess("BREM",1); // Bremstrahlung
  gMC->SetProcess("MUNU",1); // Muon nuclear interactions
  gMC->SetProcess("CKOV",1); // Cerenkov production
  gMC->SetProcess("HADR",1); // Hadronic interactions
  gMC->SetProcess("LOSS",2); // Energy loss (2=complete fluct.)
  gMC->SetProcess("MULS",1); // Multiple scattering
  gMC->SetProcess("RAYL",1); // Rayleigh scattering
  
  // Set the transport package cuts
  Float_t cut = 1.e-3;        // 1MeV cut by default
  Float_t tofmax = 1.e10;

  gMC->SetCut("CUTGAM", cut); // Cut for gammas
  gMC->SetCut("CUTELE", cut); // Cut for electrons
  gMC->SetCut("CUTNEU", cut); // Cut for neutral hadrons
  gMC->SetCut("CUTHAD", cut); // Cut for charged hadrons
  gMC->SetCut("CUTMUO", cut); // Cut for muons
  gMC->SetCut("BCUTE",  cut); // Cut for electron brems.
  gMC->SetCut("BCUTM",  cut); // Cut for muon brems.
  gMC->SetCut("DCUTE",  cut); // Cut for electron delta-rays
  gMC->SetCut("DCUTM",  cut); // Cut for muon delta-rays
  gMC->SetCut("PPCUTM", cut); // Cut for e+e- pairs by muons
  gMC->SetCut("TOFMAX", tofmax); // Time of flight cut
  
  // Set up the particle generation

  // AliGenCocktail permits to combine several different generators
  AliGenCocktail *gener = new AliGenCocktail();

  // The phi range is always inside 0-360
  gener->SetPhiRange(0, 360);

  // Set pseudorapidity range from -8 to 8.
  Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
  Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
  gener->SetThetaRange(thmin,thmax);

  gener->SetOrigin(0, 0, 0);  // vertex position
  gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetCutVertexZ(1.);     // Truncate at 1 sigma
  gener->SetVertexSmear(kPerEvent); 

  // First cocktail component: 100 ``background'' particles 
  AliGenHIJINGpara *hijingparam = new AliGenHIJINGpara(100);
  hijingparam->SetMomentumRange(0.2, 999);
  gener->AddGenerator(hijingparam,"HIJING PARAM",1);

  // Second cocktail component: one gamma in PHOS direction
  AliGenBox *genbox = new AliGenBox(1);
  genbox->SetMomentumRange(10,11.);
  genbox->SetPhiRange(270.5,270.7);
  genbox->SetThetaRange(90.5,90.7);
  genbox->SetPart(kGamma);
  gener->AddGenerator(genbox,"GENBOX GAMMA for PHOS",1);

  gener->Init(); // Initialization of the coctail generator

  // Field (the last parameter is 1 => L3 0.4 T)
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));

  // Make sure the current ROOT directory is in galice.root 
  rl->CdGAFile();

  // Build the setup and set some detector parameters

  // ALICE BODY parameters. BODY is always present
  AliBODY *BODY = new AliBODY("BODY", "ALICE envelop");

  // Start with Magnet since detector layouts may be depending
  // on the selected Magnet dimensions
  AliMAG *MAG = new AliMAG("MAG", "Magnet");

  AliABSO *ABSO = new AliABSOv0("ABSO", "Muon Absorber");       // Absorber

  AliDIPO *DIPO = new AliDIPOv2("DIPO", "Dipole version 2");    // Dipole magnet

  AliHALL *HALL = new AliHALL("HALL", "ALICE Hall");            // Hall
    
  AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");   // Space frame

  AliSHIL *SHIL = new AliSHILv2("SHIL", "Shielding Version 2"); // Shielding

  AliPIPE *PIPE = new AliPIPEv0("PIPE", "Beam Pipe");           // Beam pipe
    
  // ITS parameters
  AliITSvPPRasymmFMD *ITS  = new AliITSvPPRasymmFMD("ITS",
	       "ITS PPR detailed version with asymmetric services");
  ITS->SetMinorVersion(2);      // don't change it if you're not an ITS developer
  ITS->SetReadDet(kFALSE);      // don't change it if you're not an ITS developer
  ITS->SetThicknessDet1(200.);  // detector thickness on layer 1:[100,300] mkm
  ITS->SetThicknessDet2(200.);  // detector thickness on layer 2:[100,300] mkm
  ITS->SetThicknessChip1(150.); // chip thickness on layer 1: [150,300] mkm
  ITS->SetThicknessChip2(150.); // chip thickness on layer 2: [150,300]
  ITS->SetRails(0);	        // 1 --> rails in ; 0 --> rails out
  ITS->SetCoolingFluid(1);      // 1 --> water ; 0 --> freon
  ITS->SetEUCLID(0);            // no output for the EUCLID CAD system 

  
  AliTPC *TPC = new AliTPCv2("TPC", "Default");                 // TPC

  AliTOF *TOF = new AliTOFv5T0("TOF", "normal TOF");            // TOF

  AliHMPID *HMPID = new AliHMPIDv1("HMPID", "normal HMPID");         // HMPID

  AliZDC *ZDC = new AliZDCv2("ZDC", "normal ZDC");              // ZDC

  AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");      // TRD

  AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");              // FMD

  AliMUON *MUON = new AliMUONv1("MUON", "default");             // MUON

  AliPHOS *PHOS = new AliPHOSv1("PHOS", "IHEP");                // PHOS

  AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");              // PMD

  AliT0 *T0 = new AliT0v1("T0", "T0 Detector");  // T0

  // EMCAL
  AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETEV1");

  AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");    // VZERO
}
