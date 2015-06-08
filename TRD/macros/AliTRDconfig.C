#include <TPDGCode.h>

void Config()
{
  Int_t iField = 0;

  // libraries required by geant321
  gSystem->Load("libgeant321");

  new TGeant3TGeo("C++ Interface to Geant3");

  // Create the output file
  TFile *rootfile = new TFile("TRD_test.root","recreate");
  rootfile->SetCompressionLevel(2);

  // Define the monte carlo
  TGeant3 *geant3 = (TGeant3*) gMC;

  AliRunLoader* rl=0x0;
  cout << "AliTRDconfig.C: Creating Run Loader ..." <<endl;
  rl = AliRunLoader::Open("TRD_test.root"
                         ,AliConfig::GetDefaultEventFolderName()
                         ,"recreate");
  if (rl == 0x0) {
    gAlice->Fatal("AliTRDconfig.C","Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(3);
  gAlice->SetRunLoader(rl);

  // Set external decayer
  TVirtualMCDecayer *decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);

  // ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
  geant3->SetTRIG(1); // Number of events to be processed 
  geant3->SetSWIT(4,10);
  geant3->SetDEBU(0,0,1);
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
  geant3->SetHADR(1); // Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
  geant3->SetLOSS(2);
  geant3->SetMULS(1);
  geant3->SetRAYL(1);
  geant3->SetAUTO(1); // Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
  geant3->SetABAN(0); // Restore 3.16 behaviour for abandoned tracks
  geant3->SetOPTI(2); // Select optimisation level for GEANT geometry searches (0,1,2)
  geant3->SetERAN(5.e-7);

  Float_t cut    = 1.e-3; // 1MeV cut by default
  Float_t tofmax = 1.e10;
  //             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
  geant3->SetCUTS(cut,cut, cut, cut, cut, cut,  cut,  cut, cut,  cut, tofmax);

  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV

  AliGenCocktail *gener = new AliGenCocktail();

  AliGenBox *genEl = new AliGenBox(100);
  genEl->SetOrigin(0,0,0);        // Vertex position
  genEl->SetSigma(0,0,0);         // Sigma in (X,Y,Z) (cm) on IP position
  genEl->SetPart(kElectron);             // Only electrons 

  AliGenBox *genPi = new AliGenBox(100);
  genPi->SetOrigin(0,0,0);        // Vertex position
  genPi->SetSigma(0,0,0);         // Sigma in (X,Y,Z) (cm) on IP position
  genPi->SetPart(kPiMinus);           // Only pions 

  gener->AddGenerator(genEl,"Electrons",1);
  gener->AddGenerator(genPi,"Pions"    ,1);

  if (iField) {

    // With magnetic field on
    AliGenerator *gg = gener->FirstGenerator()->Generator();
    gg->SetMomentumRange(3.00,3.01);
    gg->SetPhiRange(76.0,92.0);
    gg->SetThetaRange(83.0,97.0);
    gg = gener->NextGenerator()->Generator();
    gg->SetMomentumRange(0.560,0.561);
    gg->SetPhiRange(62.0,78.0);
    gg->SetThetaRange(83.0,97.0);

    gener->Init();

    // Specify maximum magnetic field in Tesla (neg. ==> default field)
    // 0.4 T
    gAlice->SetField(-999,2,2.0);    

  }
  else {

    // With magnetic field off
    AliGenerator *gg = gener->FirstGenerator()->Generator();
    gg->SetMomentumRange(3.00,3.01);
    gg->SetPhiRange(82.0,98.0);
    gg->SetThetaRange(83.0,97.0);
    gg = gener->NextGenerator()->Generator();
    gg->SetMomentumRange(0.560,0.561);
    gg->SetPhiRange(82.0,98.0);
    gg->SetThetaRange(83.0,97.0);

    gener->Init();

    // Specify maximum magnetic field in Tesla (neg. ==> default field)
    // No field
    gAlice->SetField(0);    

  }

  Int_t iMAG   = 1;
  Int_t iITS   = 0;
  Int_t iTPC   = 0;
  Int_t iTRD   = 1;
  Int_t iABSO  = 1;
  Int_t iDIPO  = 1;
  Int_t iHALL  = 1;
  Int_t iFRAME = 1;
  Int_t iSHIL  = 1;
  Int_t iPIPE  = 1;

  rl->CdGAFile();

  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY","Alice envelop");

  if (iMAG) {
    //=================== MAG parameters ============================
    // --- Start with Magnet since detector layouts may be depending ---
    // --- on the selected Magnet dimensions ---
    AliMAG *MAG  = new AliMAG("MAG","Magnet");
  }

  if (iABSO) { 
    //=================== ABSO parameters ============================
    AliABSO *ABSO  = new AliABSOv0("ABSO","Muon Absorber");
  }

  if (iDIPO) {
    //=================== DIPO parameters ============================
    AliDIPO *DIPO  = new AliDIPOv2("DIPO","Dipole version 2");
  }

  if (iHALL) {
    //=================== HALL parameters ============================
    AliHALL *HALL  = new AliHALL("HALL","Alice Hall");
  }

  if (iFRAME) {
    //=================== FRAME parameters ============================
    AliFRAMEv2 *FRAME  = new AliFRAMEv2("FRAME","Space Frame");
    FRAME->SetHoles(0);
  }

  if (iSHIL) {
    //=================== SHIL parameters ============================
    AliSHIL *SHIL  = new AliSHILv0("SHIL","Shielding");
  }

  if (iPIPE) {
    //=================== PIPE parameters ============================
    AliPIPE *PIPE  = new AliPIPEv0("PIPE","Beam Pipe");
  }

  if (iITS) {
    //=================== ITS parameters ============================
    //

    AliITSvPPRasymmFMD *ITS = new AliITSvPPRasymmFMD("ITS","ITS PPR");
    ITS->SetMinorVersion(2);
    ITS->SetReadDet(kTRUE);
    ITS->SetThicknessDet1(200.);
    ITS->SetThicknessDet2(200.);
    ITS->SetThicknessChip1(200.);
    ITS->SetThicknessChip2(200.);
    ITS->SetRails(0);
    ITS->SetCoolingFluid(1);
    ITS->SetEUCLID(0);

  }

  if (iTPC) {
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

    AliTPC *TPC  = new AliTPCv2("TPC","Default");

  }

  if (iTRD) {
    //=================== TRD parameters ============================
  
    AliTRD *TRD  = new AliTRDv1("TRD","TRD slow simulator");
  
    // Set to detailed display
    TRD->SetDisplayType(1);

    // Draw TR photons
    TRD->SetDrawTR(1);

  }
        
}
