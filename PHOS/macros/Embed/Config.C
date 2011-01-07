#include <TPDGCode.h>

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
static PprGeo_t sgeo = kNoHoles;
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
  UInt_t at = (UInt_t) gSystem->Now() ;
  UInt_t seed = ((gSystem->GetPid()*111)%at)*137 ;
//  gRandom->SetSeed(seed);
  gRandom->SetSeed(12345);
  printf("MySeed: %d\n",seed) ;
  cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<endl;

  
  
  // libraries required by fluka21

  Bool_t isFluka = kFALSE;
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
  rl->SetNumberOfEventsPerFile(1000);
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
  gMC->SetProcess("DRAY",0); //AZ 1);
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

  ((AliMC*)gMC)->SetTransPar("./galice.cuts") ;
  //
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV
 
    AliGenBox *gener = new AliGenBox(5);
    gener->SetMomentumRange(0.5, 5.);
    gener->SetPhiRange(260., 280.);
    gener->SetThetaRange(82.,98.);
    gener->SetPart(kGamma);

    gener->SetOrigin(0, 0, 0);  //vertex position
    gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position
    gener->Init() ;
 
  // 
  // Activate this line if you want the vertex smearing to happen
  // track by track
  //
  //  gener->SetVertexSmear(kPerEvent) ;
 


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
        comment = comment.Append(" | Holes for PHOS/RICH");
                                                                                
    } else {
        comment = comment.Append(" | No holes for PHOS/RICH");
    }
                                                                                
    printf("\n \n Comment: %s \n \n", comment.Data());
                                                                                
                                                                                
// Field (L3 0.4 T)
    //Zero magnetic field
    AliMagF* field = new AliMagF("Maps","Maps", 0., 0., AliMagF::k5kGUniform);
    //    AliMagF* field = new AliMagF("Maps","Maps", 2, -1., -1., 10., smag);
    TGeoGlobalMagField::Instance()->SetField(field);

    rl->CdGAFile();
 
  Int_t   iABSO  = 0; 
  Int_t   iCRT   = 0; 
  Int_t   iDIPO  = 0; 
  Int_t   iFMD   = 0; 
  Int_t   iFRAME = 0; 
  Int_t   iHALL  = 0; 
  Int_t   iITS   = 0; 
  Int_t   iMAG   = 0; 
  Int_t   iMUON  = 0; 
  Int_t   iPHOS  = 1; 
  Int_t   iPIPE  = 0; 
  Int_t   iPMD   = 0; 
  Int_t   iRICH  = 0; 
  Int_t   iSHIL  = 0; 
  Int_t   iSTART = 0; 
  Int_t   iTOF   = 0; 
  Int_t   iTPC   = 0;
  Int_t   iTRD   = 0; 
  Int_t   iZDC   = 0; 
  Int_t   iEMCAL = 0; 
  Int_t   iVZERO = 0;
 
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
        if (sgeo == kHoles) {
            FRAME->SetHoles(1);
        } else {
            FRAME->SetHoles(0);
        }
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
        ITS->SetReadDet(kTRUE);   // don't touch this parameter if you're not an ITS developer
    //    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
        ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
        ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
        ITS->SetThicknessChip1(200.);  // chip thickness on layer 1 must be in the range [150,300]
        ITS->SetThicknessChip2(200.);  // chip thickness on layer 2 must be in the range [150,300]
        ITS->SetRails(0);            // 1 --> rails in ; 0 --> rails out
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
        //============================ TPC parameters ================================
//        AliTPC *TPC = new AliTPCv0("TPC", "Default");
        AliTPC *TPC = new AliTPCv2("TPC", "Default");
    }
                                                                                
                                                                                
    if (iTOF) {
        //=================== TOF parameters ============================
        AliTOF *TOF = new AliTOFv4T0("TOF", "normal TOF");
    }
                                                                                
                                                                                
    if (iRICH)
    {
        //=================== RICH parameters ===========================
        AliRICH *RICH = new AliRICHv1("RICH", "normal RICH");
                                                                                
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
                                                                                
        // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 90% Xe + 10% CO2)
        TRD->SetGasMix(1);
        if (sgeo == kHoles) {
            // With hole in front of PHOS
            TRD->SetPHOShole();
            // With hole in front of RICH
            TRD->SetRICHhole();
        }
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
       AliPHOS *PHOS = new AliPHOSv1("PHOS", "IHEP");
//         AliPHOS *PHOS = new AliPHOSv1("PHOS", "noCPV");
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
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETEV1");
    }
                                                                                
     if (iCRT)
    {
        //=================== CRT parameters ============================
        AliCRT *CRT = new AliCRTv0("CRT", "normal ACORDE");
    }
                                                                                
     if (iVZERO)
    {
        //=================== CRT parameters ============================
        AliVZERO *VZERO = new AliVZEROv3("VZERO", "normal VZERO");
    }
                                                                                
                                                                                
}
                                                                                
Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

