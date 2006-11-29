//
// Configuration for the Physics Data Challenge 2006 modified to generate
// e,mu,pi,K,p for PID studies of the TRD.
// Per event, AliGenBox is used to generate 100 particles per species
// (50 positive and 50 negative ones). Barrel detectors only.
// s.masciocchi@gsi.de
//

// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config_PDC06.C++")

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenParam.h"
#include "EVGEN/AliGenMUONlib.h"
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "PYTHIA6/AliGenPythia.h"
#include "STEER/AliMagFMaps.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"
#include "STRUCT/AliABSOv0.h"
#include "STRUCT/AliDIPOv2.h"
#include "STRUCT/AliHALL.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv2.h"
#include "STRUCT/AliPIPEv0.h"
#include "ITS/AliITSgeom.h"
#include "ITS/AliITSvPPRasymmFMD.h"
#include "TPC/AliTPCv2.h"
#include "TOF/AliTOFv5T0.h"
#include "HMPID/AliHMPIDv1.h"
#include "ZDC/AliZDCv2.h"
#include "TRD/AliTRDv1.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "CRT/AliCRTv0.h"
#include "VZERO/AliVZEROv6.h"
#endif


enum PDC06Proc_t 
{
//--- Heavy Flavour Production ---
  kCharmPbPb5500,  kCharmpPb8800,  kCharmpp14000,  kCharmpp14000wmi,
  kD0PbPb5500,     kD0pPb8800,     kD0pp14000,
  kDPlusPbPb5500,  kDPluspPb8800,  kDPluspp14000,
  kBeautyPbPb5500, kBeautypPb8800, kBeautypp14000, kBeautypp14000wmi, 
// -- Pythia Mb
  kPyMbNoHvq, kPyOmegaPlus, kPyOmegaMinus, kRunMax
}; 

const char * pprRunName[] = {
  "kCharmPbPb5500",  "kCharmpPb8800",  "kCharmpp14000",  "kCharmpp14000wmi",
  "kD0PbPb5500",     "kD0pPb8800",     "kD0pp14000",
  "kDPlusPbPb5500",  "kDPluspPb8800",  "kDPluspp14000",
  "kBeautyPbPb5500", "kBeautypPb8800", "kBeautypp14000", "kBeautypp14000wmi", 
  "kPyMbNoHvq", "kPyOmegaPlus", "kPyOmegaMinus"
};


//--- Decay Mode ---
enum DecayHvFl_t 
{
  kNature,  kHadr, kSemiEl, kSemiMu
};
//--- Rapidity Cut ---
enum YCut_t
{
  kFull, kBarrel, kMuonArm
};
//--- Magnetic Field ---
enum Mag_t
{
    k2kG, k4kG, k5kG
};
//--- Functions ---
AliGenPythia *PythiaHVQ(PDC06Proc_t proc);
AliGenerator *MbCocktail();
AliGenerator *PyMbTriggered(Int_t pdg);
void ProcessEnvironmentVars();
Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

// This part for configuration
static DecayHvFl_t   decHvFl   = kNature; 
static YCut_t        ycut      = kBarrel;
static Mag_t         mag       = k5kG; 
//========================//
// Set Random Number seed //
//========================//
TDatime dt;
static UInt_t seed    = dt.Get();

// nEvts = -1  : you get 1 QQbar pair and all the fragmentation and 
//               decay chain
// nEvts = N>0 : you get N charm / beauty Hadrons 
Int_t nEvts = -1; 
// stars = kTRUE : all heavy resonances and their decay stored
//       = kFALSE: only final heavy hadrons and their decays stored
Bool_t stars = kTRUE;

// To be used only with kCharmppMNRwmi and kBeautyppMNRwmi
// To get a "reasonable" agreement with MNR results, events have to be 
// generated with the minimum ptHard set to 2.76 GeV.
// To get a "perfect" agreement with MNR results, events have to be 
// generated in four ptHard bins with the following relative 
// normalizations:
//  CHARM
// 2.76-3 GeV: 25%
//    3-4 GeV: 40%
//    4-8 GeV: 29%
//     >8 GeV:  6%
//  BEAUTY
// 2.76-4 GeV:  5% 
//    4-6 GeV: 31%
//    6-8 GeV: 28%
//     >8 GeV: 36%
Float_t ptHardMin =  2.76;
Float_t ptHardMax = -1.;


// Comment line
static TString comment;

void Config()
{
 

  // Get settings from environment variables
  ProcessEnvironmentVars();

  gRandom->SetSeed(seed);
  cerr<<"Seed for random number generation= "<<seed<<endl; 

  // libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("libgeant321");
#endif

  new TGeant3TGeo("C++ Interface to Geant3");

  //=======================================================================
  //  Create the output file

   
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
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV


    gMC->SetProcess("DCAY",1);
    gMC->SetProcess("PAIR",1);
    gMC->SetProcess("COMP",1);
    gMC->SetProcess("PHOT",1);
    gMC->SetProcess("PFIS",0);
    gMC->SetProcess("DRAY",0);
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




  // Set External decayer //
  //======================//
  TVirtualMCDecayer* decayer = new AliDecayerPythia();
  // DECAYS
  //
  switch(decHvFl) {
  case kNature:
    decayer->SetForceDecay(kAll);
    break;
  case kHadr:
    decayer->SetForceDecay(kHadronicD);
    break;
  case kSemiEl:
    decayer->SetForceDecay(kSemiElectronic);
    break;
  case kSemiMu:
    decayer->SetForceDecay(kSemiMuonic);
    break;
  }
  decayer->Init();
  gMC->SetExternalDecayer(decayer);

  //=========================//
  // Generator Configuration //
  //=========================//
  
  AliGenCocktail *gener  = new AliGenCocktail();
  Double_t momen         = 0.6;
  Int_t Npart            = 50;
  // Set pseudorapidity range from -1 to 1.
  Float_t thmin          = EtaToTheta(1.);   // theta min. <---> eta max
  Float_t thmax          = EtaToTheta(-1.);  // theta max. <---> eta min
  
  // electron generator
  AliGenBox *eminus= new AliGenBox(Npart);
  eminus->SetPart(11);
  eminus->SetMomentumRange(momen,momen);
  eminus->SetPhiRange(0,360);
  eminus->SetThetaRange(thmin,thmax);
  eminus->SetOrigin(0., 0., 0.);    // primary vertex position
  eminus->SetSigma(0.,0.,0.);       // sigma (x,y,z) in cm, on position

  // positron generator
  AliGenBox *eplus= new AliGenBox(Npart);
  eplus->SetPart(-11);
  eplus->SetMomentumRange(momen,momen);
  eplus->SetPhiRange(0,360);
  eplus->SetThetaRange(thmin,thmax);
  eplus->SetOrigin(0., 0., 0.);    // primary vertex position
  eplus->SetSigma(0.,0.,0.);       // sigma (x,y,z) in cm, on position

  // mu- generator
  AliGenBox *muminus= new AliGenBox(Npart);
  muminus->SetPart(13);
  muminus->SetMomentumRange(momen,momen);
  muminus->SetPhiRange(0,360);
  muminus->SetThetaRange(thmin,thmax);
  muminus->SetOrigin(0., 0., 0.);    // primary vertex position
  muminus->SetSigma(0.,0.,0.);       // sigma (x,y,z) in cm, on position

  // mu+ generator
  AliGenBox *muplus= new AliGenBox(Npart);
  muplus->SetPart(-13);
  muplus->SetMomentumRange(momen,momen);
  muplus->SetPhiRange(0,360);
  muplus->SetThetaRange(thmin,thmax);
  muplus->SetOrigin(0., 0., 0.);    // primary vertex position
  muplus->SetSigma(0.,0.,0.);       // sigma (x,y,z) in cm, on position

  // pi- generator
  AliGenBox *piminus= new AliGenBox(Npart);
  piminus->SetPart(-211);
  piminus->SetMomentumRange(momen,momen);
  piminus->SetPhiRange(0,360);
  piminus->SetThetaRange(thmin,thmax);
  piminus->SetOrigin(0., 0., 0.);    // primary vertex position
  piminus->SetSigma(0.,0.,0.);       // sigma (x,y,z) in cm, on position

  // pi+ generator
  AliGenBox *piplus= new AliGenBox(Npart);
  piplus->SetPart(211);
  piplus->SetMomentumRange(momen,momen);
  piplus->SetPhiRange(0,360);
  piplus->SetThetaRange(thmin,thmax);
  piplus->SetOrigin(0., 0., 0.);    // primary vertex position
  piplus->SetSigma(0.,0.,0.);       // sigma (x,y,z) in cm, on position

  // K- generator
  AliGenBox *kminus= new AliGenBox(Npart);
  kminus->SetPart(-321);
  kminus->SetMomentumRange(momen,momen);
  kminus->SetPhiRange(0,360);
  kminus->SetThetaRange(thmin,thmax);
  kminus->SetOrigin(0., 0., 0.);    // primary vertex position
  kminus->SetSigma(0.,0.,0.);       // sigma (x,y,z) in cm, on position

  // K+ generator
  AliGenBox *kplus= new AliGenBox(Npart);
  kplus->SetPart(321);
  kplus->SetMomentumRange(momen,momen);
  kplus->SetPhiRange(0,360);
  kplus->SetThetaRange(thmin,thmax);
  kplus->SetOrigin(0., 0., 0.);    // primary vertex position
  kplus->SetSigma(0.,0.,0.);       // sigma (x,y,z) in cm, on position

  // p generator
  AliGenBox *proton= new AliGenBox(Npart);
  proton->SetPart(2212);
  proton->SetMomentumRange(momen,momen);
  proton->SetPhiRange(0,360);
  proton->SetThetaRange(thmin,thmax);
  proton->SetOrigin(0., 0., 0.);    // primary vertex position
  proton->SetSigma(0.,0.,0.);       // sigma (x,y,z) in cm, on position

  // anti-p generator
  AliGenBox *aproton= new AliGenBox(Npart);
  aproton->SetPart(-2212);
  aproton->SetMomentumRange(momen,momen);
  aproton->SetPhiRange(0,360);
  aproton->SetThetaRange(thmin,thmax);
  aproton->SetOrigin(0., 0., 0.);    // primary vertex position
  aproton->SetSigma(0.,0.,0.);       // sigma (x,y,z) in cm, on position



  gener->AddGenerator(eminus,"e-",1);
  gener->AddGenerator(eplus,"e+",1);
  gener->AddGenerator(muminus,"mu-",1);
  gener->AddGenerator(muplus,"mu+",1);
  gener->AddGenerator(piminus,"pi-",1);
  gener->AddGenerator(piplus,"pi+",1);
  gener->AddGenerator(kminus,"K-",1);
  gener->AddGenerator(kplus,"K+",1);
  gener->AddGenerator(proton,"p",1);
  gener->AddGenerator(aproton,"anti-p",1);

  gener->Init();

  // FIELD
  //    
  if (mag == k2kG) {
    comment = comment.Append(" | L3 field 0.2 T");
  } else if (mag == k4kG) {
    comment = comment.Append(" | L3 field 0.4 T");
  } else if (mag == k5kG) {
    comment = comment.Append(" | L3 field 0.5 T");
  }
  printf("\n \n Comment: %s \n \n", comment.Data());
    
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., mag);
  field->SetL3ConstField(0); //Using const. field in the barrel
  rl->CdGAFile();
  gAlice->SetField(field);    

  Int_t iABSO  = 1;
  Int_t iCRT   = 0;
  Int_t iDIPO  = 0;
  Int_t iEMCAL = 0;
  Int_t iFMD   = 0;
  Int_t iFRAME = 1;
  Int_t iHALL  = 1;
  Int_t iITS   = 1;
  Int_t iMAG   = 1;
  Int_t iMUON  = 0;
  Int_t iPHOS  = 1;
  Int_t iPIPE  = 1;
  Int_t iPMD   = 0;
  Int_t iHMPID  = 1;
  Int_t iSHIL  = 0;
  Int_t iT0 = 0;
  Int_t iTOF   = 1;
  Int_t iTPC   = 1;
  Int_t iTRD   = 1;
  Int_t iVZERO = 0;
  Int_t iZDC   = 0;


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
	FRAME->SetHoles(0);
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
	ITS->SetReadDet(kTRUE);	  // don't touch this parameter if you're not an ITS developer
    //    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
	ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
	ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
	ITS->SetThicknessChip1(200.);  // chip thickness on layer 1 must be in the range [150,300]
	ITS->SetThicknessChip2(200.);  // chip thickness on layer 2 must be in the range [150,300]
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
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "SHISH");
    }

     if (iCRT)
    {
        //=================== CRT parameters ============================
        AliCRT *CRT = new AliCRTv0("CRT", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== CRT parameters ============================
        AliVZERO *VZERO = new AliVZEROv6("VZERO", "normal VZERO");
    }
}


void ProcessEnvironmentVars()
{
    // Run type
    if (gSystem->Getenv("CONFIG_RUN_TYPE")) {
      for (Int_t iRun = 0; iRun < kRunMax; iRun++) {
	if (strcmp(gSystem->Getenv("CONFIG_RUN_TYPE"), pprRunName[iRun])==0) {
	  proc = (PDC06Proc_t)iRun;
	  cout<<"Run type set to "<<pprRunName[iRun]<<endl;
	}
      }
    }

    // Random Number seed
    if (gSystem->Getenv("CONFIG_SEED")) {
      seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    }
}



