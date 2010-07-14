//
// Configuration for the first physics production 2008
//

// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config.C++")

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "PYTHIA6/AliGenPythia.h"
#include "TDPMjet/AliGenDPMjet.h"
#include "STEER/AliMagFCheb.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"
#include "STRUCT/AliABSOv3.h"
#include "STRUCT/AliDIPOv3.h"
#include "STRUCT/AliHALLv3.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv3.h"
#include "STRUCT/AliPIPEv3.h"
#include "ITS/AliITSv11Hybrid.h"
#include "TPC/AliTPCv2.h"
#include "TOF/AliTOFv6T0.h"
#include "HMPID/AliHMPIDv3.h"
#include "ZDC/AliZDCv3.h"
#include "TRD/AliTRDv1.h"
#include "TRD/AliTRDgeometry.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PHOS/AliPHOSSimParam.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "ACORDE/AliACORDEv1.h"
#include "VZERO/AliVZEROv7.h"
#endif

enum PDC06Proc_t
{
  kPythia6, kPhojet, kTherminator, kRunMax
};

const char  * pprRunName[] =
{
  "kPythia6", "kPhojet", "kTherminator"
};

enum Mag_t
{
  kNoField, k5kG, kFieldMax
};

const char * pprField[] =
{
"kNoField", "k5kG"
};

enum PprTrigConf_t
{
kDefaultPPTrig, kDefaultPbPbTrig
};

const char * pprTrigConfName[] =
{
"p-p","Pb-Pb"
};

static PprTrigConf_t strig = kDefaultPbPbTrig;// default PP trigger configuration

//--- Functions ---
class AliGenPythia;
AliGenerator *MbPythia();
AliGenerator *MbPhojet();
AliGenerator *Therminator();
void ProcessEnvironmentVars();

// Geterator, field, beam energy
static PDC06Proc_t proc = kTherminator;
static Mag_t mag = k5kG;
static Float_t energy = 10000; // energy in CMS
//======================= =//
// Set Random Number seed //
//========================//
TDatime dt;
static UInt_t seed = dt.Get();

// Comment line
static TString comment;

void Config()
{

  // Get settings from environment variables
  ProcessEnvironmentVars();

  gRandom->SetSeed(seed);
  cerr<<"Seed for random number generation= "<<seed<<endl;

  // Libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("liblhapdf"); // Parton density functions
  gSystem->Load("libEGPythia6"); // TGenerator interface
  gSystem->Load("libpythia6"); // Pythia
  gSystem->Load("libAliPythia6"); // ALICE specific implementations
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

  // Set the trigger configuration: proton-proton
  gAlice->SetTriggerDescriptor(pprTrigConfName[strig]);
  cout<<"Trigger configuration is set to  "<<pprTrigConfName[strig]<<endl;

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

  Float_t cut = 1.e-3; // 1MeV cut by default
  Float_t tofmax = 1.e10;

  gMC->SetCut("CUTGAM", cut);
  gMC->SetCut("CUTELE", cut);
  gMC->SetCut("CUTNEU", cut);
  gMC->SetCut("CUTHAD", cut);
  gMC->SetCut("CUTMUO", cut);
  gMC->SetCut("BCUTE", cut);
  gMC->SetCut("BCUTM", cut);
  gMC->SetCut("DCUTE", cut);
  gMC->SetCut("DCUTM", cut);
  gMC->SetCut("PPCUTM", cut);
  gMC->SetCut("TOFMAX", tofmax);

  //======================//
  // Set External decayer //
  //======================//
  TVirtualMCDecayer* decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);

  //=========================//
  // Generator Configuration //
  //=========================//
  AliGenerator* gener = 0x0;

  switch (proc)
  {
    case kPythia6:
    gener = MbPythia();
    break;
    case kPhojet:
    gener = MbPhojet();
    break;
    case kTherminator:
    gener = Therminator();
    break;
  }

  // PRIMARY VERTEX
  //
  gener->SetOrigin(0., 0., 0.); // vertex position
  //
  // Size of the interaction diamond
  // Longitudinal
  Float_t sigmaz = 5.4 / TMath::Sqrt(2.); // [cm]
  if (energy == 900)
  sigmaz = 10.5 / TMath::Sqrt(2.); // [cm]
  //
  // Transverse
  Float_t betast = 10; // beta* [m]
  Float_t eps = 3.75e-6; // emittance [m]
  Float_t gamma = energy / 2.0 / 0.938272; // relativistic gamma [1]
  Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.; // [cm]
  printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);

  gener->SetSigma(sigmaxy, sigmaxy, sigmaz); // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetCutVertexZ(3.); // Truncate at 3 sigma
  gener->SetVertexSmear(kPerEvent);

  gener->Init();

  // FIELD
  AliMagF* field = 0x0;

  if (mag == kNoField)
  {
    comment = comment.Append(" | L3 field 0.0 T");
    field = new AliMagF("Maps","Maps", -1., -1., AliMagF::k2kG);
  }
  else if (mag == k5kG)
  {
    comment = comment.Append(" | L3 field 0.5 T");
    field = new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG);
  }
  printf("\n \n Comment: %s \n \n", comment.Data());

  TGeoGlobalMagField::Instance()->SetField(field);

  rl->CdGAFile();

  Int_t iABSO = 1;
  Int_t iACORDE= 0;
  Int_t iDIPO = 1;
  Int_t iEMCAL = 0;
  Int_t iFMD = 1;
  Int_t iFRAME = 1;
  Int_t iHALL = 1;
  Int_t iITS = 1;
  Int_t iMAG = 1;
  Int_t iMUON = 1;
  Int_t iPHOS = 1;
  Int_t iPIPE = 1;
  Int_t iPMD = 0;
  Int_t iHMPID = 1;
  Int_t iSHIL = 1;
  Int_t iT0 = 1;
  Int_t iTOF = 1;
  Int_t iTPC = 1;
  Int_t iTRD = 1;
  Int_t iVZERO = 1;
  Int_t iZDC = 1;

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

    AliITS *ITS = new AliITSv11Hybrid("ITS","ITS v11Hybrid");
  }

  if (iTPC)
  {
    //============================ TPC parameters =====================

    AliTPC *TPC = new AliTPCv2("TPC", "Default");
  }

  if (iTOF)
  {
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
    AliTRDgeometry *geoTRD = TRD->GetGeometry();
    // Partial geometry: modules at 0,1,7,8,9,10,17
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

  if (iPHOS)
  {
    //=================== PHOS parameters ===========================

    AliPHOS *PHOS = new AliPHOSv1("PHOS", "IHEP");
    //Set simulation parameters different from the default ones.
    AliPHOSSimParam* simEmc = AliPHOSSimParam::GetInstance();

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
AliT0  *T0 = new AliT0v1("T0", "T0 Detector");
}

if (iEMCAL)
{
  //=================== EMCAL parameters ============================

  AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETE");
}

if (iACORDE)
{
  //  =================== ACORDE parameters ============================

AliACORDE *ACORDE = new AliACORDEv1("ACORDE", "normal ACORDE");
}

if (iVZERO)
{
//=================== ACORDE parameters ============================

AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
}
}
//
//           PYTHIA
//

AliGenerator* MbPythia()
{
comment = comment.Append(" pp at 14 TeV: Pythia low-pt");
//
//    Pythia
AliGenPythia* pythia = new AliGenPythia(-1);
pythia->SetMomentumRange(0, 999999.);
pythia->SetThetaRange(0., 180.);
pythia->SetYRange(-12.,12.);
pythia->SetPtRange(0,1000.);
pythia->SetProcess(kPyMb);
pythia->SetEnergyCMS(energy);

return pythia;
}

AliGenerator* MbPhojet()
{
comment = comment.Append(" pp at 14 TeV: Phojet low-pt");
//
//    DPMJET
#if defined(__CINT__)
gSystem->Load("libdpmjet"); // Parton density functions
gSystem->Load("libTDPMjet"); // Parton density functions
#endif
AliGenDPMjet* dpmjet = new AliGenDPMjet(-1);
dpmjet->SetMomentumRange(0, 999999.);
dpmjet->SetThetaRange(0., 180.);
dpmjet->SetYRange(-12.,12.);
dpmjet->SetPtRange(0,1000.);
dpmjet->SetProcess(kDpmMb);
dpmjet->SetEnergyCMS(energy);

return dpmjet;
}

AliGenerator* Therminator()
{
comment = comment.Append(" pp at 14 TeV: Therminator");

#if defined(__CINT__)
gSystem->Load("libTTherminator"); // Therminator library
#endif
AliGenTherminator *genther = new AliGenTherminator();
genther->SetFileName("event.out");
genther->SetEventNumberInFile(1);
genther->SetTemperature(0.145);
genther->SetMiuI(-0.0009);
genther->SetMiuS(0.000);
genther->SetMiuB(0.0008);
genther->SetAlfaRange(8.0);
genther->SetRapRange(4.0);
genther->SetRhoMax(7.74);
genther->SetTau(9.74);
genther->SetModel("Lhyquid3D");
genther->SetLhyquidSet("LHC500C2030");

return genther;
}

void ProcessEnvironmentVars()
{
// Run type
if (gSystem->Getenv("CONFIG_RUN_TYPE"))
{
  for (Int_t iRun = 0; iRun < kRunMax; iRun++)
  {
    if (strcmp(gSystem->Getenv("CONFIG_RUN_TYPE"), pprRunName[iRun])==0)
    {
      proc = (PDC06Proc_t)iRun;
      cout<<"Run type set to "<<pprRunName[iRun]<<endl;
    }
  }
}

// Field
if (gSystem->Getenv("CONFIG_FIELD"))
{
  for (Int_t iField = 0; iField < kFieldMax; iField++)
  {
    if (strcmp(gSystem->Getenv("CONFIG_FIELD"), pprField[iField])==0)
    {
      mag = (Mag_t)iField;
      cout<<"Field set to "<<pprField[iField]<<endl;
    }
  }
}

// Energy
if (gSystem->Getenv("CONFIG_ENERGY"))
{
  energy = atoi(gSystem->Getenv("CONFIG_ENERGY"));
  cout<<"Energy set to "<<energy<<" GeV"<<endl;
}

// Random Number seed
if (gSystem->Getenv("CONFIG_SEED"))
{
  seed = atoi(gSystem->Getenv("CONFIG_SEED"));
}
}
