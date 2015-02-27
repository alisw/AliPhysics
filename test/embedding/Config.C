// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [1] gSystem->Load("libpythia6.4.25.so");
// root [2] gSystem->Load("libqpythia.so");
// root [3] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE/geant3/TGeant3");
// root [4] AliSimulation sim
// root [5] sim.SetConfigFile("Config.C++")
// root [6] sim.Run()

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include <TGeoGlobalMagField.h>
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliConfig.h"
#include "AliDecayerPythia.h"
#include "AliGenPythia.h"
#include "AliGenDPMjet.h"
#include "AliMagF.h"
#include "AliBODY.h"
#include "AliMAG.h"
#include "AliABSOv3.h"
#include "AliDIPOv3.h"
#include "AliHALLv3.h"
#include "AliFRAMEv2.h"
#include "AliSHILv3.h"
#include "AliPIPEv3.h"
#include "AliITSv11.h"
#include "AliTPCv2.h"
#include "AliTOFv6T0.h"
#include "AliHMPIDv3.h"
#include "AliZDCv4.h"
#include "AliTRDv1.h"
#include "AliTRDgeometry.h"
#include "AliFMDv1.h"
#include "AliMUONv1.h"
#include "AliPHOSv1.h"
#include "AliPMDv1.h"
#include "AliT0v1.h"
#include "AliEMCALv2.h"
#include "AliACORDEv1.h"
#include "AliVZEROv7.h"
#include "AliSimulation.h"
#include "AliGenCocktail.h"
#include "AliGenTherminator.h"
#include "AliGenerator.h"
#endif


enum PDC06Proc_t 
{
  kPythia6, kPhojet, kRunMax
};

const char * pprRunName[] = {
  "kPythia6", "kPhojet"
};

//--- Magnetic Field ---
enum Mag_t
{
  kNoField, k5kG, kFieldMax
};

const char * pprField[] = {
  "kNoField", "k5kG"
};

//--- Embedding run ---
enum EmbedRun_t
{
  kBackground, kMerged, kSignal, kEmbRunMax
};

const char *embedRun[] = {
  "kBackground", "kMerged", "kSignal"
};

//--- Functions ---
class AliGenPythia;
AliGenerator *MbPythia();
AliGenerator *MbPhojet();
void ProcessEnvironmentVars();

// Geterator, field, beam energy
static PDC06Proc_t   proc     = kPhojet;
static Mag_t         mag      = k5kG;
static Float_t       energy   = 10000; // energy in CMS
static EmbedRun_t    embedrun = kBackground;
//========================//
// Set Random Number seed //
//========================//
TDatime dt;
static UInt_t seed    = dt.Get();

// Comment line
static TString comment;

Float_t EtaToTheta(Float_t arg);

void Config()
{
    

  // Get settings from environment variables
  ProcessEnvironmentVars();

  gRandom->SetSeed(seed);
  cerr<<"Seed for random number generation= "<<seed<<endl; 

  // Libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libgeant321");
  gSystem->Load("libTTherminator");
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
  rl->SetNumberOfEventsPerFile(3);
  gAlice->SetRunLoader(rl);
  
  
  // Set the trigger configuration
  if ((embedrun == kBackground) || (embedrun == kMerged)) {
    AliSimulation::Instance()->SetTriggerConfig("Pb-Pb");
    cout<<"Trigger configuration is set to  Pb-Pb"<<endl;
  }
  else {
    // Set the trigger configuration: proton-proton
    AliSimulation::Instance()->SetTriggerConfig("p-p");
  }

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
  
  //======================//
  // Set External decayer //
  //======================//
  TVirtualMCDecayer* decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
  
  AliGenerator* gener = 0x0;
  if ((embedrun == kMerged) || (embedrun == kSignal)) {
    //=========================//
    // Generator Configuration //
    //=========================//
    
    if (proc == kPythia6) {
      gener = MbPythia();
    } else if (proc == kPhojet) {
      gener = MbPhojet();
    }
  }
  else {
    AliGenCocktail *gener = new AliGenCocktail();
    gener->SetPhiRange(0, 360);
    // Set pseudorapidity range from -8 to 8.
    Float_t thmin = EtaToTheta(1);   // theta min. <---> eta max
    Float_t thmax = EtaToTheta(-1);  // theta max. <---> eta min 
    gener->SetThetaRange(thmin,thmax);
    gener->SetProjectile("A",208,82);
    gener->SetTarget("A",208,82);

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

    gener->AddGenerator(genther, "THERMINATOR LHYQUID3D", 1);
  }
  
  

  // PRIMARY VERTEX
  //
  gener->SetOrigin(0., 0., 0.);    // vertex position
  //
  //
  // Size of the interaction diamond
  // Longitudinal
  Float_t sigmaz;

  if (embedrun == kBackground) {
    sigmaz  = 7.55 / TMath::Sqrt(2.); // [cm]
  }
  else {
    Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
    if (energy == 900)
      sigmaz  = 10.5 / TMath::Sqrt(2.); // [cm]
  }

  //
  // Transverse
  Float_t betast  = 10;                 // beta* [m]
  Float_t eps     = 3.75e-6;            // emittance [m]
  Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]
  Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
  printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);
    
  gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetCutVertexZ(3.);        // Truncate at 3 sigma
  gener->SetVertexSmear(kPerEvent);

  gener->Init();

  // FIELD
  //
// Field

    //  AliMagF* field = 0x0;
  if (mag == kNoField) {
    comment = comment.Append(" | L3 field 0.0 T");
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 0., 0., AliMagF::k5kGUniform));
  } else if (mag == k5kG) {
    comment = comment.Append(" | L3 field 0.5 T");
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));
  }
  printf("\n \n Comment: %s \n \n", comment.Data());
  //  TGeoGlobalMagField::Instance()->SetField(field);
    
  rl->CdGAFile();
  
  Int_t iABSO  = 1;
  Int_t iACORDE= 0;
  Int_t iDIPO  = 1;
  Int_t iEMCAL = 1;
  Int_t iFMD   = 1;
  Int_t iFRAME = 1;
  Int_t iHALL  = 1;
  Int_t iITS   = 1;
  Int_t iMAG   = 1;
  Int_t iMUON  = 1;
  Int_t iPHOS  = 1;
  Int_t iPIPE  = 1;
  Int_t iPMD   = 1;
  Int_t iHMPID = 1;
  Int_t iSHIL  = 1;
  Int_t iT0    = 1;
  Int_t iTOF   = 1;
  Int_t iTPC   = 1;
  Int_t iTRD   = 1;
  Int_t iVZERO = 1;
  Int_t iZDC   = 1;
  

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

    if (iPHOS)
    {
        //=================== PHOS parameters ===========================
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
      comment = comment.Append(" pp at 14 TeV: Pythia low-pt");
//
//    DPMJET
#if defined(__CINT__)
#endif
      gSystem->Load("libDPMJET");      // Parton density functions
      gSystem->Load("libTDPMjet");      // Parton density functions

      AliGenDPMjet* dpmjet = new AliGenDPMjet(-1); 
      dpmjet->SetMomentumRange(0, 999999.);
      dpmjet->SetThetaRange(0., 180.);
      dpmjet->SetYRange(-12.,12.);
      dpmjet->SetPtRange(0,1000.);
      dpmjet->SetProcess(kDpmMb);
      dpmjet->SetEnergyCMS(energy);

      return dpmjet;
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

    // Field
    if (gSystem->Getenv("CONFIG_FIELD")) {
      for (Int_t iField = 0; iField < kFieldMax; iField++) {
	if (strcmp(gSystem->Getenv("CONFIG_FIELD"), pprField[iField])==0) {
	  mag = (Mag_t)iField;
	  cout<<"Field set to "<<pprField[iField]<<endl;
	}
      }
    }

    // Energy
    if (gSystem->Getenv("CONFIG_ENERGY")) {
      energy = atoi(gSystem->Getenv("CONFIG_ENERGY"));
      cout<<"Energy set to "<<energy<<" GeV"<<endl;
    }

    // Random Number seed
    if (gSystem->Getenv("CONFIG_SEED")) {
      seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    }

    // Embedding job type
    if (gSystem->Getenv("CONFIG_EMBEDDING")) {
      for (Int_t iEmb = 0; iEmb < kEmbRunMax; iEmb++) {
	if (strcmp(gSystem->Getenv("CONFIG_EMBEDDING"), embedRun[iEmb])==0) {
	  embedrun = (EmbedRun_t)iEmb;
	  cout<<"Embedding run set to "<<embedRun[embedrun]<<endl;
	}
      }
      
    }

}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
