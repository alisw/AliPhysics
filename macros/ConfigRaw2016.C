//
// Configuration file for generating 2016 geometry
//

// To use this macro in compiled mode one first would need to do:
// root [0] gSystem->Load("libgeant321");
// root [1] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -Igeant3_src/TGeant3");
//
// However it is meant to be called in the following sequence:
// root [0] AliCDBManager* cdb = AliCDBManager::Instance();
// root [0] cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
// root [0] cdb->SetRun(200000);
// root [0] gROOT->LoadMacro("ConfigRaw2016.C");
// root [0] gInterpreter->ProcessLine(gAlice->GetConfigFunction());
// root [0] gAlice->GetMCApp()->Init();

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeoGlobalMagField.h>
#include <TGeant3TGeo.h>
#include "STEER/STEER/AliRunLoader.h"
#include "STEER/STEER/AliRun.h"
#include "STEER/STEER/AliConfig.h"
#include "STEER/STEER/AliSimulation.h"
#include "PYTHIA6/AliPythia6/AliDecayerPythia.h"
#include "PYTHIA6/AliPythia6/AliGenPythia.h"
#include "TDPMjet/AliGenDPMjet.h"
#include "STEER/STEERBase/AliMagF.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"
#include "STRUCT/AliABSOv3.h"
#include "STRUCT/AliDIPOv3.h"
#include "STRUCT/AliHALLv3.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv3.h"
#include "STRUCT/AliPIPEv3.h"
#include "ITS/ITSsim/AliITSv11.h"
#include "TPC/TPCsim/AliTPCv2.h"
#include "TOF/TOFsim/AliTOFv6T0.h"
#include "HMPID/HMPIDsim/AliHMPIDv3.h"
#include "ZDC/ZDCsim/AliZDCv4.h"
#include "TRD/TRDsim/AliTRDv1.h"
#include "TRD/TRDbase/AliTRDgeometry.h"
#include "FMD/FMDsim/AliFMDv1.h"
#include "MUON/MUONsim/AliMUONv1.h"
#include "PHOS/PHOSsim/AliPHOSv1.h"
#include "PHOS/PHOSbase/AliPHOSSimParam.h"
#include "PMD/PMDsim/AliPMDv1.h"
#include "T0/T0sim/AliT0v1.h"
#include "EMCAL/EMCALsim/AliEMCALv2.h"
#include "ACORDE/ACORDEsim/AliACORDEv1.h"
#include "VZERO/VZEROsim/AliVZEROv7.h"
#include "AD/ADsim/AliADv1.h"
#endif

enum PDC06Proc_t 
{
  kPythia6, kPythia6D6T, kPhojet, kRunMax
};

const char * pprRunName[] = {
    "kPythia6", "kPythia6D6T", "kPhojet" 
};

enum Mag_t
{
  kNoField, k5kG, kFieldMax
};

const char * pprField[] = {
  "kNoField", "k5kG"
};

//--- Functions ---
class AliGenPythia;
AliGenerator *MbPythia();
AliGenerator *MbPythiaTuneD6T();
AliGenerator *MbPhojet();
void ProcessEnvironmentVars();

// Geterator, field, beam energy
static PDC06Proc_t   proc     = kPhojet;
static Mag_t         mag      = k5kG;
static Float_t       energy   = 10000; // energy in CMS
//========================//
// Set Random Number seed //
//========================//
TDatime dt;
static UInt_t seed    = dt.Get();

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
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  if (proc != kPythia6D6T) {
      gSystem->Load("libpythia6");     // Pythia 6.2
  } else {
      gSystem->Load("libqpythia");     // Pythia 6.4
  }
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
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
  // gAlice->SetGeometryFromFile("geometry.root");
  // gAlice->SetGeometryFromCDB();
  
  // Set the trigger configuration: proton-proton
  AliSimulation::Instance()->SetTriggerConfig("p-p");

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
  
  if (proc == kPythia6) {
      gener = MbPythia();
  } else if (proc == kPythia6D6T) {
      gener = MbPythiaTuneD6T();
  } else if (proc == kPhojet) {
      gener = MbPhojet();
  }
  
  

  // PRIMARY VERTEX
  //
  gener->SetOrigin(0., 0., 0.);    // vertex position
  //
  //
  // Size of the interaction diamond
  // Longitudinal
  Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
  if (energy == 900)
    sigmaz  = 10.5 / TMath::Sqrt(2.); // [cm]
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
  AliMagF* field = 0x0;
  if (mag == kNoField) {
    comment = comment.Append(" | L3 field 0.0 T");
    field = new AliMagF("Maps","Maps", 0., 0., AliMagF::k5kGUniform,AliMagF::kBeamTypepp, energy/2.0);
  } else if (mag == k5kG) {
    comment = comment.Append(" | L3 field 0.5 T");
    field = new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG,	AliMagF::kBeamTypepp, energy/2.0);
  }
  printf("\n \n Comment: %s \n \n", comment.Data());

  TGeoGlobalMagField::Instance()->SetField(field);
    
  rl->CdGAFile();
  
  Int_t iABSO  = 1;
  Int_t iACORDE= 1;
  Int_t iAD    = 1;
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

        AliTPC *TPC = new AliTPCv2("TPC", "Ne-CO2");
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
	// SM indexing starts at 3h in positive direction. We assume all SM
        // are installed (i.e. all SM missing last year: 4, 5, 12, 13, 14)
        // AliTRDgeometry *geoTRD = TRD->GetGeometry();
	// geoTRD->SetSMstatus(4,0);
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

        AliPHOS *PHOS = new AliPHOSv1("PHOS", "Run2");
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

        AliEMCAL* EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETE12SMV1_DCAL_8SM", kFALSE); 
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

    if (iAD)
    {
        //=================== AD parameters ============================

        AliAD *AD = new AliADv1("AD", "The AD geometry");
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

AliGenerator* MbPythiaTuneD6T()
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
//    Tune
//    109     D6T : Rick Field's CDF Tune D6T (NB: needs CTEQ6L pdfs externally)
      pythia->SetTune(109);
      pythia->SetStrucFunc(kCTEQ6l);
//
      return pythia;
}

AliGenerator* MbPhojet()
{
      comment = comment.Append(" pp at 14 TeV: Pythia low-pt");
//
//    DPMJET
#if defined(__CINT__)
  gSystem->Load("libdpmjet");      // Parton density functions
  gSystem->Load("libTDPMjet");      // Parton density functions
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
}
