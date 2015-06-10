//
// A Config.C for ITSU pileup studies in pp. 
// ( based of the Config.C from $ALICE_ROOT/test/pileup/ )
//

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>

#include "../geant3/TGeant3/TGeant3TGeo.h"

#include "AliMagF.h"
#include "STEER/STEER/AliRunLoader.h"
#include "STEER/STEER/AliRun.h"
#include "STEER/STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenHIJINGpara.h"
#include "EVGEN/AliGenFixed.h"
#include "EVGEN/AliGenBox.h"
#include "EVGEN/AliGenPileup.h"
#include "PYTHIA6/AliGenPythia.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"
#include "STRUCT/AliABSOv3.h"
#include "STRUCT/AliDIPOv3.h"
#include "STRUCT/AliHALLv3.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv3.h"
#include "STRUCT/AliPIPEv3.h"
#include "STRUCT/AliPIPEupgrade.h"
#include "ITS/AliITSv11.h"
#include "ITS/UPGRADE/AliITSUv1.h"
#include "TPC/Sim/AliTPCv2.h"
#include "TOF/AliTOFv6T0.h"
#include "HMPID/AliHMPIDv3.h"
#include "ZDC/AliZDCv3.h"
#include "TRD/AliTRDv1.h"
#include "TRD/AliTRDgeometry.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "ACORDE/AliACORDEv1.h"
#include "VZERO/AliVZEROv7.h"
#endif

static Float_t energy   = 14000; // energy in CMS

AliGenerator* MbPythia(Float_t energy)
{
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1); 
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-1.5,1.5);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
      pythia->SetProjectile("p", 1, 1); 
      pythia->SetTarget("p", 1, 1); 
      
      return pythia;
}

void Config()
{
  // Set Random Number seed
  gRandom->SetSeed(1); // Set 0 to use the currecnt time


  // libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("liblhapdf");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libgeant321");
  gSystem->Load("libhijing");	
  gSystem->Load("libTHijing");
#endif
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");


  new     TGeant3TGeo("C++ Interface to Geant3");

  AliRunLoader* rl=0x0;


  rl = AliRunLoader::Open("galice.root",
			  AliConfig::GetDefaultEventFolderName(),
			  "recreate");
  if (rl == 0x0)
    {
      gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
      return;
    }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(2000);
  gAlice->SetRunLoader(rl);

  //
  // Set External decayer
  TVirtualMCDecayer *decayer = new AliDecayerPythia();

  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
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
  gMC->SetProcess("HADR",0);
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

// *********************** Pile up generation *********************
  // Minimum bias Pythia with pileup
  AliGenPileup *gener = new AliGenPileup();
  AliGenerator* pyth = MbPythia(energy);

  // Set the pileup interaction generator
  // The second argument is the pileup rate
  // in terms of event rate per bunch crossing
  gener->SetGenerator(pyth,0.05);
  // Set the beam time structure
  // Details on the syntax in STEER/AliTriggerBCMask
  gener->SetBCMask("1782(1H1L)"); // Full LHC

   // Generate the trigger interaction
  gener->GenerateTrigInteraction(kTRUE);

  //Set the ITS integration time
  gener->SetPileUpTimeWindow(30e-6);

  // PRIMARY VERTEX
  //
  //gener->SetOrigin(0., 0., 0.);    // vertex position
  //
  //
  // Size of the interaction diamond
  // Longitudinal
  Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
  //
  // Transverse
  Float_t betast  = 10;                 // beta* [m]
  Float_t eps     = 3.75e-6;            // emittance [m]
  Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]
  Float_t sigmaxy = TMath::Sqrt(eps*betast/gamma)/TMath::Sqrt(2.)*100.; // [cm]
  printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);
    
  gener->SetSigma(sigmaxy, sigmaxy, sigmaz);// Sigma (cm) of the IP position
  gener->SetCutVertexZ(3.);        // Truncate at 3 sigma
  gener->SetVertexSmear(kPerEvent);

  gener->Init();

// *********************** End of the pile up generation *********************

  // 
  // Activate this line if you want the vertex smearing to happen
  // track by track
  //
  //VertexSmear_t perTrack;
  //gener->SetVertexSmear(perTrack); 
  // Field (L3 0.5 T)
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));

  Int_t   iABSO  =  0;
  Int_t   iDIPO  =  0;
  Int_t   iFMD   =  0;
  Int_t   iFRAME =  0;
  Int_t   iHALL  =  0;
  Int_t   iITS   =  1;
  Int_t   iMAG   =  0;
  Int_t   iMUON  =  0;
  Int_t   iPHOS  =  0;
  Int_t   iPIPE  =  1;
  Int_t   iPMD   =  0;
  Int_t   iHMPID =  0;
  Int_t   iSHIL  =  0;
  Int_t   iT0    =  0;
  Int_t   iTOF   =  0;
  Int_t   iTPC   =  0;
  Int_t   iTRD   =  0;
  Int_t   iZDC   =  0;
  Int_t   iEMCAL =  0;
  Int_t   iACORDE=  0;
  Int_t   iVZERO =  0;
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

      AliPIPE *PIPE = new AliPIPEupgrade("PIPE", "Beam Pipe");
      //AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
    }
 
  if (iITS)
    {
      //=================== ITS parameters ============================
//      gROOT->ProcessLine(".x CreateITSU.C");
      gROOT->ProcessLine(".x CreateITSUv1.C");
//      gROOT->ProcessLine(".x CreateITSUv0.C");
//      gROOT->ProcessLine(".x CreateITSU_MS.C");
      //    CreateITSU();

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
      AliTRDgeometry *geoTRD = TRD->GetGeometry();
      // Partial geometry: modules at 2,3,4,6,11,12,14,15
      // starting at 6h in positive direction
      geoTRD->SetSMstatus(0,0);
      geoTRD->SetSMstatus(1,0);
      geoTRD->SetSMstatus(5,0);
      geoTRD->SetSMstatus(7,0);
      geoTRD->SetSMstatus(8,0);
      geoTRD->SetSMstatus(9,0);
      geoTRD->SetSMstatus(10,0);
      geoTRD->SetSMstatus(13,0);
      geoTRD->SetSMstatus(16,0);
      geoTRD->SetSMstatus(17,0);
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
      AliMUON *MUON = new AliMUONv1("MUON","default");
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
      AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETEV1");
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


}


