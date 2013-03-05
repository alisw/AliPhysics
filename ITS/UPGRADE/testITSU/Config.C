/// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config.C++"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenHIJINGpara.h"
#include "EVGEN/AliGenFixed.h"
#include "EVGEN/AliGenBox.h"
#include "STEER/AliMagWrapCheb.h"
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
#include "ITS/UPGRADE/AliITSUv11.h"
#include "TPC/AliTPCv2.h"
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
#include <TVirtualMagField.h>
#endif

Int_t generatorFlag = 0;

/* $Id: Config.C 47147 2011-02-07 11:46:44Z amastros $ */
enum PprTrigConf_t
{
  kDefaultPPTrig, kDefaultPbPbTrig
};

const char * pprTrigConfName[] = {
  "p-p","Pb-Pb"
};

Float_t EtaToTheta(Float_t arg);

static PprTrigConf_t strig = kDefaultPPTrig;// default PP trigger configuration

void Config()
{
  // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
  // Theta range given through pseudorapidity limits 22/6/2001

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
  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeSim.so");


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

  // Set the trigger configuration
  // gAlice->SetTriggerDescriptor(pprTrigConfName[strig]);
  //cout<<"Trigger configuration is set to  "<<pprTrigConfName[strig]<<endl;
  AliSimulation::Instance()->SetTriggerConfig(pprTrigConfName[strig]);
  cout<<"Trigger configuration is set to  pprTrigConfName[strig] "<<endl;

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

  // Special generation for Valgrind tests
  // Each detector is fired by few particles selected 
  // to cover specific cases


  // The cocktail itself
  
  if (generatorFlag==0) {
    // Fast generator with parametrized pi,kaon,proton distributions
    
    int     nParticles = 30;//14022;
    AliGenHIJINGpara *gener = new AliGenHIJINGpara(nParticles);
    gener->SetMomentumRange(0.1, 10.);
    gener->SetPhiRange(0., 360.);
    Float_t thmin = EtaToTheta(2.5);   // theta min. <---> eta max
    Float_t thmax = EtaToTheta(-2.5);  // theta max. <---> eta min
    gener->SetThetaRange(thmin,thmax);
    gener->SetOrigin(0, 0, 0);  //vertex position
    gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position
    gener->Init();
    
  } else if (generatorFlag==1) {
    
    // Pure HiJing generator adapted to ~2000dNdy at highest energy
    
    AliGenHijing *generHijing = new AliGenHijing(-1);
    generHijing->SetEnergyCMS(5500.); // GeV
    generHijing->SetImpactParameterRange(0,2);
    generHijing->SetReferenceFrame("CMS");
    generHijing->SetProjectile("A", 208, 82);
    generHijing->SetTarget    ("A", 208, 82);
    generHijing->KeepFullEvent();
    generHijing->SetJetQuenching(1);
    generHijing->SetShadowing(1);
    generHijing->SetSpectators(0);
    generHijing->SetSelectAll(0);
    generHijing->SetPtHardMin(4.5);
    
    AliGenerator*  gener = generHijing;
    gener->SetSigma(0, 0, 6);      // Sigma in (X,Y,Z) (cm) on IP position
    gener->SetVertexSmear(kPerEvent);
    gener->Init();
      
  }

  // 
  // Activate this line if you want the vertex smearing to happen
  // track by track
  //
  //VertexSmear_t perTrack;
  //gener->SetVertexSmear(perTrack); 
  // Field (L3 0.5 T)
  //AliMagF* field = new AliMagF("map","map",2, -1.,1., 15, AliMagF::k5kGUniform);
  //TGeoGlobalMagField::Instance()->SetField(field);
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
  Int_t   iTPC   =  1;
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

      AliPIPE *PIPE = new AliPIPEupgrade("PIPE", "Beam Pipe",0,1.8,0.08,40.0);
      //AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
    }
 
  if (iITS)
    {
      //=================== ITS parameters ============================
      //
      // create segmentations:
      AliITSUSegmentationPix* seg0 = new AliITSUSegmentationPix(0,    // segID (0:9)
								5,    // chips per module
								1500, // ncols (total for module)
								350,  //835,  // nrows
								33.e-4,  // default row pitch in cm
								20.e-4,  // default col pitch in cm
								18.e-4  // sensor thickness in cm
								);    // see AliITSUSegmentationPix.h for extra options
      seg0->Store(AliITSUGeomTGeo::GetITSsegmentationFileName());
      AliITSUSegmentationPix* seg1 = new AliITSUSegmentationPix(1,    // segID (0:9)
								5*2,    // chips per module
								1500, // ncols (total for module)
								700,//835,  // nrows
								33.e-4,  // default row pitch in cm
								20.e-4,  // default col pitch in cm
								18.e-4  // sensor thickness in cm
								);    // see AliITSUSegmentationPix.h for extra options
      seg1->Store(AliITSUGeomTGeo::GetITSsegmentationFileName());
      AliITSUSegmentationPix* seg2 = new AliITSUSegmentationPix(2,    // segID (0:9)
								5*2,    // chips per module
								1500, // ncols (total for module)
								700,//835,  // nrows
								33.e-4,  // default row pitch in cm
								20.e-4,  // default col pitch in cm
								18.e-4   // sensor thickness in cm
								);    // see AliITSUSegmentationPix.h for extra options
      seg2->Store(AliITSUGeomTGeo::GetITSsegmentationFileName());
      //
      int nmod,nlad; // modules per ladded, n ladders
      // sum of insensitive boarder around module (in cm)
      Float_t deadX = 0.1;  // on each side
      Float_t deadZ = 0.01; // on each side
      double tilt = 10.;     double thickLr = 0.05;
      //      virtual void   DefineLayerTurbo(const Int_t nlay, const Double_t r,  const Double_t zlen, const Int_t nladd,   const Int_t nmod, const Double_t width,
      //				  const Double_t tilt,   const Double_t lthick = 0.,    const Double_t dthick = 0.,   const UInt_t detType=0);
      AliITSUv11 *ITS  = new AliITSUv11("ITS Upgrade",7);
      nmod = 9;
      nlad = 12;
      ITS->DefineLayerTurbo(0,0., 2.2,  nmod*(seg0->Dz()+deadZ*2), nlad, nmod, seg0->Dx()+deadX*2, tilt, thickLr, seg0->Dy(), seg0->GetDetTypeID());
      nmod = 9;
      nlad = 16;
      ITS->DefineLayerTurbo(1,0., 2.8,  nmod*(seg0->Dz()+deadZ*2), nlad, nmod, seg0->Dx()+deadX*2, tilt, thickLr, seg0->Dy(), seg0->GetDetTypeID());
      nmod = 9;
      nlad = 20;
      ITS->DefineLayerTurbo(2,0., 3.6,  nmod*(seg0->Dz()+deadZ*2), nlad, nmod, seg0->Dx()+deadX*2, tilt, thickLr, seg0->Dy(), seg0->GetDetTypeID());
      nmod = 29;
      nlad = 48;
      ITS->DefineLayerTurbo(3,0., 20.0, nmod*(seg1->Dz()+deadZ*2), nlad, nmod, seg1->Dx()+deadX*2, tilt, thickLr, seg1->Dy(), seg1->GetDetTypeID());
      nmod = 29;
      nlad = 48;
      ITS->DefineLayerTurbo(4,0., 22.0, nmod*(seg1->Dz()+deadZ*2), nlad, nmod, seg1->Dx()+deadX*2, tilt, thickLr, seg1->Dy(), seg1->GetDetTypeID());
      nmod = 50;
      nlad = 94;
      ITS->DefineLayerTurbo(5,0., 40.0, nmod*(seg2->Dz()+deadZ*2), nlad, nmod, seg2->Dx()+deadX*2, tilt, thickLr, seg2->Dy(), seg2->GetDetTypeID()); //41 creates ovl!
      nmod = 50;
      nlad = 94;
      ITS->DefineLayerTurbo(6,0., 43.0, nmod*(seg2->Dz()+deadZ*2), nlad, nmod, seg2->Dx()+deadX*2, tilt, thickLr, seg2->Dy(), seg2->GetDetTypeID()); 
      //

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

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
