// $Id$
//
// Configuration for the Geant4 production 2010
// By E. Sicking, CERN
//

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
#include "ITS/AliITSv11.h"
#include "TPC/AliTPCv2.h"
#include "TOF/AliTOFv6T0.h"
#include "HMPID/AliHMPIDv3.h"
#include "ZDC/AliZDCv4.h"
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
    kPythia6, kPythia6D6T, kPythia6ATLAS, kPythia6ATLAS_Flat, kPythiaPerugia0, kPhojet, kRunMax
  };

const char * pprRunName[] = {
  "kPythia6", "kPythia6D6T", "kPythia6ATLAS", "kPythia6ATLAS_Flat", "kPythiaPerugia0", "kPhojet" 
};

enum Mag_t
  {
    kNoField, k5kG, kFieldMax
  };

const char * pprField[] = {
  "kNoField", "k5kG"
};

enum PhysicsList_t 
  {
    QGSP_BERT_EMV,         CHIPS,         QGSP_BERT_CHIPS,         QGSP_FTFP_BERT,         FTFP_BERT,         FTFP_BERT_EMV,
    QGSP_BERT_EMV_OPTICAL, CHIPS_OPTICAL, QGSP_BERT_CHIPS_OPTICAL, QGSP_FTFP_BERT_OPTICAL, FTFP_BERT_OPTICAL, FTFP_BERT_EMV_OPTICAL, kListMax
  };

const char * physicsListName[] = {
  "QGSP_BERT_EMV",         "CHIPS",         "QGSP_BERT_CHIPS",         "QGSP_FTFP_BERT",         "FTFP_BERT",         "FTFP_BERT_EMV",
  "QGSP_BERT_EMV_OPTICAL", "CHIPS_OPTICAL", "QGSP_BERT_CHIPS_OPTICAL", "QGSP_FTFP_BERT_OPTICAL", "FTFP_BERT_OPTICAL", "FTFP_BERT_EMV_OPTICAL",
};

enum PprTrigConf_t
  {
    kDefaultPPTrig, kDefaultPbPbTrig
  };

const char * pprTrigConfName[] = {
  "p-p","Pb-Pb"
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
static Int_t         runNumber = 0;
static PprTrigConf_t strig = kDefaultPPTrig; // default pp trigger configuration
static PhysicsList_t physicslist  = QGSP_BERT_EMV;

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
  if (proc == kPythia6 || proc == kPhojet) {
    gSystem->Load("libpythia6");        // Pythia 6.2
  } else {
    gSystem->Load("libpythia6.4.21");   // Pythia 6.4
  }
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  // gSystem->Load("libgeant321");
#endif

  // new TGeant3TGeo("C++ Interface to Geant3");

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

  AliSimulation::Instance()->SetTriggerConfig(pprTrigConfName[strig]);
  cout <<"Trigger configuration is set to  "<<pprTrigConfName[strig]<<endl;


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
      TPC->SetPrimaryIonisation();// not used with Geant3
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
      //Collimators aperture
      ZDC->SetVCollSideCAperture(0.85);
      ZDC->SetVCollSideCCentre(0.);
      ZDC->SetVCollSideAAperture(0.75);
      ZDC->SetVCollSideACentre(0.);
      //Detector position
      ZDC->SetYZNC(1.6);
      ZDC->SetYZNA(1.6);
      ZDC->SetYZPC(1.6);
      ZDC->SetYZPA(1.6);
    }

  if (iTRD)
    {
      //=================== TRD parameters ============================

      AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
      AliTRDgeometry *geoTRD = TRD->GetGeometry();
      // Partial geometry: modules at 0,1,7,8,9,16,17
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
      // activate trigger efficiency by cells
      MUON->SetTriggerEffCells(1);
    }

  if (iPHOS)
    {
      //=================== PHOS parameters ===========================

      AliPHOS *PHOS = new AliPHOSv1("PHOS", "noCPV_Modules123");

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

      AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_FIRSTYEARV1");
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



  // Load Geant4 + Geant4 VMC libraries
  //
  std::string g4libmacro("$G4VMCINSTALL/share/examples/macro/g4libs.C");
  if (gClassTable->GetID("TGeant4") == -1) {
    // Load Geant4 libraries
    if (!gInterpreter->IsLoaded(g4libmacro.c_str())) {
      gROOT->LoadMacro(g4libmacro.c_str());
      gInterpreter->ProcessLine("g4libs()");
    }
  }


  // Create Geant4 VMC
  //  
  TGeant4 *geant4 = 0;
  if ( ! gMC ) {
    TG4RunConfiguration* runConfiguration=0x0;
    for (Int_t iList = 0; iList < kListMax; iList++) {
      if(iList<kListMax/2){
	if(physicslist == iList){
	  runConfiguration = 
	    new TG4RunConfiguration("geomRoot", 
				    physicsListName[iList], 
				    "specialCuts+stackPopper+stepLimiter",
				    true);
	}
      }
      else if(iList>=kListMax/2){//add "optical" PL to HadronPhysicsList
	if(physicslist == iList){
	  runConfiguration = 
	    new TG4RunConfiguration("geomRoot", 
				    Form("%s+optical",physicsListName[iList-kListMax/2]), 
				    "specialCuts+stackPopper+stepLimiter",
				    true);
	}
      }
    }
    geant4 = new TGeant4("TGeant4", "The Geant4 Monte Carlo", runConfiguration);
    cout << "Geant4 has been created." << endl;
  } 
  else {
    cout << "Monte Carlo has been already created." << endl;
  }  
  
  
  
  // Customization of Geant4 VMC
  //
  geant4->ProcessGeantCommand("/mcVerbose/all 1");  
  geant4->ProcessGeantCommand("/mcVerbose/geometryManager 1");  
  geant4->ProcessGeantCommand("/mcVerbose/opGeometryManager 1");  
  geant4->ProcessGeantCommand("/mcTracking/loopVerbose 1");     
  geant4->ProcessGeantCommand("/mcPhysics/rangeCuts 0.01 mm"); 
  // for Geant4 <= 9.4.p03
  //geant4->ProcessGeantCommand("/mcPhysics/selectOpProcess Scintillation");
  //geant4->ProcessGeantCommand("/mcPhysics/setOpProcessActivation false");
  // for Geant4 >= 9.5
  geant4->ProcessGeantCommand("/optics_engine/selectOpProcess Scintillation");
  geant4->ProcessGeantCommand("/optics_engine/setOpProcessUse false");
  geant4->ProcessGeantCommand("/optics_engine/selectOpProcess OpWLS");
  geant4->ProcessGeantCommand("/optics_engine/setOpProcessUse false");
  geant4->ProcessGeantCommand("/optics_engine/selectOpProcess OpMieHG");
  geant4->ProcessGeantCommand("/optics_engine/setOpProcessUse false");
  
  geant4->ProcessGeantCommand("/mcVerbose/composedPhysicsList 2");  
  geant4->ProcessGeantCommand("/mcTracking/skipNeutrino true");
  // geant4->ProcessGeantCommand("/mcDet/setMaxStepInLowDensityMaterials 1 cm");


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
  } else if (proc == kPythia6ATLAS) {
    gener = MbPythiaTuneATLAS();
  } else if (proc == kPythiaPerugia0) {
    gener = MbPythiaTunePerugia0();
  } else if (proc == kPythia6ATLAS_Flat) {
    gener = MbPythiaTuneATLAS_Flat();
  } else if (proc == kPhojet) {
    gener = MbPhojet();
  }
  
  
  //
  //
  // Size of the interaction diamond
  // Longitudinal
  Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
  if (energy == 900)
    //sigmaz  = 10.5 / TMath::Sqrt(2.); // [cm]
    //sigmaz = 3.7;
    if (energy == 7000)
      sigmaz  = 6.3 / TMath::Sqrt(2.); // [cm]
  
  //
  // Transverse

  // beta*
  Float_t betast                  = 10.0;  // beta* [m]
  if (runNumber >= 117048) betast =  2.0;
  if (runNumber >  122375) betast =  3.5;  // starting with fill 1179
  //	
  //
  Float_t eps     = 5.0e-6;            // emittance [m]
  Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]
  Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
  printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);
    
  gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetVertexSmear(kPerEvent);
  gener->Init();

  printf("\n \n Comment: %s \n \n", comment.Data());


}
//
//           PYTHIA
//

AliGenerator* MbPythia()
{
  comment = comment.Append(" pp: Pythia low-pt");
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
  comment = comment.Append(" pp: Pythia low-pt");
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
  pythia->SetTune(109); // F I X 
  pythia->SetStrucFunc(kCTEQ6l);
  //
  return pythia;
}

AliGenerator* MbPythiaTunePerugia0()
{
  comment = comment.Append(" pp: Pythia low-pt (Perugia0)");
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
  //    320     Perugia 0
  pythia->SetTune(320); 
  pythia->UseNewMultipleInteractionsScenario();
  pythia->SetCrossingAngle(0,0.000280);

  //
  return pythia;
}


AliGenerator* MbPythiaTuneATLAS()
{
  comment = comment.Append(" pp: Pythia low-pt");
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
  //    C   306 ATLAS-CSC: Arthur Moraes' (new) ATLAS tune (needs CTEQ6L externally)
  pythia->SetTune(306);
  pythia->SetStrucFunc(kCTEQ6l);
  //
  return pythia;
}

AliGenerator* MbPythiaTuneATLAS_Flat()
{
  AliGenPythia* pythia = MbPythiaTuneATLAS();
      
  comment = comment.Append("; flat multiplicity distribution");
      
  // set high multiplicity trigger
  // this weight achieves a flat multiplicity distribution
  Double_t cont[] =
    {0, 
     0.234836, 0.103484, 0.00984802, 0.0199906, 0.0260018, 0.0208481, 0.0101477, 0.00146998, -0.00681513, -0.0114928,
     -0.0113352, -0.0102012, -0.00895238, -0.00534961, -0.00261648, -0.000819048, 0.00130816, 0.00177978, 0.00373838, 0.00566255,
     0.00628156, 0.00687458, 0.00668017, 0.00702917, 0.00810848, 0.00876167, 0.00832783, 0.00848518, 0.0107709, 0.0106849,
     0.00933702, 0.00912525, 0.0106553, 0.0102785, 0.0101756, 0.010962, 0.00957103, 0.00970448, 0.0117133, 0.012271,
     0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113,
     0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113,
     0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113,
     0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113,
     0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113,
     0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     0};

  Double_t err[] =
    {0, 
     0.00168216, 0.000743548, 0.00103125, 0.00108605, 0.00117101, 0.00124577, 0.00129119, 0.00128341, 0.00121421, 0.00112431,
     0.00100588, 0.000895078, 0.000790314, 0.000711673, 0.000634547, 0.000589133, 0.000542763, 0.000509552, 0.000487375, 0.000468906, 
     0.000460196, 0.000455806, 0.00044843, 0.000449317, 0.00045007, 0.000458016, 0.000461167, 0.000474742, 0.00050161, 0.00051637, 
     0.000542336, 0.000558854, 0.000599169, 0.000611982, 0.000663855, 0.000696322, 0.000722825, 0.000771686, 0.000838023, 0.000908317, 
     0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003,
     0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003,
     0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003,
     0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003,
     0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003,
     0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0};

  TH1F *weight = new TH1F("newweight","newweight",120,-0.5,119.5);

  weight->SetContent(cont);
  weight->SetError(err);
        
  Int_t limit = weight->GetRandom();
  pythia->SetTriggerChargedMultiplicity(limit, 1.4);
      
  comment = comment.Append(Form("; multiplicity threshold set to %d in |eta| < 1.4", limit));

  return pythia;
}

AliGenerator* MbPhojet()
{
  comment = comment.Append(" pp: Pythia low-pt");
  //
  //    DPMJET
#if defined(__CINT__)
  gSystem->Load("libDPMJET");      // Parton density functions
  gSystem->Load("libTDPMjet");      // Parton density functions
#endif
  AliGenDPMjet* dpmjet = new AliGenDPMjet(-1); 
  dpmjet->SetMomentumRange(0, 999999.);
  dpmjet->SetThetaRange(0., 180.);
  dpmjet->SetYRange(-12.,12.);
  dpmjet->SetPtRange(0,1000.);
  dpmjet->SetProcess(kDpmMb);
  dpmjet->SetEnergyCMS(energy);
  dpmjet->SetCrossingAngle(0,0.000280);
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

  // Run number
  if (gSystem->Getenv("DC_RUN")) {
    runNumber = atoi(gSystem->Getenv("DC_RUN"));
  }

  // Physics lists
  if (gSystem->Getenv("CONFIG_PHYSICSLIST")) {
    for (Int_t iList = 0; iList < kListMax; iList++) {
      if (strcmp(gSystem->Getenv("CONFIG_PHYSICSLIST"), physicsListName[iList])==0){
	physicslist = (PhysicsList_t)iList;
	cout<<"Physics list set to "<< physicsListName[iList]<<endl;
      }
    }
  }

}
