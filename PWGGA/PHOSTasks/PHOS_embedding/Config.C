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
  rl->SetNumberOfEventsPerFile(3000);
  gAlice->SetRunLoader(rl);
  // gAlice->SetGeometryFromFile("geometry.root");
  // gAlice->SetGeometryFromCDB();
  
  // Set the trigger configuration: proton-proton
 // gAlice->SetTriggerDescriptor("p-p");

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
/*
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
*/  
/*
   AliGenCocktail *gener  = new AliGenCocktail();

   AliGenPHOSlib * lib = new AliGenPHOSlib() ;
   //4pi0
   AliGenParam *genPHOS = new AliGenParam(4,AliGenPHOSlib::kPion,lib->GetPt(AliGenPHOSlib::kPi0),lib->GetY(AliGenPHOSlib::kPi0Flat),lib->GetIp(AliGenPHOSlib::kPi0)) ;
   genPHOS->SetPhiRange(250.,330.) ;
   genPHOS->SetYRange(-0.15.,0.15) ;
   genPHOS->SetPtRange(0.5,30.) ;
   gener->AddGenerator(genPHOS,"PHOS",1.) ;
*/

printf("Creating FLAT generator \n") ;

   gener = new AliGenBox(2) ;
   gener->SetPhiRange(250.,330.) ;
   gener->SetYRange(-0.15.,0.15) ;
   gener->SetEtaRange(-0.15,0.15) ;
   gener->SetPtRange(0.5,30.) ;
   gener->SetPart(111) ;

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
  Float_t betast  = 10.;                 // beta* [m]
  if (runNumber >= 117048) betast = 2.;
  printf("beta* for run# %8d is %13.3f", runNumber, betast);
  //
  Float_t eps     = 3.75e-6;            // emittance [m]
  Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]
  Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
  printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);
    
  gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetVertexSmear(kPerEvent);
  gener->Init();

  printf("\n \n Comment: %s \n \n", comment.Data());

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

	AliITS *ITS  = new AliITSv11Hybrid("ITS","ITS v11Hybrid");
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

        AliZDC *ZDC = new AliZDCv3("ZDC", "normal ZDC");
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

        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_FIRSTYEAR");
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
      TH1 *weight = new TH1D("weight","weight",201,-0.5,200.5);
      weight->SetBinContent(1,5.49443);
      weight->SetBinContent(2,8.770816);
      weight->SetBinContent(6,0.4568624);
      weight->SetBinContent(7,0.2919915);
      weight->SetBinContent(8,0.6674189);
      weight->SetBinContent(9,0.364737);
      weight->SetBinContent(10,0.8818444);
      weight->SetBinContent(11,0.531885);
      weight->SetBinContent(12,1.035197);
      weight->SetBinContent(13,0.9394057);
      weight->SetBinContent(14,0.9643193);
      weight->SetBinContent(15,0.94543);
      weight->SetBinContent(16,0.9426507);
      weight->SetBinContent(17,0.9423649);
      weight->SetBinContent(18,0.789456);
      weight->SetBinContent(19,1.149026);
      weight->SetBinContent(20,1.100491);
      weight->SetBinContent(21,0.6350525);
      weight->SetBinContent(22,1.351941);
      weight->SetBinContent(23,0.03233504);
      weight->SetBinContent(24,0.9574557);
      weight->SetBinContent(25,0.868133);
      weight->SetBinContent(26,1.030998);
      weight->SetBinContent(27,1.08897);
      weight->SetBinContent(28,1.251382);
      weight->SetBinContent(29,0.1391099);
      weight->SetBinContent(30,1.192876);
      weight->SetBinContent(31,0.448944);
      weight->SetBinContent(32,1);
      weight->SetBinContent(33,1);
      weight->SetBinContent(34,1);
      weight->SetBinContent(35,1);
      weight->SetBinContent(36,0.9999997);
      weight->SetBinContent(37,0.9999997);
      weight->SetBinContent(38,0.9999996);
      weight->SetBinContent(39,0.9999996);
      weight->SetBinContent(40,0.9999995);
      weight->SetBinContent(41,0.9999993);
      weight->SetBinContent(42,1);
      weight->SetBinContent(43,1);
      weight->SetBinContent(44,1);
      weight->SetBinContent(45,1);
      weight->SetBinContent(46,1);
      weight->SetBinContent(47,0.9999999);
      weight->SetBinContent(48,0.9999998);
      weight->SetBinContent(49,0.9999998);
      weight->SetBinContent(50,0.9999999);
      weight->SetBinContent(51,0.9999999);
      weight->SetBinContent(52,0.9999999);
      weight->SetBinContent(53,0.9999999);
      weight->SetBinContent(54,0.9999998);
      weight->SetBinContent(55,0.9999998);
      weight->SetBinContent(56,0.9999998);
      weight->SetBinContent(57,0.9999997);
      weight->SetBinContent(58,0.9999996);
      weight->SetBinContent(59,0.9999995);
      weight->SetBinContent(60,1);
      weight->SetBinContent(61,1);
      weight->SetBinContent(62,1);
      weight->SetBinContent(63,1);
      weight->SetBinContent(64,1);
      weight->SetBinContent(65,0.9999999);
      weight->SetBinContent(66,0.9999998);
      weight->SetBinContent(67,0.9999998);
      weight->SetBinContent(68,0.9999999);
      weight->SetBinContent(69,1);
      weight->SetBinContent(70,1);
      weight->SetBinContent(71,0.9999997);
      weight->SetBinContent(72,0.9999995);
      weight->SetBinContent(73,0.9999994);
      weight->SetBinContent(74,1);
      weight->SetBinContent(75,1);
      weight->SetBinContent(76,1);
      weight->SetBinContent(77,1);
      weight->SetBinContent(78,0.9999999);
      weight->SetBinContent(79,1);
      weight->SetBinContent(80,1);
      weight->SetEntries(526);
        
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

    // Run number
    if (gSystem->Getenv("DC_RUN")) {
      runNumber = atoi(gSystem->Getenv("DC_RUN"));
    }
}
