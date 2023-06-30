#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ALICE_ROOT)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
R__ADD_INCLUDE_PATH(./)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliConfig.h"
#include "Pythia.h"
#include "AliGenLib.h"
#include "AliGenCocktail.h"
#include "AliGenParam.h"
#include "AliGenPythia.h"
#include "AliDecayerPythia.h"
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
#include "AliZDCv3.h"
#include "AliTRDv1.h"
#include "AliTRDgeometry.h"
#include "AliFMDv1.h"
#include "AliMUONv1.h"
#include "AliPHOSv1.h"
#include "AliPHOSSimParam.h"
#include "AliPMDv1.h"
#include "AliT0v1.h"
#include "AliEMCALv2.h"
#include "AliACORDEv1.h"
#include "AliVZEROv7.h"
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

void ProcessEnvironmentVars();

class MyPHOSlib :
public AliGenLib
{
 public:
    GenFunc   GetPt(Int_t param, const char* tname=0) const { GenFunc func = PtGamma; return func;}
    GenFunc   GetY (Int_t param, const char* tname=0) const { GenFunc func = YGamma; return func;}
    GenFunc   GetV2 (Int_t param, const char* tname=0) const { GenFunc func = V2Gamma; return func;}
    GenFuncIp GetIp(Int_t param, const char* tname=0) const {GenFuncIp func = IpGamma; return func;}   
 private:
   static Double_t YGamma( const Double_t *, const Double_t *) { return 1;}
   static Double_t V2Gamma( const Double_t *, const Double_t *) { return 0.0001;}
   static Double_t PtGamma( const Double_t *px, const Double_t *) { return exp(-3.*(*px));}
   static Int_t IpGamma(TRandom */*ran*/) { return 22;}
};

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
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  // if (proc == kPythia6 || proc == kPhojet) {
    gSystem->Load("libpythia6");        // Pythia 6.2
  // } else {
    // gSystem->Load("libpythia6.4.21");   // Pythia 6.4
  // }

printf("========AAAAAA========\n") ;
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libgeant321");

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

   printf("Creating custom generator \n") ;
   AliGenCocktail* gener = new AliGenCocktail();
   MyPHOSlib *mpl = new MyPHOSlib();
   //double yPHOSmin[3]={-0.125,-0.042, 0.042} ;
   //double yPHOSmax[3]={-0.042, 0.042, 0.125} ;
   //double phiPHOSmin[14]={250.,255.,260.,265.,270.,275.,280.,285.,290.,295.,300.,305.,310.,315.} ;
   //double phiPHOSmax[14]={255.,260.,265.,270.,275.,280.,285.,290.,295.,300.,305.,310.,315.,320.} ;
   //double phiPHOSmin[4]={240.,260.,280.,300.} ;
   //double phiPHOSmax[4]={260.,280.,300.,320.} ;

        AliGenParam *genPHOS = new AliGenParam(35,1,mpl->GetPt(1),mpl->GetY(1),mpl->GetV2(1),mpl->GetIp(1));
        genPHOS->SetPhiRange(240.,320.);
        genPHOS->SetYRange(-0.125,0.125) ;
        genPHOS->SetPtRange(0.1,3.) ; 
        genPHOS->SetForceDecay(kNoDecay);
        gener->AddGenerator(genPHOS,"PHOS",1.) ;

  // Longitudinal
  Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
  
  // Transverse
  // beta*
  Float_t betast  = 2.;                 // beta* [m]
  Float_t eps     = 3.75e-6;            // emittance [m]
  Float_t gamma   = 5.02 / 2.0 / 0.938272;  // relativistic gamma [1]
  Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
  printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);
    
  gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetVertexSmear(kPerEvent);
  gener->Init();

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
       ITS->DisableStepManager();
    }

    if (iTPC)
    {
      //============================ TPC parameters =====================

        AliTPC *TPC = new AliTPCv2("TPC", "Default");
        TPC->DisableStepManager();
    }


    if (iTOF) {
        //=================== TOF parameters ============================
  	AliTOF *TOF = new AliTOFv6T0("TOF", "normal TOF");
        TOF->DisableStepManager();
    }


    if (iHMPID)
    {
        //=================== HMPID parameters ===========================

        AliHMPID *HMPID = new AliHMPIDv3("HMPID", "normal HMPID");
        HMPID->DisableStepManager();

    }


    if (iZDC)
    {
        //=================== ZDC parameters ============================

        AliZDC *ZDC = new AliZDCv3("ZDC", "normal ZDC");
        ZDC->DisableStepManager();
    }

    if (iTRD)
    {
        //=================== TRD parameters ============================

        AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
        TRD->DisableStepManager();
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
        FMD->DisableStepManager();
   }

    if (iMUON)
    {
        //=================== MUON parameters ===========================
        // New MUONv1 version (geometry defined via builders)

        AliMUON *MUON = new AliMUONv1("MUON", "default");
        MUON->DisableStepManager();
        
        // activate trigger efficiency by cells

        MUON->SetTriggerEffCells(1);     

     }

    if (iPHOS)
    {
        //=================== PHOS parameters ===========================

     //AliPHOS *PHOS = new AliPHOSv1("PHOS", "noCPV_Modules123");
     AliPHOS *PHOS = new AliPHOSv1("PHOS", "Run2");

    }


    if (iPMD)
    {
        //=================== PMD parameters ============================

        AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");
        PMD->DisableStepManager();
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
        EMCAL->DisableStepManager();
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
