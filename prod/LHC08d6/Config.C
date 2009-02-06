// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config.C++")

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
#include "STEER/AliMagF.h"
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
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "ACORDE/AliACORDEv1.h"
#include "VZERO/AliVZEROv7.h"
#endif

//--- Functions ---
void ProcessEnvironmentVars();
Float_t EtaToTheta(Float_t arg);
void    LoadPythia();

//--- Trigger config ---
enum TrigConf_t
{
    kDefaultPPTrig, kDefaultPbPbTrig
};
const char * TrigConfName[] = {
    "p-p","Pb-Pb"
};

// This part for configuration
static AliMagF::BMap_t         mag         = AliMagF::k5kG;
static TrigConf_t    trig        = kDefaultPbPbTrig; // default pp trigger configuration
static Int_t         runNumber   = 0;
static Int_t         eventNumber = 0;

//========================//
// Set Random Number seed //
//========================//
TDatime dt;
static UInt_t seed    = dt.Get();
static Int_t   runNumber= 0;
static Float_t bMin = 0.;
static Float_t bMax = 5.;
static UInt_t  quench = 1;
static UInt_t  shad = 1;
static Float_t etaMin = -8.0;
static Float_t etaMax = 8.0;
static Float_t phiMin = 0.;
static Float_t phiMax = 360.;
static TString comment;

void Config()
{
  // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
  // Theta range given through pseudorapidity limits 22/6/2001

  // Set Random Number seed
  //gRandom->SetSeed(123456); // Set 0 to use the current time

  AliLog::Message(AliLog::kInfo, Form("Seed for random number generation = %d",gRandom->GetSeed()), "Config.C", "Config.C", "Config()","Config.C", __LINE__);

  // Get settings from environment variables
  ProcessEnvironmentVars();

  // Load Pythia libraries
  LoadPythia();
  // Libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("libgeant321");
#endif

  new     TGeant3TGeo("C++ Interface to Geant3");

  // Output every 100 tracks
  ((TGeant3*)gMC)->SetSWIT(4,100);

  AliRunLoader* rl=0x0;

    AliLog::Message(AliLog::kInfo, "Creating Run Loader", "Config.C", "Config.C", "Config()"," Config.C", __LINE__);

    rl = AliRunLoader::Open("galice.root",
			    AliConfig::GetDefaultEventFolderName(),
			    "recreate");
    if (rl == 0x0)
      {
	gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
	return;
      }
    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(2);
    gAlice->SetRunLoader(rl);

    // gAlice->SetGeometryFromFile("geometry.root");

    // Uncomment if you want to load geometry from OCDB!   >>>>
/*
    if(!AliCDBManager::Instance()->IsDefaultStorageSet()){
	 cout << "#####################################################" << endl;
	 cout << "#                                                   #" << endl;
	 cout << "#     WARNING: CDB DEFAULT STORAGE NOT SET !!!      #" << endl;
	 cout << "#     SETTING IT TO local://$ALICE_ROOT !!!         #" << endl;
	 cout << "#                                                   #" << endl;
	 cout << "#####################################################" << endl;

         AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
    }

    if(AliCDBManager::Instance()->GetRun() < 0){
	 cout << "#####################################################" << endl;
	 cout << "#                                                   #" << endl;
	 cout << "#     WARNING: RUN NUMBER NOT SET !!!               #" << endl;
	 cout << "#     SETTING IT TO 0 !!!                           #" << endl;
	 cout << "#                                                   #" << endl;
	 cout << "#####################################################" << endl;

         AliCDBManager::Instance()->SetRun(0);
    }
    gAlice->SetGeometryFromCDB();
*/
    // Uncomment if you want to load geometry from OCDB!   <<<<

    // Set the trigger configuration
    gAlice->SetTriggerDescriptor(TrigConfName[trig]);
    cout<<"Trigger configuration is set to  "<<TrigConfName[trig]<<endl;


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


    AliGenHijing *gener = new AliGenHijing(-1);
    // centre of mass energy
    gener->SetEnergyCMS(5500);
    // reference frame
    gener->SetReferenceFrame("CMS     ");
    // projectile
    gener->SetProjectile("A       ", 208, 82);
    gener->SetTarget    ("A       ", 208, 82);
    // impact parameter range
    gener->SetImpactParameterRange(bMin, bMax); // bMin = 0 - bMax = 3
    // evaluate cross section before run
    gener->SetEvaluate(0);
    // tell hijing to keep the full parent child chain
    gener->KeepFullEvent();
    // enable jet quenching
    gener->SetJetQuenching(quench); // 1
    // enable shadowing
    gener->SetShadowing(shad); // 1
    // neutral pion and heavy particle decays switched off
    gener->SetDecaysOff(1);
    // Don't track spectators
    gener->SetSpectators(0);
    // trigger
    //	gener->SetTrigger(0);
    // kinematic selection
    gener->SetSelectAll(0);
    // momentum range
    gener->SetMomentumRange(0,999);
    // No resytriction on phi, theta
    Float_t thmin = EtaToTheta(etaMax); // Theta min <---> eta max 2.
    Float_t thmax = EtaToTheta(etaMin); // Theta max <---> eta min -2.
    gener->SetPhiRange(phiMin, phiMax); // 0 - 360
    gener->SetThetaRange(thmin,thmax);
    // PRIMARY VERTEX
    gener->SetOrigin(0, 0, 0);  //vertex position
//    gener->SetSigma(0, 0, 5.3);   //Sigma in (X,Y,Z) (cm) on IP position
//    gener->SetCutVertexZ(3.);        // Truncate at 3 sigma
//    gener->SetVertexSmear(kPerEvent);

    // Size of the interaction diamond
    // Longitudinal
    Float_t sigmaz  = 7.55 / TMath::Sqrt(2.); // [cm]
    // Transverse
    Float_t betast  = 10;                 // beta* [m]
    Float_t eps     = 3.75e-6;            // emittance [m]
//    Float_t gamma   = 7000. / 0.938272;   // relativistic gamma [1]
    Float_t gamma   = 2750. / 0.938272;   // relativistic gamma [1]
    Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
    printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);

    gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
    gener->SetCutVertexZ(3.);        // Truncate at 3 sigma
    gener->SetVertexSmear(kPerEvent);

    //
    // Activate this line if you want the vertex smearing to happen
    // track by track
    //
    //gener->SetVertexSmear(perTrack);

    gener->Init();

    // Field (L3 0.5 T)
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 2, 1., 10., mag));

    Int_t   iABSO   = 1;
    Int_t   iDIPO   = 1;
    Int_t   iFMD    = 1;
    Int_t   iFRAME  = 1;
    Int_t   iHALL   = 1;
    Int_t   iITS    = 1;
    Int_t   iMAG    = 1;
    Int_t   iMUON   = 1;
    Int_t   iPHOS   = 1;
    Int_t   iPIPE   = 1;
    Int_t   iPMD    = 1;
    Int_t   iHMPID  = 1;
    Int_t   iSHIL   = 1;
    Int_t   iT0     = 1;
    Int_t   iTOF    = 1;
    Int_t   iTPC    = 1;
    Int_t   iTRD    = 1;
    Int_t   iZDC    = 1;
    Int_t   iEMCAL  = 1;
    Int_t   iACORDE = 1;
    Int_t   iVZERO  = 1;
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

        AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
    }

    if (iITS)
    {
        //=================== ITS parameters ============================

	AliITS *ITS  = new AliITSv11Hybrid("ITS","ITS v11Hybrid");
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
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETE");
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

     AliLog::Message(AliLog::kInfo, "End of Config", "Config.C", "Config.C", "Config()"," Config.C", __LINE__);

}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

void ProcessEnvironmentVars()
{
    cout << "Processing environment variables" << endl;
    // Random Number seed
    if (gSystem->Getenv("CONFIG_SEED")) {
      seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    }

    gRandom->SetSeed(seed);
    cout<<"Seed for random number generation= "<<seed<<endl;

    // Run Number
    if (gSystem->Getenv("DC_RUN")) {
      runNumber = atoi(gSystem->Getenv("DC_RUN"));
    }
    cout<<"Run number "<<runNumber<<endl;

    // Impact param
    if (gSystem->Getenv("CONFIG_BMIN")) {
      bMin = atof(gSystem->Getenv("CONFIG_BMIN"));
    }

    if (gSystem->Getenv("CONFIG_BMAX")) {
      bMax = atof(gSystem->Getenv("CONFIG_BMAX"));
    }
    cout<<"Impact parameter in ["<<bMin<<","<<bMax<<"]"<<endl;

    // Quenching scenario
    if (gSystem->Getenv("QUENCH")) {
      quench = atoi(gSystem->Getenv("QUENCH"));
    }
    if(quench==1) cout<<"With quenching "<<endl;
    else cout<<"Without quenching "<<endl;

    // Quenching scenario
    if (gSystem->Getenv("QHAT")) {
      shad = atoi(gSystem->Getenv("QHAT"));
    }
    if(shad==1) cout<<"With shadowing "<<endl;
    else cout<<"Without shadowing "<<endl;

    // Acceptance param
    if (gSystem->Getenv("CONFIG_ETAMIN")) {
      etaMin = atof(gSystem->Getenv("CONFIG_ETAMIN"));
    }

    if (gSystem->Getenv("CONFIG_ETAMAX")) {
      etaMax = atof(gSystem->Getenv("CONFIG_ETAMAX"));
    }
    cout<<"Eta acceptance ["<<etaMin<<","<<etaMax<<"]"<<endl;

    if (gSystem->Getenv("CONFIG_PHIMIN")) {
      phiMin = atof(gSystem->Getenv("CONFIG_PHIMIN"));
    }

    if (gSystem->Getenv("CONFIG_PHIMAX")) {
      phiMax = atof(gSystem->Getenv("CONFIG_PHIMAX"));
    }
    cout<<"Phi acceptance ["<<phiMin<<","<<phiMax<<"]"<<endl;


}

void LoadPythia()
{
    // Load Pythia related libraries
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
}
