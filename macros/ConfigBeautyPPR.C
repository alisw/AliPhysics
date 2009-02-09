#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3.h>
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "PYTHIA6/AliGenPythia.h"
#include "STEER/AliMagFCM.h"
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
#include "FMD/AliFMDv0.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "ACORDE/AliACORDEv1.h"
#endif

void LoadPythia();

void Config()
{
 
  //
  // Set Random Number seed
  TDatime dt;
  UInt_t curtime=dt.Get();
  UInt_t procid=gSystem->GetPid();
  UInt_t seed=curtime-procid;

  //  gRandom->SetSeed(seed);
  gRandom->SetSeed(12345);
  cerr<<"Seed for random number generation= "<<seed<<endl; 
  // Load Pythia libraries
  LoadPythia();
  // libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("libgeant321");
#endif

  new TGeant3TGeo("C++ Interface to Geant3");

  //=======================================================================
  //  Create the output file
   
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()){
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    AliCDBManager::Instance()->SetRun(0);
  }
  
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
  // gAlice->SetGeometryFromFile("geometry.root");
  // gAlice->SetGeometryFromCDB();

  // Set the trigger configuration
  gAlice->SetTriggerDescriptor("Pb-Pb");
  cout<<"Trigger configuration is set to  Pb-Pb"<<endl;

  //
  // Set External decayer
  AliDecayer* decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);


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


  // AliGenPythia *gener = new AliGenPythia(ntracks);
  AliGenPythia *gener = new AliGenPythia(-1);

  gener->SetMomentumRange(0,999);
  gener->SetPhiRange(0.,360.);
  gener->SetThetaRange(0,180);
  gener->SetYRange(-999,999);
  //gener->SetPtRange(0,100);
  gener->SetOrigin(0,0,0);          // vertex position
  //gener->SetVertexSmear(kPerEvent); 
  gener->SetSigma(0,0,5.3);  // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetTrackingFlag(0);
  // gener->SetForceDecay(kHadronicD);

  //
  // The following settings select the Pythia parameters tuned to agree
  // with beauty NLO calculation for Pb-Pb @ 5.5 TeV with MNR code.
  //
  gener->SetProcess(kPyBeautyPbPbMNR);
  gener->SetStrucFunc(kCTEQ4L);
  gener->SetPtHard(2.75,-1.0);
  gener->SetEnergyCMS(5500.);
  gener->SetNuclei(208,208); // Pb-Pb collisions
  // Force no decay heavy quark mesons
  gener->SetForceDecay(kNoDecayHeavy);

  gener->Init();
  // 
  // Activate this line if you want the vertex smearing to happen
  // track by track
  //
  //gener->SetVertexSmear(perTrack); 

  // Field (L3 0.4 T)
  AliMagFCM* field = new AliMagFCM(
	      "Map2","$(ALICE_ROOT)/data/field01.dat", 2, 1., 10.);
  field->SetSolenoidField(4.);
  gAlice->SetField(field);    

  Int_t iABSO=0;
  Int_t iACORDE=0;
  Int_t iDIPO=0;
  Int_t iFMD=0;
  Int_t iFRAME=0;
  Int_t iHALL=0;
  Int_t iITS=0;
  Int_t iMAG=0;
  Int_t iMUON=0;
  Int_t iPHOS=0;
  Int_t iPIPE=0;
  Int_t iPMD=0;
  Int_t iHMPID=0;
  Int_t iSHIL=0;
  Int_t iT0=0;
  Int_t iTOF=0;
  Int_t iTPC=0;
  Int_t iTRD=0;
  Int_t iZDC=0;

  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY","Alice envelop");

  if(iMAG) {
    //=================== MAG parameters ============================
    // --- Start with Magnet since detector layouts may be depending ---
    // --- on the selected Magnet dimensions ---
    AliMAG *MAG  = new AliMAG("MAG","Magnet");
  }


  if(iABSO) {
    //=================== ABSO parameters ============================
    AliABSO *ABSO  = new AliABSOv3("ABSO","Muon Absorber");
  }

  if(iDIPO) {
    //=================== DIPO parameters ============================

    AliDIPO *DIPO  = new AliDIPOv3("DIPO","Dipole version 3");
  }

  if(iHALL) {
    //=================== HALL parameters ============================

    AliHALL *HALL  = new AliHALLv3("HALL","Alice Hall");
  }


  if(iFRAME) {
    //=================== FRAME parameters ============================

    AliFRAME *FRAME  = new AliFRAMEv2("FRAME","Space Frame");
    FRAME->SetHoles(1);
  }

  if(iSHIL) {
    //=================== SHIL parameters ============================

    AliSHIL *SHIL  = new AliSHILv3("SHIL","Shielding");
  }


  if(iPIPE) {
    //=================== PIPE parameters ============================

    AliPIPE *PIPE  = new AliPIPEv3("PIPE","Beam Pipe");
  }


  if(iITS) {
    //=================== ITS parameters ============================

    AliITS *ITS  = new AliITSv11Hybrid("ITS","ITS v11Hybrid");
  }
  

  if(iTPC) {
    //============================ TPC parameters ===================
    AliTPC *TPC  = new AliTPCv2("TPC","Default");
  }


  if(iTOF) {
    //=================== TOF parameters ============================
    AliTOF *TOF  = new AliTOFv6T0("TOF","normal TOF");
  }

  if(iHMPID) {
    //=================== HMPID parameters ===========================
    AliHMPID *HMPID  = new AliHMPIDv3("HMPID","normal HMPID");    

  }


  if(iZDC) {
    //=================== ZDC parameters ============================

    AliZDC *ZDC  = new AliZDCv3("ZDC","normal ZDC");
  }

  if(iACORDE) {
    //=================== ACORDE parameters ============================

    AliACORDE *ACORDE  = new AliACORDEv1("ACORDE","normal ACORDE");
  }

  if(iTRD) {
    //=================== TRD parameters ============================
  
    AliTRD *TRD  = new AliTRDv1("TRD","TRD slow simulator");
  }

  if(iFMD) {
    //=================== FMD parameters ============================

    AliFMD *FMD  = new AliFMDv0("FMD","normal FMD");
  }

  if(iMUON) {
    //=================== MUON parameters ===========================
    AliMUON *MUON  = new AliMUONv1("MUON","default");
  }
 
  //=================== PHOS parameters ===========================

  if(iPHOS) {
    AliPHOS *PHOS  = new AliPHOSv1("PHOS","GPS2");
  }


  if(iPMD) {
    //=================== PMD parameters ============================

    AliPMD *PMD  = new AliPMDv1("PMD","normal PMD");
    PMD->SetPAR(1., 1., 0.8, 0.02);
    PMD->SetIN(6., 18., -580., 27., 27.);
    PMD->SetGEO(0.0, 0.2, 4.);
    PMD->SetPadSize(0.8, 1.0, 1.0, 1.5);

  }

  if(iT0) {
    //=================== T0 parameters ============================
    AliT0 *T0  = new AliT0v1("T0","T0 Detector");
  }

         
}


void LoadPythia()
{
    // Load Pythia related libraries
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
}
