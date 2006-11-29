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
#include "STRUCT/AliABSOv0.h"
#include "STRUCT/AliDIPOv2.h"
#include "STRUCT/AliHALL.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv2.h"
#include "STRUCT/AliPIPEv0.h"
#include "ITS/AliITSvPPRasymm.h"
#include "TPC/AliTPCv2.h"
#include "TOF/AliTOFv5T0.h"
#include "HMPID/AliHMPIDv1.h"
#include "ZDC/AliZDCv1.h"
#include "TRD/AliTRDv1.h"
#include "FMD/AliFMDv0.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "ACORDE/AliACORDEv1.h"
#endif

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

  // libraries required by geant321
#if defined(__CINT__)
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
  rl->SetNumberOfEventsPerFile(3);
  gAlice->SetRunLoader(rl);

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
    AliABSO *ABSO  = new AliABSOv0("ABSO","Muon Absorber");
  }

  if(iDIPO) {
    //=================== DIPO parameters ============================

    AliDIPO *DIPO  = new AliDIPOv2("DIPO","Dipole version 2");
  }

  if(iHALL) {
    //=================== HALL parameters ============================

    AliHALL *HALL  = new AliHALL("HALL","Alice Hall");
  }


  if(iFRAME) {
    //=================== FRAME parameters ============================

    AliFRAME *FRAME  = new AliFRAMEv2("FRAME","Space Frame");

  }

  if(iSHIL) {
    //=================== SHIL parameters ============================

    AliSHIL *SHIL  = new AliSHILv2("SHIL","Shielding");
  }


  if(iPIPE) {
    //=================== PIPE parameters ============================

    AliPIPE *PIPE  = new AliPIPEv0("PIPE","Beam Pipe");
  }


  if(iITS) {

//=================== ITS parameters ============================
    //
    // As the innermost detector in ALICE, the Inner Tracking System "impacts" on
    // almost all other detectors. This involves the fact that the ITS geometry
    // still has several options to be followed in parallel in order to determine
    // the best set-up which minimizes the induced background. All the geometries
    // available to date are described in the following. Read carefully the comments
    // and use the default version (the only one uncommented) unless you are making
    // comparisons and you know what you are doing. In this case just uncomment the
    // ITS geometry you want to use and run Aliroot.
    //
    // Detailed geometries:         
    //
    //
    //AliITS *ITS  = new AliITSv5symm("ITS","Updated ITS TDR detailed version with symmetric services");
    //
    //AliITS *ITS  = new AliITSv5asymm("ITS","Updates ITS TDR detailed version with asymmetric services");
    //
    AliITSvPPRasymm *ITS  = new AliITSvPPRasymm("ITS","New ITS PPR detailed version with asymmetric services");
    ITS->SetMinorVersion(2);                                         // don't touch this parameter if you're not an ITS developer
    ITS->SetReadDet(kFALSE);                                         // don't touch this parameter if you're not an ITS developer
    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
    ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [150,300]
    ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [150,300]
    ITS->SetThicknessChip1(150.);  // chip thickness on layer 1 must be in the range [100,300]
    ITS->SetThicknessChip2(150.);  // chip thickness on layer 2 must be in the range [100,300]
    ITS->SetRails(1);	           // 1 --> rails in ; 0 --> rails out
    ITS->SetCoolingFluid(1);       // 1 --> water ; 0 --> freon
    //
    // Coarse geometries (warning: no hits are produced with these coarse geometries and they unuseful 
    // for reconstruction !):
    //                                                     
    //
    //AliITSvPPRcoarseasymm *ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS PPR coarse version with asymmetric services");
    //ITS->SetRails(1);                // 1 --> rails in ; 0 --> rails out
    //ITS->SetSupportMaterial(0);      // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    //
    //AliITS *ITS  = new AliITSvPPRcoarsesymm("ITS","New ITS PPR coarse version with symmetric services");
    //ITS->SetRails(1);                // 1 --> rails in ; 0 --> rails out
    //ITS->SetSupportMaterial(0);      // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    //                      
    //
    //
    // Geant3 <-> EUCLID conversion
    // ============================
    //
    // SetEUCLID is a flag to output (=1) or not to output (=0) both geometry and
    // media to two ASCII files (called by default ITSgeometry.euc and
    // ITSgeometry.tme) in a format understandable to the CAD system EUCLID.
    // The default (=0) means that you dont want to use this facility.
    //
    ITS->SetEUCLID(0);  
  }
  

  if(iTPC) {
    //============================ TPC parameters ===================
    AliTPC *TPC  = new AliTPCv2("TPC","Default");
  }


  if(iTOF) {
    //=================== TOF parameters ============================
    AliTOF *TOF  = new AliTOFv5T0("TOF","normal TOF");
  }

  if(iHMPID) {
    //=================== HMPID parameters ===========================
    AliHMPID *HMPID  = new AliHMPIDv1("HMPID","normal HMPID");    

  }


  if(iZDC) {
    //=================== ZDC parameters ============================

    AliZDC *ZDC  = new AliZDCv1("ZDC","normal ZDC");
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
