//
// Configuration for the Physics Data Challenge 2006
//
// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config_PDC06.C++")
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenParam.h"
#include "EVGEN/AliGenMUONlib.h"
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "PYTHIA6/AliGenPythia.h"
#include "STEER/AliMagF.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"
#include "STRUCT/AliABSOv0.h"
#include "STRUCT/AliDIPOv2.h"
#include "STRUCT/AliHALL.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv2.h"
#include "STRUCT/AliPIPEv0.h"
#include "ITS/AliITSgeom.h"
#include "ITS/AliITSvPPRasymmFMD.h"
#include "TPC/AliTPCv2.h"
#include "TOF/AliTOFv5T0.h"
#include "HMPID/AliHMPIDv1.h"
#include "ZDC/AliZDCv2.h"
#include "TRD/AliTRDv1.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "CRT/AliCRTv0.h"
#include "VZERO/AliVZEROv7.h"
#endif
//--- Trigger config ---
enum TrigConf_t
{
    kDefaultPPTrig, kDefaultPbPbTrig
};
const char * TrigConfName[] = {
    "p-p","Pb-Pb"
};
//--- Decay Mode ---
enum DecayHvFl_t 
{
  kNature,  kHadr, kSemiEl, kSemiMu
};
//--- Rapidity Cut ---
enum YCut_t
{
  kFull, kBarrel, kMuonArm
};
void ProcessEnvironmentVars();
// This part for configuration
static DecayHvFl_t   decHvFl  = kNature; 
static YCut_t        ycut     = kFull;
static AliMagF::BMap_t mag    = AliMagF::k5kG; 
static TrigConf_t    trig     = kDefaultPPTrig; 
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
  //  gRandom->SetSeed(12345);
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
  rl->SetNumberOfEventsPerFile(1000);
  gAlice->SetRunLoader(rl);
  
  // Set the trigger configuration
  gAlice->SetTriggerDescriptor(TrigConfName[trig]);
  cout<<"Trigger configuration is set to  "<<TrigConfName[trig]<<endl;
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
    //    ((TGeant3 *) gMC)->SetSWIT(2,2); 
    //    ((TGeant3 *) gMC)->SetDEBU(1,999,1); 
  // Set External decayer //
  //======================//
  TVirtualMCDecayer* decayer = new AliDecayerPythia();
  // DECAYS
  //
  switch(decHvFl) {
  case kNature:
    decayer->SetForceDecay(kAll);
    break;
  case kHadr:
    decayer->SetForceDecay(kHadronicD);
    break;
  case kSemiEl:
    decayer->SetForceDecay(kSemiElectronic);
    break;
  case kSemiMu:
    decayer->SetForceDecay(kSemiMuonic);
    break;
  }
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
  //=========================//
  // Generator Configuration //
  //=========================//
  AliGenBox* gener = new AliGenBox();
  gener->SetThetaRange(45,135);
  gener->SetPhiRange(30,150);
  gener->SetMomentumRange(9.8,10.2);
  gener->SetPart(kMuonMinus);
  gener->SetNumberParticles(20);
  gener->Init();
  // FIELD
  //    
  if (mag == AliMagF::k2kG) {
    comment = comment.Append(" | L3 field 0.2 T");
  } else if (mag == AliMagF::k5kG) {
    comment = comment.Append(" | L3 field 0.5 T");
  }
  printf("\n \n Comment: %s \n \n", comment.Data());
    
  AliMagF* field = new AliMagF("Maps","Maps",2, -1.,-1., 10.,mag);
  TGeoGlobalMagField::Instance()->SetField(field);

  rl->CdGAFile();

  Int_t iABSO  = 0;
  Int_t iACORDE = 1;
  Int_t iDIPO  = 0;
  Int_t iEMCAL = 0;
  Int_t iFMD   = 0;
  Int_t iFRAME = 0;
  Int_t iHALL  = 0;
  Int_t iITS   = 0;
  Int_t iMAG   = 0;
  Int_t iMUON  = 0;
  Int_t iPHOS  = 0;
  Int_t iPIPE  = 0;
  Int_t iPMD   = 0;
  Int_t iHMPID  = 0;
  Int_t iSHIL  = 0;
  Int_t iT0 = 0;
  Int_t iTOF   = 0;
  Int_t iTPC   = 0;
  Int_t iTRD   = 0;
  Int_t iVZERO = 0;
  Int_t iZDC   = 0;
  
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
        AliABSO *ABSO = new AliABSOv0("ABSO", "Muon Absorber");
    }
    if (iDIPO)
    {
        //=================== DIPO parameters ============================
        AliDIPO *DIPO = new AliDIPOv2("DIPO", "Dipole version 2");
    }
    if (iHALL)
    {
        //=================== HALL parameters ============================
        AliHALL *HALL = new AliHALL("HALL", "Alice Hall");
    }
    if (iFRAME)
    {
        //=================== FRAME parameters ============================
        AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
    }
    if (iSHIL)
    {
        //=================== SHIL parameters ============================
        AliSHIL *SHIL = new AliSHILv2("SHIL", "Shielding Version 2");
    }
    if (iPIPE)
    {
        //=================== PIPE parameters ============================
        AliPIPE *PIPE = new AliPIPEv0("PIPE", "Beam Pipe");
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
	AliITSvPPRasymmFMD *ITS  = new AliITSvPPRasymmFMD("ITS","New ITS PPR detailed version with asymmetric services");
	ITS->SetMinorVersion(2);  // don't touch this parameter if you're not an ITS developer
	ITS->SetReadDet(kFALSE);	  // don't touch this parameter if you're not an ITS developer
    //    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
	ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
	ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
	ITS->SetThicknessChip1(150.);  // chip thickness on layer 1 must be in the range [150,300]
	ITS->SetThicknessChip2(150.);  // chip thickness on layer 2 must be in the range [150,300]
	ITS->SetRails(0);	     // 1 --> rails in ; 0 --> rails out
	ITS->SetCoolingFluid(1);   // 1 --> water ; 0 --> freon
    // Coarse geometries (warning: no hits are produced with these coarse geometries and they unuseful 
    // for reconstruction !):
    //                                                     
    //
    //AliITSvPPRcoarseasymm *ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS PPR coarse version with asymmetric services");
    //ITS->SetRails(0);                // 1 --> rails in ; 0 --> rails out
    //ITS->SetSupportMaterial(0);      // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    //
    //AliITS *ITS  = new AliITSvPPRcoarsesymm("ITS","New ITS PPR coarse version with symmetric services");
    //ITS->SetRails(0);                // 1 --> rails in ; 0 --> rails out
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
    if (iTPC)
    {
      //============================ TPC parameters =====================
        AliTPC *TPC = new AliTPCv2("TPC", "Default");
    }
    if (iTOF) {
        //=================== TOF parameters ============================
	AliTOF *TOF = new AliTOFv5T0("TOF", "normal TOF");
	// Partial geometry: modules at 2,3,4,6,7,11,12,14,15,16
	// starting at 6h in positive direction
	//	Int_t TOFSectors[18]={-1,-1,0,0,0,-1,0,0,-1,-1,-1,0,0,-1,0,0,0,0};
	// Partial geometry: modules at 1,2,6,7,9,10,11,12,15,16,17
	// (ALICE numbering convention)
       	Int_t TOFSectors[18]={-1,0,0,-1,-1,-1,0,0,-1,0,0,0,0,-1,-1,0,0,0};
	TOF->SetTOFSectors(TOFSectors);
    }
    if (iHMPID)
    {
        //=================== HMPID parameters ===========================
        AliHMPID *HMPID = new AliHMPIDv1("HMPID", "normal HMPID");
    }
    if (iZDC)
    {
        //=================== ZDC parameters ============================
        AliZDC *ZDC = new AliZDCv2("ZDC", "normal ZDC");
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
	// ACORDE->SetITSGeometry(kTRUE);
	// ACORDE->SetCreateCavern(kFALSE);
    }
     if (iVZERO)
    {
        //=================== VZERO parameters =============================
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
    // Random Number seed
    if (gSystem->Getenv("CONFIG_SEED")) {
      seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    }
}
