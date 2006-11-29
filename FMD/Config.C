//____________________________________________________________________
//
// $Id$
//
// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"ConfigPPR.C++")
//
/** @file    Config.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:50:29 2006
    @brief   Simulation configuration script
*/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3.h>
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "STEER/AliGenerator.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "EVGEN/AliGenHIJINGpara.h"
#include "THijing/AliGenHijing.h"
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenSlowNucleons.h"
#include "EVGEN/AliSlowNucleonModelExp.h"
#include "EVGEN/AliGenParam.h"
#include "EVGEN/AliGenMUONlib.h"
#include "EVGEN/AliGenMUONCocktail.h"
#include "PYTHIA6/AliGenPythia.h"
#include "STEER/AliMagFMaps.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"
#include "STRUCT/AliABSOv0.h"
#include "STRUCT/AliDIPOv2.h"
#include "STRUCT/AliHALL.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv2.h"
#include "STRUCT/AliPIPEv0.h"
#include "ITS/AliITSvPPRasymmFMD.h"
#include "TPC/AliTPCv2.h"
#include "TOF/AliTOFv4T0.h"
#include "HMPID/AliHMPIDv1.h"
#include "ZDC/AliZDCv2.h"
#include "TRD/AliTRDv1.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "MUON/AliMUONSt1GeometryBuilderV2.h"
#include "MUON/AliMUONSt2GeometryBuilder.h"
#include "MUON/AliMUONSlatGeometryBuilder.h"
#include "MUON/AliMUONTriggerGeometryBuilder.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "START/AliSTARTv1.h"
#include "EMCAL/AliEMCALv1.h"
#include "CRT/AliCRTv0.h"
#include "VZERO/AliVZEROv2.h"
#endif

//____________________________________________________________________
// 
// Generator types 
//
enum EG_t {
  test50,
  kParam_8000,			//
  kParam_4000,			//
  kParam_2000,			//
  kParam_fmd,			//
  kHijing_cent1,		//
  kHijing_cent2,		//
  kHijing_per1,			//
  kHijing_per2,			//
  kHijing_per3,			//
  kHijing_per4,			//
  kHijing_per5,			//
  kHijing_jj25,			//
  kHijing_jj50,			//
  kHijing_jj75,			//
  kHijing_jj100,		//
  kHijing_jj200,		//
  kHijing_gj25,			//
  kHijing_gj50,			//
  kHijing_gj75,			//
  kHijing_gj100,		//
  kHijing_gj200,		//
  kHijing_pA,			//
  kPythia6,			//
  kPythia6Jets20_24,		//
  kPythia6Jets24_29,		//
  kPythia6Jets29_35,		//
  kPythia6Jets35_42,		//
  kPythia6Jets42_50,		//
  kPythia6Jets50_60,		//
  kPythia6Jets60_72,		//
  kPythia6Jets72_86,		//
  kPythia6Jets86_104,		//
  kPythia6Jets104_125,		//
  kPythia6Jets125_150,		//
  kPythia6Jets150_180,		//
  kD0PbPb5500,			//
  kCharmSemiElPbPb5500,		//
  kBeautySemiElPbPb5500,	//
  kCocktailTRD,			//
  kPyJJ,			//
  kPyGJ,			//
  kMuonCocktailCent1,		//
  kMuonCocktailPer1,		//
  kMuonCocktailPer4,		//
  kMuonCocktailCent1HighPt,	//
  kMuonCocktailPer1HighPt,	//
  kMuonCocktailPer4HighPt,	//
  kMuonCocktailCent1Single,	//
  kMuonCocktailPer1Single,	//
  kMuonCocktailPer4Single,
  kFMD1Flat, 
  kFMD2Flat, 
  kFMD3Flat,
  kFMDFlat,
  kEgMax
};

//____________________________________________________________________
// 
// Generator types names
//
const char* egName[kEgMax] = {
  "test50",
  "kParam_8000",		//
  "kParam_4000",		//
  "kParam_2000",		//
  "kParam_fmd",			//
  "kHijing_cent1",		//
  "kHijing_cent2",		//
  "kHijing_per1",		//
  "kHijing_per2",		//
  "kHijing_per3",		//
  "kHijing_per4",		//
  "kHijing_per5",		//
  "kHijing_jj25",		//
  "kHijing_jj50",		//
  "kHijing_jj75",		//
  "kHijing_jj100",		//
  "kHijing_jj200",		//
  "kHijing_gj25",		//
  "kHijing_gj50",		//
  "kHijing_gj75",		//
  "kHijing_gj100",		//
  "kHijing_gj200",		//
  "kHijing_pA",			//
  "kPythia6",			//
  "kPythia6Jets20_24",		//
  "kPythia6Jets24_29",		//
  "kPythia6Jets29_35",		//
  "kPythia6Jets35_42",		//
  "kPythia6Jets42_50",		//
  "kPythia6Jets50_60",		//
  "kPythia6Jets60_72",		//
  "kPythia6Jets72_86",		//
  "kPythia6Jets86_104",		//
  "kPythia6Jets104_125",	//
  "kPythia6Jets125_150",	//
  "kPythia6Jets150_180",	//
  "kD0PbPb5500",		//
  "kCharmSemiElPbPb5500",	//
  "kBeautySemiElPbPb5500",	//
  "kCocktailTRD",		//
  "kPyJJ",			//
  "kPyGJ",			//
  "kMuonCocktailCent1",		//
  "kMuonCocktailPer1",		//
  "kMuonCocktailPer4",		//
  "kMuonCocktailCent1HighPt",	//
  "kMuonCocktailPer1HighPt",	//
  "kMuonCocktailPer4HighPt",	//
  "kMuonCocktailCent1Single",	//
  "kMuonCocktailPer1Single",	//
  "kMuonCocktailPer4Single",
  "kFMD1Flat",
  "kFMD2Flat",
  "kFMD3Flat",
  "kFMDFlat"
};

//____________________________________________________________________
enum Geo_t {
  kHoles,			//
  kNoHoles			//
};

//____________________________________________________________________
enum Rad_t {
  kGluonRadiation,		//
  kNoGluonRadiation		//
};

//____________________________________________________________________
enum Mag_t {
  k2kG,				//
  k4kG,				//
  k5kG				//
};

//____________________________________________________________________
enum MC_t {
  kFLUKA, 
  kGEANT3, 
  kGEANT4, 
  kGEANT3TGEO,
};

//____________________________________________________________________
// Functions
Float_t       EtaToTheta(Float_t eta);
Eg_t          LookupEG(const Char_t* name);
AliGenerator* GeneratorFactory(EG_t eg, Rad_t rad, TString& comment);
AliGenHijing* HijingStandard();
void          ProcessEnvironmentVars(EG_t& eg, Int_t& seed);

//____________________________________________________________________
void 
Config()
{
  //____________________________________________________________________
  // This part for configuration    
  // EG_t  eg   = test50;
  // EG_t  eg   = kParam_fmd;
  // EG_t  eg   = kParam_2000; // kPythia;
  EG_t  eg   = kFMDFlat;
  Geo_t geo  = kNoHoles;
  Rad_t rad  = kGluonRadiation;
  Mag_t mag  = k5kG;
  Int_t seed = 12345; //Set 0 to use the current time
  MC_t  mc   = kGEANT3TGEO;
  
  //____________________________________________________________________
  // Comment line 
  static TString  comment;
  
  //____________________________________________________________________
  // Get settings from environment variables
  ProcessEnvironmentVars(eg, seed);

  //____________________________________________________________________
  // Set Random Number seed
  gRandom->SetSeed(seed);
  cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<endl; 


  //__________________________________________________________________
  switch (mc) {
  case kFLUKA: 
    // 
    // libraries required by fluka21
    // 
    gSystem->Load("libGeom");
    cout << "\t* Loading TFluka..." << endl;  
    gSystem->Load("libTFluka");    
    gSystem->MakeDirectory("peg");
    // 
    // FLUKA MC
    //
    cout << "\t* Instantiating TFluka..." << endl;
    new TFluka("C++ Interface to Fluka", 0/*verbosity*/);
    break;
  case kGEANT3: 
    {
      //
      // Libraries needed by GEANT 3.21 
      //
      gSystem->Load("libgeant321");
      
      // 
      // GEANT 3.21 MC 
      // 
      TGeant3* gmc = new TGeant3("C++ Interface to Geant3");
      gmc->SetSWIT(4, 1000);
    }
    break;
  case kGEANT3TGEO:
    {
      //
      // Libraries needed by GEANT 3.21 
      //
      gSystem->Load("libgeant321");
    
      // 
      // GEANT 3.21 MC 
      // 
      TGeant3TGeo* gmc  = new TGeant3TGeo("C++ Interface to Geant3");
      gmc->SetSWIT(4, 1000);
      Printf("Making a TGeant3TGeo objet");
    }
    break;
  default:
    gAlice->Fatal("Config.C", "No MC type chosen");
    return;
  }

  //__________________________________________________________________
  AliRunLoader* rl = 0;

  cout<<"Config.C: Creating Run Loader ..."<<endl;
  rl = AliRunLoader::Open("galice.root",
			  AliConfig::GetDefaultEventFolderName(),
			  "recreate");
  if (!rl) {
    gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(3);
  gAlice->SetRunLoader(rl);

  //__________________________________________________________________
  // For FLUKA 
  switch (mc) {
  case kFLUKA: 
    {
      //
      // Use kTRUE as argument to generate alice.pemf first
      //
      TString alice_pemf(gSystem->Which(".", "peg/mat17.pemf"));
      if (!alice_pemf.IsNull()) 
	((TFluka*)gMC)->SetGeneratePemf(kFALSE);
      else
	((TFluka*)gMC)->SetGeneratePemf(kTRUE);
      TString flupro(gSystem->Getenv("FLUPRO"));
      if (flupro.IsNull()) 
	Fatal("Config.C", "Environment variable FLUPRO not set");
#if 0
      char* files[] = { "brems_fin.bin", "cohff.bin", "elasct.bin", 
			"gxsect.bin", "nuclear.bin", "sigmapi.bin", 
			0 };
      char* file = files[0];
      while (file) {
	TString which(gSystem->Which(".", file));
	if (which.IsNull()) {
	  if (gSystem->Symlink(Form("%s/%s", flupro.Data(), file), file)!=0) 
	    Fatal("Config.C", "Couldn't link $(FLUPRO)/%s -> .", file);
	}
	file++;
      }
#endif
      TString neuxsc(gSystem->Which(".", "neuxsc.bin"));
      if (neuxsc.IsNull()) 
	gSystem->Symlink(Form("%s/neuxsc_72.bin", flupro.Data()), 
			 "neuxsc.bin"); 
      gSystem->CopyFile("$(FLUPRO)/random.dat", "old.seed", kTRUE);
    }
    break;
  }

  //__________________________________________________________________
  //
  // Set External decayer
#if 0
  AliDecayer *decayer = new AliDecayerPythia();
  switch (eg) {
  case kD0PbPb5500:           decayer->SetForceDecay(kHadronicD);      break;
  case kCharmSemiElPbPb5500:  decayer->SetForceDecay(kSemiElectronic); break;
  case kBeautySemiElPbPb5500: decayer->SetForceDecay(kSemiElectronic); break;
  default:                    decayer->SetForceDecay(kAll);            break;
  }
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
#endif

  //__________________________________________________________________
  // *********** STEERING parameters FOR ALICE SIMULATION ************
  // - Specify event type to be tracked through the ALICE setup
  // - All positions are in cm, angles in degrees, and P and E in GeV 
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

  
  //__________________________________________________________________
  // Generator Configuration
  AliGenerator* gener = GeneratorFactory(eg, rad, comment);
  gener->SetOrigin(0, 0, 0);    // vertex position
  gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetCutVertexZ(1.);     // Truncate at 1 sigma
  gener->SetVertexSmear(kPerEvent); 
  gener->SetTrackingFlag(1);
  gener->Init();
    
  //__________________________________________________________________
  // 
  // Comments 
  // 
  switch (mag) {
  case k2kG: comment = comment.Append(" | L3 field 0.2 T"); break;
  case k4kG: comment = comment.Append(" | L3 field 0.4 T"); break;
  case k5kG: comment = comment.Append(" | L3 field 0.5 T"); break;
  }

  switch (rad) {
  case kGluonRadiation: 
    comment = comment.Append(" | Gluon Radiation On");  break;
  default:
    comment = comment.Append(" | Gluon Radiation Off"); break;
  }

  switch(geo) {
  case kHoles: comment = comment.Append(" | Holes for PHOS/HMPID"); break;
  default:     comment = comment.Append(" | No holes for PHOS/HMPID"); break;
  }

  std::cout << "\n\n Comment: " << comment << "\n" << std::endl;

  //__________________________________________________________________
  // Field (L3 0.4 T)
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., mag);
  field->SetL3ConstField(0); //Using const. field in the barrel
  rl->CdGAFile();
  gAlice->SetField(field);    

  //__________________________________________________________________
  // 
  // Used detectors 
  // 
  Bool_t useABSO  = kFALSE; 
  Bool_t useCRT   = kFALSE; 
  Bool_t useDIPO  = kFALSE; 
  Bool_t useFMD   = kTRUE; 
  Bool_t useFRAME = kFALSE; 
  Bool_t useHALL  = kFALSE; 
  Bool_t useITS   = kFALSE;
  Bool_t useMAG   = kFALSE; 
  Bool_t useMUON  = kFALSE; 
  Bool_t usePHOS  = kFALSE; 
  Bool_t usePIPE  = kFALSE; 
  Bool_t usePMD   = kFALSE; 
  Bool_t useHMPID  = kFALSE; 
  Bool_t useSHIL  = kFALSE; 
  Bool_t useSTART = kFALSE; 
  Bool_t useTOF   = kFALSE; 
  Bool_t useTPC   = kFALSE;
  Bool_t useTRD   = kFALSE; 
  Bool_t useZDC   = kFALSE; 
  Bool_t useEMCAL = kFALSE; 
  Bool_t useVZERO = kFALSE;

  cout << "\t* Creating the detectors ..." << endl;
  // ================= Alice BODY parameters =========================
  AliBODY *BODY = new AliBODY("BODY", "Alice envelop");
  
  
  if (useMAG) {
    // =================== MAG parameters ============================
    // Start with Magnet since detector layouts may be depending on
    // the selected Magnet dimensions 
    AliMAG *MAG = new AliMAG("MAG", "Magnet");
  }

  if (useABSO) {
    // =================== ABSO parameters ===========================
    AliABSO *ABSO = new AliABSOv0("ABSO", "Muon Absorber");
  }

  if (useDIPO) {
    // =================== DIPO parameters ===========================
    AliDIPO *DIPO = new AliDIPOv2("DIPO", "Dipole version 2");
  }

  if (useHALL) {
    // =================== HALL parameters ===========================
    AliHALL *HALL = new AliHALL("HALL", "Alice Hall");
  }


  if (useFRAME) {
    // ================== FRAME parameters ===========================
    AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
    switch (geo) {
    case kHoles: FRAME->SetHoles(1); break;
    default:     FRAME->SetHoles(0); break;
    }
  }

  if (useSHIL) {
    // ================== SHIL parameters ============================
    AliSHIL *SHIL = new AliSHILv2("SHIL", "Shielding Version 2");
  }


  if (usePIPE) {
    // ================== PIPE parameters ============================
    AliPIPE *PIPE = new AliPIPEv0("PIPE", "Beam Pipe");
  }
  
  if (useITS) {
    // =================== ITS parameters ============================
    //
    // As the innermost detector in ALICE, the Inner Tracking System
    // "impacts" on almost all other detectors. This involves the fact
    // that the ITS geometry still has several options to be followed
    // in parallel in order to determine the best set-up which
    // minimizes the induced background. All the geometries available
    // to date are described in the following. Read carefully the
    // comments and use the default version (the only one uncommented)
    // unless you are making comparisons and you know what you are
    // doing. In this case just uncomment the ITS geometry you want to
    // use and run Aliroot.
    //
    // Detailed geometries:
    //
    //
    // AliITS *ITS = 
    //   new AliITSv5symm("ITS", "Updated ITS TDR detailed version "
    //  		  "with symmetric services");
    // AliITS *ITS  = 
    //   new AliITSv5asymm("ITS","Updates ITS TDR detailed version "
    // 			   "with asymmetric services");
    //
    AliITSvPPRasymmFMD *ITS  = 
      new AliITSvPPRasymmFMD("ITS","New ITS PPR detailed version "
			     "with asymmetric services");
     // don't touch this parameter if you're not an ITS developer
    ITS->SetMinorVersion(2); 
    // don't touch this parameter if you're not an ITS developer
    ITS->SetReadDet(kTRUE);
    // don't touch this parameter if you're not an ITS developer
    // ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  
    // detector thickness on layer 1 must be in the range [100,300]
    ITS->SetThicknessDet1(200.);   
    // detector thickness on layer 2 must be in the range [100,300]
    ITS->SetThicknessDet2(200.);   
    // chip thickness on layer 1 must be in the range [150,300]
    ITS->SetThicknessChip1(200.);  
    // chip thickness on layer 2 must be in the range [150,300]
    ITS->SetThicknessChip2(200.);
    // 1 --> rails in ; 0 --> rails out
    ITS->SetRails(0);          
    // 1 --> water ; 0 --> freon
    ITS->SetCoolingFluid(1);   

    // Coarse geometries (warning: no hits are produced with these
    // coarse geometries and they unuseful for reconstruction !):
    //
    //
    // AliITSvPPRcoarseasymm *ITS  = 
    //   new AliITSvPPRcoarseasymm("ITS","New ITS PPR coarse version "
    //                             "with asymmetric services");
    // 1 --> rails in ; 0 --> rails out
    // ITS->SetRails(0);
    // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    // ITS->SetSupportMaterial(0);      
    //
    // AliITS *ITS  = 
    //  new AliITSvPPRcoarsesymm("ITS","New ITS PPR coarse version "
    //                           "with symmetric services");
    // 1 --> rails in ; 0 --> rails out
    // ITS->SetRails(0);                
    // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    // ITS->SetSupportMaterial(0);      
    //
    // Geant3 <-> EUCLID conversion
    // ============================
    //
    // SetEUCLID is a flag to output (=1) or not to output (=0) both
    // geometry and media to two ASCII files (called by default
    // ITSgeometry.euc and ITSgeometry.tme) in a format understandable
    // to the CAD system EUCLID.  The default (=0) means that you dont
    // want to use this facility.
    //
    ITS->SetEUCLID(0);
  }

  if (useTPC) {
    // =================== TPC parameters ============================
    //
    // This allows the user to specify sectors for the SLOW (TPC
    // geometry 2) Simulator. SecAL (SecAU) <0 means that ALL lower
    // (upper) sectors are specified, any value other than that
    // requires at least one sector (lower or upper)to be specified!
    //
    // Reminder: 
    //   sectors 1-24 are lower sectors (1-12 -> z>0, 13-24 -> z<0)
    //   sectors 25-72 are the upper ones (25-48 -> z>0, 49-72 -> z<0)
    //
    //   SecLows - number of lower sectors specified (up to 6)
    //   SecUps  - number of upper sectors specified (up to 12)
    //   Sens    - sensitive strips for the Slow Simulator !!!
    //
    // This does NOT work if all S or L-sectors are specified, i.e.
    // if SecAL or SecAU < 0
    //
    //
    //----------------------------------------------------------------
    //  gROOT->LoadMacro("SetTPCParam.C");
    //  AliTPCParam *param = SetTPCParam();
    AliTPC *TPC = new AliTPCv2("TPC", "Default");
  }

  if (useTOF) {
    // ================== TOF parameters =============================
    AliTOF *TOF = new AliTOFv4T0("TOF", "normal TOF");
  }

  if (useHMPID) {
    // ================== HMPID parameters ============================
    AliHMPID *HMPID = new AliHMPIDv1("HMPID", "normal HMPID");

  }

  if (useZDC) {
    // ================== ZDC parameters =============================
    AliZDC *ZDC = new AliZDCv2("ZDC", "normal ZDC");
  }

  if (useTRD) {
    // ================== TRD parameters =============================
    AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");

    // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 90% Xe + 10% CO2)
    TRD->SetGasMix(1);
    if (geo == kHoles) {
      // With hole in front of PHOS
      TRD->SetPHOShole();
      // With hole in front of HMPID
      TRD->SetHMPIDhole();
    }
    // Switch on TR
    AliTRDsim *TRDsim = TRD->CreateTR();
  }

  if (useFMD) {
    // =================== FMD parameters ============================
    AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
    // FMD->UseDetailed(kFALSE);
    // FMD->UseAssembly();
    // FMD->UseOld();
  }

  if (useMUON) {
    // =================== MUON parameters ===========================
    AliMUON *MUON = new AliMUONv1("MUON", "default");
    // MUON->AddGeometryBuilder(new AliMUONSt1GeometryBuilder(MUON));
    // MUON->AddGeometryBuilder(new AliMUONSt2GeometryBuilder(MUON));
    // MUON->AddGeometryBuilder(new AliMUONSlatGeometryBuilder(MUON));
    // MUON->AddGeometryBuilder(new AliMUONTriggerGeometryBuilder(MUON));
  }

  if (usePHOS) {
    // =================== PHOS parameters ===========================
    AliPHOS *PHOS = new AliPHOSv1("PHOS", "IHEP");
  }

  if (usePMD) {
    // =================== PMD parameters ============================
    AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");
  }

  if (useSTART) {
    // =================== START parameters ==========================
    AliSTART *START = new AliSTARTv1("START", "START Detector");
  }

  if (useEMCAL) {
    // =================== EMCAL parameters ==========================
    AliEMCAL *EMCAL = new AliEMCALv1("EMCAL", "EMCAL_55_25");
  }

  if (useCRT) {
    // =================== CRT parameters ============================
    AliCRT *CRT = new AliCRTv0("CRT", "normal ACORDE");
  }

  if (useVZERO) {
    // =================== V0 parameters =============================
    AliVZERO *VZERO = new AliVZEROv3("VZERO", "normal VZERO");
  }
}

//____________________________________________________________________
Float_t EtaToTheta(Float_t arg)
{
  return (180./TMath::Pi())*2.*TMath::ATan(TMath::Exp(-arg));
}

//____________________________________________________________________
Int_t 
LookupEG(const Char_t* name) 
{
  TString n(name);
  for (Int_t i = 0; i < kEgMax; i++) {
    if (n == egName[i]) return i;
  }
  return -1;
}

//____________________________________________________________________  
AliGenerator* 
GeneratorFactory(EG_t eg, Rad_t rad, TString& comment)  
{
  Int_t isw = 3;
  if (rad == kNoGluonRadiation) isw = 0;
  
  
  AliGenerator * gGener = 0;
  switch (eg) {
  case test50:
    {
      comment = comment.Append(":HIJINGparam test 50 particles");
      AliGenHIJINGpara *gener = new AliGenHIJINGpara(50);
      gener->SetMomentumRange(0, 999999.);
      gener->SetPhiRange(0., 360.);
      // Set pseudorapidity range from -8 to 8.
      Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
      Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
      gener->SetThetaRange(thmin,thmax);
      gGener=gener;
    }
    break;
  case kParam_8000:
    {
      comment = comment.Append(":HIJINGparam N=8000");
      AliGenHIJINGpara *gener = new AliGenHIJINGpara(86030);
      gener->SetMomentumRange(0, 999999.);
      gener->SetPhiRange(0., 360.);
      // Set pseudorapidity range from -8 to 8.
      Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
      Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
      gener->SetThetaRange(thmin,thmax);
      gGener=gener;
    }
    break;
  case kParam_4000:
    {
      comment = comment.Append("HIJINGparam N=4000");
      AliGenHIJINGpara *gener = new AliGenHIJINGpara(43015);
      gener->SetMomentumRange(0, 999999.);
      gener->SetPhiRange(0., 360.);
      // Set pseudorapidity range from -8 to 8.
      Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
      Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
      gener->SetThetaRange(thmin,thmax);
      gGener=gener;
    }
    break;
  case kParam_2000:
    {
      comment = comment.Append("HIJINGparam N=2000");
      AliGenHIJINGpara *gener = new AliGenHIJINGpara(21507);
      gener->SetMomentumRange(0, 999999.);
      gener->SetPhiRange(0., 360.);
      // Set pseudorapidity range from -8 to 8.
      Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
      Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
      gener->SetThetaRange(thmin,thmax);
      gGener=gener;
    }
    break;
  case kParam_fmd:
    {
      comment = comment.Append("HIJINGparam N=100");
      AliGenHIJINGpara *gener = new AliGenHIJINGpara(500);
      gener->SetMomentumRange(0, 999999.);
      gener->SetPhiRange(0., 360.);
      // Set pseudorapidity range from -8 to 8.
      Float_t thmin = EtaToTheta(6);   // theta min. <---> eta max
      Float_t thmax = EtaToTheta(2);  // theta max. <---> eta min 
      gener->SetThetaRange(thmin,thmax);
      gGener=gener;
    }
    break;
    //
    //  Hijing Central
    //
  case kHijing_cent1:
    {
      comment = comment.Append("HIJING cent1");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 5.);
      gGener=gener;
    }
    break;
  case kHijing_cent2:
    {
      comment = comment.Append("HIJING cent2");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 2.);
      gGener=gener;
    }
    break;
    //
    // Hijing Peripheral 
    //
  case kHijing_per1:
    {
      comment = comment.Append("HIJING per1");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(5., 8.6);
      gGener=gener;
    }
    break;
  case kHijing_per2:
    {
      comment = comment.Append("HIJING per2");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(8.6, 11.2);
      gGener=gener;
    }
    break;
  case kHijing_per3:
    {
      comment = comment.Append("HIJING per3");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(11.2, 13.2);
      gGener=gener;
    }
    break;
  case kHijing_per4:
    {
      comment = comment.Append("HIJING per4");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(13.2, 15.);
      gGener=gener;
    }
    break;
  case kHijing_per5:
    {
      comment = comment.Append("HIJING per5");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(15., 100.);
      gGener=gener;
    }
    break;
    //
    //  Jet-Jet
    //
  case kHijing_jj25:
    {
      comment = comment.Append("HIJING Jet 25 GeV");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 5.);
      // trigger
      gener->SetTrigger(1);
      gener->SetPtJet(25.);
      gener->SetRadiation(isw);
      gener->SetSimpleJets(!isw);
      gener->SetJetEtaRange(-0.3,0.3);
      gener->SetJetPhiRange(75., 165.);   
      gGener=gener;
    }
    break;

  case kHijing_jj50:
    {
      comment = comment.Append("HIJING Jet 50 GeV");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 5.);
      // trigger
      gener->SetTrigger(1);
      gener->SetPtJet(50.);
      gener->SetRadiation(isw);
      gener->SetSimpleJets(!isw);
      gener->SetJetEtaRange(-0.3,0.3);
      gener->SetJetPhiRange(75., 165.);   
      gGener=gener;
    }
    break;

  case kHijing_jj75:
    {
      comment = comment.Append("HIJING Jet 75 GeV");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 5.);
      // trigger
      gener->SetTrigger(1);
      gener->SetPtJet(75.);
      gener->SetRadiation(isw);
      gener->SetSimpleJets(!isw);
      gener->SetJetEtaRange(-0.3,0.3);
      gener->SetJetPhiRange(75., 165.);   
      gGener=gener;
    }
    break;

  case kHijing_jj100:
    {
      comment = comment.Append("HIJING Jet 100 GeV");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 5.);
      // trigger
      gener->SetTrigger(1);
      gener->SetPtJet(100.);
      gener->SetRadiation(isw);
      gener->SetSimpleJets(!isw);
      gener->SetJetEtaRange(-0.3,0.3);
      gener->SetJetPhiRange(75., 165.);   
      gGener=gener;
    }
    break;

  case kHijing_jj200:
    {
      comment = comment.Append("HIJING Jet 200 GeV");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 5.);
      // trigger
      gener->SetTrigger(1);
      gener->SetPtJet(200.);
      gener->SetRadiation(isw);
      gener->SetSimpleJets(!isw);
      gener->SetJetEtaRange(-0.3,0.3);
      gener->SetJetPhiRange(75., 165.);   
      gGener=gener;
    }
    break;
    //
    // Gamma-Jet
    //
  case kHijing_gj25:
    {
      comment = comment.Append("HIJING Gamma 25 GeV");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 5.);
      // trigger
      gener->SetTrigger(2);
      gener->SetPtJet(25.);
      gener->SetRadiation(isw);
      gener->SetSimpleJets(!isw);
      gener->SetJetEtaRange(-0.12, 0.12);
      gener->SetJetPhiRange(220., 320.);
      gGener=gener;
    }
    break;

  case kHijing_gj50:
    {
      comment = comment.Append("HIJING Gamma 50 GeV");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 5.);
      // trigger
      gener->SetTrigger(2);
      gener->SetPtJet(50.);
      gener->SetRadiation(isw);
      gener->SetSimpleJets(!isw);
      gener->SetJetEtaRange(-0.12, 0.12);
      gener->SetJetPhiRange(220., 320.);
      gGener=gener;
    }
    break;

  case kHijing_gj75:
    {
      comment = comment.Append("HIJING Gamma 75 GeV");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 5.);
      // trigger
      gener->SetTrigger(2);
      gener->SetPtJet(75.);
      gener->SetRadiation(isw);
      gener->SetSimpleJets(!isw);
      gener->SetJetEtaRange(-0.12, 0.12);
      gener->SetJetPhiRange(220., 320.);
      gGener=gener;
    }
    break;

  case kHijing_gj100:
    {
      comment = comment.Append("HIJING Gamma 100 GeV");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 5.);
      // trigger
      gener->SetTrigger(2);
      gener->SetPtJet(100.);
      gener->SetRadiation(isw);
      gener->SetSimpleJets(!isw);
      gener->SetJetEtaRange(-0.12, 0.12);
      gener->SetJetPhiRange(220., 320.);
      gGener=gener;
    }
    break;

  case kHijing_gj200:
    {
      comment = comment.Append("HIJING Gamma 200 GeV");
      AliGenHijing *gener = HijingStandard();
      // impact parameter range
      gener->SetImpactParameterRange(0., 5.);
      // trigger
      gener->SetTrigger(2);
      gener->SetPtJet(200.);
      gener->SetRadiation(isw);
      gener->SetSimpleJets(!isw);
      gener->SetJetEtaRange(-0.12, 0.12);
      gener->SetJetPhiRange(220., 320.);
      gGener=gener;
    }
    break;
  case kHijing_pA:
    {
      comment = comment.Append("HIJING pA");

      AliGenCocktail *gener  = new AliGenCocktail();

      AliGenHijing   *hijing = new AliGenHijing(-1);
      // centre of mass energy 
      hijing->SetEnergyCMS(TMath::Sqrt(82./208.) * 14000.);
      // impact parameter range
      hijing->SetImpactParameterRange(0., 15.);
      // reference frame
      hijing->SetReferenceFrame("CMS");
      hijing->SetBoostLHC(1);
      // projectile
      hijing->SetProjectile("P", 1, 1);
      hijing->SetTarget    ("A", 208, 82);
      // tell hijing to keep the full parent child chain
      hijing->KeepFullEvent();
      // enable jet quenching
      hijing->SetJetQuenching(0);
      // enable shadowing
      hijing->SetShadowing(1);
      // Don't track spectators
      hijing->SetSpectators(0);
      // kinematic selection
      hijing->SetSelectAll(0);
      //
      AliGenSlowNucleons*  gray    = new AliGenSlowNucleons(1);
      AliSlowNucleonModel* model   = new AliSlowNucleonModelExp();
      gray->SetSlowNucleonModel(model);
      gray->SetDebug(1);
      gener->AddGenerator(hijing,"Hijing pPb", 1);
      gener->AddGenerator(gray,  "Gray Particles",1);
      gGener=gener;
    }
    break;
  case kPythia6:
    {
      comment = comment.Append(":Pythia p-p @ 14 TeV");
      AliGenPythia *gener = new AliGenPythia(-1); 
      gener->SetMomentumRange(0,999999);
      gener->SetThetaRange(0., 180.);
      gener->SetYRange(-12,12);
      gener->SetPtRange(0,1000);
      gener->SetProcess(kPyMb);
      gener->SetEnergyCMS(14000.);
      gGener=gener;
    }
    break;
  case kPythia6Jets20_24:
    {
      comment = comment.Append(":Pythia jets 20-24 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(20., 24.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kPythia6Jets24_29:
    {
      comment = comment.Append(":Pythia jets 24-29 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(24., 29.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kPythia6Jets29_35:
    {
      comment = comment.Append(":Pythia jets 29-35 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(29., 35.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kPythia6Jets35_42:
    {
      comment = comment.Append(":Pythia jets 35-42 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(35., 42.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kPythia6Jets42_50:
    {
      comment = comment.Append(":Pythia jets 42-50 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(42., 50.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kPythia6Jets50_60:
    {
      comment = comment.Append(":Pythia jets 50-60 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(50., 60.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kPythia6Jets60_72:
    {
      comment = comment.Append(":Pythia jets 60-72 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(60., 72.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kPythia6Jets72_86:
    {
      comment = comment.Append(":Pythia jets 72-86 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(72., 86.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kPythia6Jets86_104:
    {
      comment = comment.Append(":Pythia jets 86-104 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(86., 104.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kPythia6Jets104_125:
    {
      comment = comment.Append(":Pythia jets 105-125 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(104., 125.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kPythia6Jets125_150:
    {
      comment = comment.Append(":Pythia jets 125-150 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(125., 150.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kPythia6Jets150_180:
    {
      comment = comment.Append(":Pythia jets 150-180 GeV @ 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);//        Centre of mass energy
      gener->SetProcess(kPyJets);//        Process type
      gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
      gener->SetJetPhiRange(0., 360.);
      gener->SetJetEtRange(10., 1000.);
      gener->SetGluonRadiation(1,1);
      //    gener->SetPtKick(0.);
      //   Structure function
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(150., 180.);// Pt transfer of the hard scattering
      gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
      gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
      gGener=gener;
    }
    break;
  case kD0PbPb5500:
    {
      comment = comment.Append(" D0 in Pb-Pb at 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(10);
      gener->SetProcess(kPyD0PbPbMNR);
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(2.1,-1.0);
      gener->SetEnergyCMS(5500.);
      gener->SetNuclei(208,208);
      gener->SetForceDecay(kHadronicD);
      gener->SetYRange(-2,2);
      gener->SetFeedDownHigherFamily(kFALSE);
      gener->SetStackFillOpt(AliGenPythia::kParentSelection);
      gener->SetCountMode(AliGenPythia::kCountParents);
      gGener=gener;
    }
    break;
  case kCharmSemiElPbPb5500:
    {
      comment = comment.Append(" Charm in Pb-Pb at 5.5 TeV");
      AliGenPythia * gener = new AliGenPythia(10);
      gener->SetProcess(kPyCharmPbPbMNR);
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(2.1,-1.0);
      gener->SetEnergyCMS(5500.);
      gener->SetNuclei(208,208);
      gener->SetForceDecay(kSemiElectronic);
      gener->SetYRange(-2,2);
      gener->SetFeedDownHigherFamily(kFALSE);
      gener->SetCountMode(AliGenPythia::kCountParents);
      gGener=gener;
    }
    break;
  case kBeautySemiElPbPb5500:
    {
      comment = comment.Append(" Beauty in Pb-Pb at 5.5 TeV");
      AliGenPythia *gener = new AliGenPythia(10);
      gener->SetProcess(kPyBeautyPbPbMNR);
      gener->SetStrucFunc(kCTEQ4L);
      gener->SetPtHard(2.75,-1.0);
      gener->SetEnergyCMS(5500.);
      gener->SetNuclei(208,208);
      gener->SetForceDecay(kSemiElectronic);
      gener->SetYRange(-2,2);
      gener->SetFeedDownHigherFamily(kFALSE);
      gener->SetCountMode(AliGenPythia::kCountParents);
      gGener=gener;
    }
    break;
  case kCocktailTRD:
    {
      comment = comment.Append(" Cocktail for TRD at 5.5 TeV");
      AliGenCocktail *gener  = new AliGenCocktail();

      AliGenParam *jpsi = new AliGenParam(10,
					  new AliGenMUONlib(),
					  AliGenMUONlib::kJpsiFamily,
					  "Vogt PbPb");

      jpsi->SetPtRange(0, 100);
      jpsi->SetYRange(-1., +1.);
      jpsi->SetForceDecay(kDiElectron);

      AliGenParam *ups = new AliGenParam(10,
					 new AliGenMUONlib(),
					 AliGenMUONlib::kUpsilonFamily,
					 "Vogt PbPb");
      ups->SetPtRange(0, 100);
      ups->SetYRange(-1., +1.);
      ups->SetForceDecay(kDiElectron);
	
      AliGenParam *charm = new AliGenParam(10,
					   new AliGenMUONlib(), 
					   AliGenMUONlib::kCharm,
					   "central");
      charm->SetPtRange(0, 100);
      charm->SetYRange(-1.5, +1.5);
      charm->SetForceDecay(kSemiElectronic);
	
	
      AliGenParam *beauty = new AliGenParam(10,
					    new AliGenMUONlib(), 
					    AliGenMUONlib::kBeauty,
					    "central");
      beauty->SetPtRange(0, 100);
      beauty->SetYRange(-1.5, +1.5);
      beauty->SetForceDecay(kSemiElectronic);

      gener->AddGenerator(jpsi,"J/psi",1);
      gener->AddGenerator(ups,"Upsilon",1);
      gener->AddGenerator(charm,"Charm",1);
      gener->AddGenerator(beauty,"Beauty",1);
      gGener=gener;
    }
    break;
  case kPyJJ:
    {
      comment = comment.Append(" Jet-jet at 5.5 TeV");
      AliGenPythia *gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);
      gener->SetProcess(kPyJets);
      Double_t ptHardMin=10.0, ptHardMax=-1.0;
      gener->SetPtHard(ptHardMin,ptHardMax);
      gener->SetYHard(-0.7,0.7);
      gener->SetJetEtaRange(-0.2,0.2);
      gener->SetEventListRange(0,1);
      gGener=gener;
    }
    break;
  case kPyGJ:
    {
      comment = comment.Append(" Gamma-jet at 5.5 TeV");
      AliGenPythia *gener = new AliGenPythia(-1);
      gener->SetEnergyCMS(5500.);
      gener->SetProcess(kPyDirectGamma);
      Double_t ptHardMin=10.0, ptHardMax=-1.0;
      gener->SetPtHard(ptHardMin,ptHardMax);
      gener->SetYHard(-1.0,1.0);
      gener->SetGammaEtaRange(-0.13,0.13);
      gener->SetGammaPhiRange(210.,330.);
      gener->SetEventListRange(0,1);
      gGener=gener;
    }
    break;
  case kMuonCocktailCent1:
    {
      comment = comment.Append(" Muon Cocktail Cent1");
      AliGenMUONCocktail * gener = new AliGenMUONCocktail();
      gener->SetPtRange(1.0,100.);       // Transverse momentum range   
      gener->SetPhiRange(0.,360.);    // Azimuthal angle range  
      gener->SetYRange(-4.0,-2.4);
      gener->SetMuonPtCut(0.8);
      gener->SetMuonThetaCut(171.,178.);
      gener->SetMuonMultiplicity(2);
      gener->SetNumberOfCollisions(1626.);  //Centrality class Cent1 for PDC04
      gener->SetNumberOfParticipants(359.4);//Centrality class Cent1 for PDC04
      gGener=gener;
    }
    break;
  case kMuonCocktailPer1:
    {
      comment = comment.Append(" Muon Cocktail Per1");
      AliGenMUONCocktail * gener = new AliGenMUONCocktail();
      gener->SetPtRange(1.0,100.);       // Transverse momentum range   
      gener->SetPhiRange(0.,360.);    // Azimuthal angle range  
      gener->SetYRange(-4.0,-2.4);
      gener->SetMuonPtCut(0.8);
      gener->SetMuonThetaCut(171.,178.);
      gener->SetMuonMultiplicity(2);
      gener->SetNumberOfCollisions(820.0);//Centrality class Per1 for PDC04
      gener->SetNumberOfParticipants(229.3);//Centrality class Per1 for PDC04
      gGener=gener;
    }
    break;
  case kMuonCocktailPer4:
    {
      comment = comment.Append(" Muon Cocktail Per4");
      AliGenMUONCocktail * gener = new AliGenMUONCocktail();
      gener->SetPtRange(1.0,100.);       // Transverse momentum range   
      gener->SetPhiRange(0.,360.);    // Azimuthal angle range  
      gener->SetYRange(-4.0,-2.4);
      gener->SetMuonPtCut(0.8);
      gener->SetMuonThetaCut(171.,178.);
      gener->SetMuonMultiplicity(2);
      gener->SetNumberOfCollisions(13.6);//Centrality class Per4 for PDC04
      gener->SetNumberOfParticipants(13.3);//Centrality class Per4 for PDC04
      gGener=gener;
    }
    break;
  case kMuonCocktailCent1HighPt:
    {
      comment = comment.Append(" Muon Cocktail HighPt Cent1");
      AliGenMUONCocktail * gener = new AliGenMUONCocktail();
      gener->SetPtRange(1.0,100.);       // Transverse momentum range   
      gener->SetPhiRange(0.,360.);    // Azimuthal angle range  
      gener->SetYRange(-4.0,-2.4);
      gener->SetMuonPtCut(2.5);
      gener->SetMuonThetaCut(171.,178.);
      gener->SetMuonMultiplicity(2);
      gener->SetNumberOfCollisions(1626.);  //Centrality class Cent1 for PDC04
      gener->SetNumberOfParticipants(359.4);//Centrality class Cent1 for PDC04
      gGener=gener;
    }
    break;
  case kMuonCocktailPer1HighPt :
    {
      comment = comment.Append(" Muon Cocktail HighPt Per1");
      AliGenMUONCocktail * gener = new AliGenMUONCocktail();
      gener->SetPtRange(1.0,100.);       // Transverse momentum range   
      gener->SetPhiRange(0.,360.);    // Azimuthal angle range  
      gener->SetYRange(-4.0,-2.4);
      gener->SetMuonPtCut(2.5);
      gener->SetMuonThetaCut(171.,178.);
      gener->SetMuonMultiplicity(2);
      gener->SetNumberOfCollisions(820.0);//Centrality class Per1 for PDC04
      gener->SetNumberOfParticipants(229.3);//Centrality class Per1 for PDC04
      gGener=gener;
    }
    break;
  case kMuonCocktailPer4HighPt:
    {
      comment = comment.Append(" Muon Cocktail HighPt Per4");
      AliGenMUONCocktail * gener = new AliGenMUONCocktail();
      gener->SetPtRange(1.0,100.);       // Transverse momentum range   
      gener->SetPhiRange(0.,360.);    // Azimuthal angle range  
      gener->SetYRange(-4.0,-2.4);
      gener->SetMuonPtCut(2.5);
      gener->SetMuonThetaCut(171.,178.);
      gener->SetMuonMultiplicity(2);
      gener->SetNumberOfCollisions(13.6);//Centrality class Per4 for PDC04
      gener->SetNumberOfParticipants(13.3);//Centrality class Per4 for PDC04
      gGener=gener;
    }
    break;
  case kMuonCocktailCent1Single:
    {
      comment = comment.Append(" Muon Cocktail Single Cent1");
      AliGenMUONCocktail * gener = new AliGenMUONCocktail();
      gener->SetPtRange(1.0,100.);       // Transverse momentum range   
      gener->SetPhiRange(0.,360.);    // Azimuthal angle range  
      gener->SetYRange(-4.0,-2.4);
      gener->SetMuonPtCut(0.8);
      gener->SetMuonThetaCut(171.,178.);
      gener->SetMuonMultiplicity(1);
      gener->SetNumberOfCollisions(1626.);  //Centrality class Cent1 for PDC04
      gener->SetNumberOfParticipants(359.4);//Centrality class Cent1 for PDC04
      gGener=gener;
    }
    break;
  case kMuonCocktailPer1Single :
    {
      comment = comment.Append(" Muon Cocktail Single Per1");
      AliGenMUONCocktail * gener = new AliGenMUONCocktail();
      gener->SetPtRange(1.0,100.);       // Transverse momentum range   
      gener->SetPhiRange(0.,360.);    // Azimuthal angle range  
      gener->SetYRange(-4.0,-2.4);
      gener->SetMuonPtCut(0.8);
      gener->SetMuonThetaCut(171.,178.);
      gener->SetMuonMultiplicity(1);
      gener->SetNumberOfCollisions(820.0);//Centrality class Per1 for PDC04
      gener->SetNumberOfParticipants(229.3);//Centrality class Per1 for PDC04
      gGener=gener;
    }
    break;
  case kMuonCocktailPer4Single:
    {
      comment = comment.Append(" Muon Cocktail Single Per4");
      AliGenMUONCocktail * gener = new AliGenMUONCocktail();
      gener->SetPtRange(1.0,100.);       // Transverse momentum range   
      gener->SetPhiRange(0.,360.);    // Azimuthal angle range  
      gener->SetYRange(-4.0,-2.4);
      gener->SetMuonPtCut(0.8);
      gener->SetMuonThetaCut(171.,178.);
      gener->SetMuonMultiplicity(1);
      gener->SetNumberOfCollisions(13.6);//Centrality class Per4 for PDC04
      gener->SetNumberOfParticipants(13.3);//Centrality class Per4 for PDC04
      gGener=gener;
    }
    break;
  case kFMD1Flat: 
    {
      comment = comment.Append(" Flat in FMD1 range");
      AliGenBox* gener = new AliGenBox(2000);
      gener->SetPart(211);
      gener->SetMomentumRange(3,4);
      gener->SetPhiRange(0, 360);
      gener->SetThetaRange(0.77, 3.08);
      gGener = gener;
    }
    break;
  case kFMD2Flat: 
    {
      comment = comment.Append(" Flat in FMD2 range");
      AliGenBox* gener = new AliGenBox(2000);
      gener->SetPart(211);
      gener->SetMomentumRange(3,4);
      gener->SetPhiRange(0, 360);
      gener->SetThetaRange(2.95, 20.42);
      gGener = gener;
    }
    break;
  case kFMD3Flat: 
    {
      comment = comment.Append(" Flat in FMD3 range");
      AliGenBox* gener = new AliGenBox(2000);
      gener->SetPart(211);
      gener->SetMomentumRange(3,4);
      gener->SetPhiRange(0, 360);
      gener->SetThetaRange(155.97, 176.73);
      gGener = gener;
    }
    break;
  case kFMDFlat:
    {
      comment = comment.Append(" Flat in FMD range");
      AliGenCocktail* gener = new AliGenCocktail();
      gener->SetMomentumRange(3,4);
      gener->SetPhiRange(0, 360);
      AliGenBox* gener3 = new AliGenBox(2000);
      gener3->SetThetaRange(155.97, 176.73);
      gener3->SetPart(211);
      gener->AddGenerator(gener3, "FMD3", .33);
      AliGenBox* gener2 = new AliGenBox(2000);
      gener2->SetThetaRange(2.95, 20.42);
      gener2->SetPart(211);
      gener->AddGenerator(gener2, "FMD2", .33);
      AliGenBox* gener1 = new AliGenBox(2000);
      gener1->SetThetaRange(0.77, 3.08);
      gener1->SetPart(211);
      gener->AddGenerator(gener1, "FMD1", .34);
      gGener = gener;
    }
    break;
    
  default: break;
  }
  return gGener;
}

//____________________________________________________________________
AliGenHijing* 
HijingStandard()
{
  AliGenHijing *gener = new AliGenHijing(-1);
  // centre of mass energy 
  gener->SetEnergyCMS(5500.);
  // reference frame
  gener->SetReferenceFrame("CMS");
  // projectile
  gener->SetProjectile("A", 208, 82);
  gener->SetTarget    ("A", 208, 82);
  // tell hijing to keep the full parent child chain
  gener->KeepFullEvent();
  // enable jet quenching
  gener->SetJetQuenching(1);
  // enable shadowing
  gener->SetShadowing(1);
  // neutral pion and heavy particle decays switched off
  gener->SetDecaysOff(1);
  // Don't track spectators
  gener->SetSpectators(0);
  // kinematic selection
  gener->SetSelectAll(0);
  return gener;
}


//____________________________________________________________________
void 
ProcessEnvironmentVars(EG_t& eg, Int_t& seed)
{
  // Run type
  if (gSystem->Getenv("CONFIG_RUN_TYPE")) {
    Int_t eg1 = LookupEG(gSystem->Getenv("CONFIG_RUN_TYPE"));
    if  (eg1 >= 0) eg = EG_t(eg1);
  }
  // Random Number seed
  if (gSystem->Getenv("CONFIG_SEED")) {
    seed = atoi(gSystem->Getenv("CONFIG_SEED"));
  }
}

//____________________________________________________________________
//
// EOF
// 
