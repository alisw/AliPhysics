// Macro to generate e, pi, mu, K, p 200 each with box generator
// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"ConfigPPR.C++")
// Prashant Shukla (shukla@physi.uni-heidelberg.de)

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include <TPDGCode.h>
#include <TF1.h>
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "STEER/AliGenerator.h"
#include "STEER/AliLog.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "EVGEN/AliGenHIJINGpara.h"
#include "THijing/AliGenHijing.h"
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenSlowNucleons.h"
#include "EVGEN/AliSlowNucleonModelExp.h"
#include "EVGEN/AliGenParam.h"
#include "EVGEN/AliGenMUONlib.h"
#include "EVGEN/AliGenSTRANGElib.h"
#include "EVGEN/AliGenMUONCocktail.h"
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenGeVSim.h"
#include "EVGEN/AliGeVSimParticle.h"
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
#include "RICH/AliRICHv1.h"
#include "ZDC/AliZDCv2.h"
#include "TRD/AliTRDv1.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "MUON/AliMUONSt1GeometryBuilderV2.h"
#include "MUON/AliMUONSt2GeometryBuilderV2.h"
#include "MUON/AliMUONSlatGeometryBuilder.h"
#include "MUON/AliMUONTriggerGeometryBuilder.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "START/AliSTARTv1.h"
#include "EMCAL/AliEMCALv1.h"
#include "CRT/AliCRTv0.h"
#include "VZERO/AliVZEROv5.h"
#endif

enum PprRun_t 
{
    test50,
    kParam_8000,   kParam_4000,  kParam_2000, 
    kHijing_cent1, kHijing_cent2, 
    kHijing_per1,  kHijing_per2, kHijing_per3, kHijing_per4,  kHijing_per5,
    kCocktailTRD,  kpieTRD
};

const char* pprRunName[] = {
    "test50",
    "kParam_8000", "kParam_4000",  "kParam_2000", 
    "kHijing_cent1", "kHijing_cent2", 
    "kHijing_per1", "kHijing_per2", "kHijing_per3", "kHijing_per4",  
    "kHijing_per5",
    "kCocktailTRD", "kpieTRD"
};

enum PprGeo_t 
{
    kHoles, kNoHoles
};

enum PprRad_t
{
    kGluonRadiation, kNoGluonRadiation
};

enum PprMag_t
{
    k2kG, k4kG, k5kG
};


// This part for configuration    
//static PprRun_t srun = test50;
static PprRun_t srun = kpieTRD;
static PprGeo_t sgeo = kNoHoles;
static PprRad_t srad = kGluonRadiation;
static PprMag_t smag = k5kG;
static Int_t    sseed = 0; //Set 0 to use the current time

// Comment line 
static TString  comment;

// Functions
Float_t EtaToTheta(Float_t arg);
AliGenerator* GeneratorFactory(PprRun_t srun);
AliGenHijing* HijingStandard();
void ProcessEnvironmentVars();

void Config()
{
    // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
    // Theta range given through pseudorapidity limits 22/6/2001

    // Get settings from environment variables
    ProcessEnvironmentVars();

    // Set Random Number seed
    gRandom->SetSeed(sseed);
    cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<endl; 


   // libraries required by geant321
#if defined(__CINT__)
    gSystem->Load("libgeant321");
#endif

    new     TGeant3TGeo("C++ Interface to Geant3");

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

    //
    // Set External decayer
    AliDecayer *decayer = new AliDecayerPythia();

    decayer->SetForceDecay(kAll);

    decayer->Init();
    gMC->SetExternalDecayer(decayer);
    //
    //
    //=======================================================================
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

    // Debug and log level
    //    AliLog::SetGlobalDebugLevel(0);
    //    AliLog::SetGlobalLogLevel(AliLog::kError);

    // Generator Configuration
    AliGenerator* gener = GeneratorFactory(srun);
    gener->SetOrigin(0, 0, 0);    // vertex position
    gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position
    gener->SetCutVertexZ(1.);     // Truncate at 1 sigma
    gener->SetVertexSmear(kPerEvent); 
    gener->SetTrackingFlag(1);
    gener->Init();
    
    if (smag == k2kG) {
	comment = comment.Append(" | L3 field 0.2 T");
    } else if (smag == k4kG) {
	comment = comment.Append(" | L3 field 0.4 T");
    } else if (smag == k5kG) {
	comment = comment.Append(" | L3 field 0.5 T");
    }
    
    
    if (srad == kGluonRadiation)
    {
	comment = comment.Append(" | Gluon Radiation On");
	
    } else {
	comment = comment.Append(" | Gluon Radiation Off");
    }

    if (sgeo == kHoles)
    {
	comment = comment.Append(" | Holes for PHOS/RICH");
	
    } else {
	comment = comment.Append(" | No holes for PHOS/RICH");
    }

    printf("\n \n Comment: %s \n \n", comment.Data());
    
    
// Field (L3 0.4 T)
    AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., smag);
    //    AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 0, 10., smag);  // B = 0
    field->SetL3ConstField(0); //Using const. field in the barrel
    rl->CdGAFile();
    gAlice->SetField(field);    
//
    Int_t   iABSO   = 1;
    Int_t   iDIPO   = 1;
    Int_t   iFMD    = 1;
    Int_t   iFRAME  = 1;
    Int_t   iHALL   = 1;
    Int_t   iITS    = 1;
    Int_t   iMAG    = 1;
    Int_t   iMUON   = 0;
    Int_t   iPHOS   = 1;
    Int_t   iPIPE   = 1;
    Int_t   iPMD    = 0;
    Int_t   iRICH   = 0;
    Int_t   iSHIL   = 1;
    Int_t   iSTART  = 1;
    Int_t   iTOF    = 0;
    Int_t   iTPC    = 1;
    Int_t   iTRD    = 1;
    Int_t   iZDC    = 1;
    Int_t   iEMCAL  = 0;
    Int_t   iVZERO  = 1;
    Int_t   iCRT    = 0;

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
	if (sgeo == kHoles) {
	    FRAME->SetHoles(1);
	} else {
	    FRAME->SetHoles(0);
	}
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
	ITS->SetReadDet(kTRUE);	  // don't touch this parameter if you're not an ITS developer
    //    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
	ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
	ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
	ITS->SetThicknessChip1(200.);  // chip thickness on layer 1 must be in the range [150,300]
	ITS->SetThicknessChip2(200.);  // chip thickness on layer 2 must be in the range [150,300]
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
        //============================ TPC parameters ================================
        // --- This allows the user to specify sectors for the SLOW (TPC geometry 2)
        // --- Simulator. SecAL (SecAU) <0 means that ALL lower (upper)
        // --- sectors are specified, any value other than that requires at least one 
        // --- sector (lower or upper)to be specified!
        // --- Reminder: sectors 1-24 are lower sectors (1-12 -> z>0, 13-24 -> z<0)
        // ---           sectors 25-72 are the upper ones (25-48 -> z>0, 49-72 -> z<0)
        // --- SecLows - number of lower sectors specified (up to 6)
        // --- SecUps - number of upper sectors specified (up to 12)
        // --- Sens - sensitive strips for the Slow Simulator !!!
        // --- This does NOT work if all S or L-sectors are specified, i.e.
        // --- if SecAL or SecAU < 0
        //
        //
        //-----------------------------------------------------------------------------

        //  gROOT->LoadMacro("SetTPCParam.C");
        //  AliTPCParam *param = SetTPCParam();
        AliTPC *TPC = new AliTPCv2("TPC", "Default");

        // All sectors included 
        TPC->SetSecAL(-1);
        TPC->SetSecAU(-1);

    }


    if (iTOF) {
        //=================== TOF parameters ============================
	AliTOF *TOF = new AliTOFv4T0("TOF", "normal TOF");
    }


    if (iRICH)
    {
        //=================== RICH parameters ===========================
        AliRICH *RICH = new AliRICHv1("RICH", "normal RICH");

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

        // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 90% Xe + 10% CO2)
        TRD->SetGasMix(1);
	if (sgeo == kHoles) {
	    // With hole in front of PHOS
	    TRD->SetPHOShole();
	    // With hole in front of RICH
	    TRD->SetRICHhole();
	}
	    // Switch on TR
	    AliTRDsim *TRDsim = TRD->CreateTR();
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
	((AliMUONv1*)MUON)->SetStepManagerVersionDE(true);
	MUON->AddGeometryBuilder(new AliMUONSt1GeometryBuilderV2(MUON));
	MUON->AddGeometryBuilder(new AliMUONSt2GeometryBuilderV2(MUON));
	MUON->AddGeometryBuilder(new AliMUONSlatGeometryBuilder(MUON));
	MUON->AddGeometryBuilder(new AliMUONTriggerGeometryBuilder(MUON));
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

    if (iSTART)
    {
        //=================== START parameters ============================
        AliSTART *START = new AliSTARTv1("START", "START Detector");
    }

    if (iEMCAL)
    {
        //=================== EMCAL parameters ============================
        AliEMCAL *EMCAL = new AliEMCALv1("EMCAL", "EMCAL_55_25");
    }

     if (iCRT)
    {
        //=================== CRT parameters ============================
        AliCRT *CRT = new AliCRTv0("CRT", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== CRT parameters ============================
        AliVZERO *VZERO = new AliVZEROv5("VZERO", "normal VZERO");
    }
 
             
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}



AliGenerator* GeneratorFactory(PprRun_t srun) {
    Int_t isw = 3;
    if (srad == kNoGluonRadiation) isw = 0;
    

    AliGenerator * gGener = 0x0;
    switch (srun) {
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
    case kCocktailTRD:
      {
	comment = comment.Append(" Cocktail for TRD at 5.5 TeV");
	AliGenCocktail *gener  = new AliGenCocktail();
	
	AliGenParam *phi = new AliGenParam(10,
                                           new AliGenMUONlib(),
                                           AliGenMUONlib::kPhi,
                                           "Vogt PbPb");

	phi->SetPtRange(0, 100);
	phi->SetYRange(-1., +1.);
	phi->SetForceDecay(kDiElectron);

	AliGenParam *omega = new AliGenParam(10,
					     new AliGenMUONlib(),
					     AliGenMUONlib::kOmega,
					     "Vogt PbPb");

	omega->SetPtRange(0, 100);
	omega->SetYRange(-1., +1.);
	omega->SetForceDecay(kDiElectron);
	
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

	AliGenParam *beautyJ = new AliGenParam(10,
					       new AliGenMUONlib(), 
					       AliGenMUONlib::kBeauty,
					       "central");
	beautyJ->SetPtRange(0, 100);
	beautyJ->SetYRange(-1.5, +1.5);
	beautyJ->SetForceDecay(kBJpsiDiElectron);

	gener->AddGenerator(phi,"Phi",1);
	gener->AddGenerator(omega,"Omega",1);
	gener->AddGenerator(jpsi,"J/psi",1);
	gener->AddGenerator(ups,"Upsilon",1);
	gener->AddGenerator(charm,"Charm",1);
	gener->AddGenerator(beauty,"Beauty",1);
	gener->AddGenerator(beautyJ,"J/Psi from Beauty",1);
	gGener=gener;
      }
      break;
    case kpieTRD:
      {
	comment = comment.Append("e, pi for TRD");
	AliGenCocktail *gener  = new AliGenCocktail();
	Double_t momen=2.0;
	Int_t Npart=200;
	AliGenBox *electron = new AliGenBox(Npart);
	electron->SetMomentumRange(momen,momen);
	electron->SetPhiRange(0,360);
	Float_t thmin = EtaToTheta(.9);   // theta min. <---> eta max
	Float_t thmax = EtaToTheta(-.9);  // theta max. <---> eta min 
	electron->SetThetaRange(thmin,thmax);
	electron->SetOrigin(0,0,0);      //vertex position
	electron->SetSigma(0,0,0);       //Sigma in (X,Y,Z) (cm) on IP position
	electron->SetPart(11);              //GEANT particle type

	AliGenBox *pion = new AliGenBox(Npart);
	pion->SetMomentumRange(momen,momen);
	pion->SetPhiRange(0,360);
	pion->SetThetaRange(thmin,thmax);
	pion->SetOrigin(0,0,0);   	//vertex position
	pion->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
	pion->SetPart(211);              //GEANT particle type

	AliGenBox *muon = new AliGenBox(Npart);
	muon->SetMomentumRange(momen,momen);
	muon->SetPhiRange(0,360);
	muon->SetThetaRange(thmin,thmax);
	muon->SetOrigin(0,0,0);   	//vertex position
	muon->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
	muon->SetPart(13);              //GEANT particle type

	AliGenBox *kaon = new AliGenBox(Npart);
	kaon->SetMomentumRange(momen,momen);
	kaon->SetPhiRange(0,360);
	kaon->SetThetaRange(thmin,thmax);
	kaon->SetOrigin(0,0,0);   	//vertex position
	kaon->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
	kaon->SetPart(321);              //GEANT particle type

	AliGenBox *proton = new AliGenBox(Npart);
	proton->SetMomentumRange(momen,momen);
	proton->SetPhiRange(0,360);
	proton->SetThetaRange(thmin,thmax);
	proton->SetOrigin(0,0,0);   	//vertex position
	proton->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
	proton->SetPart(2212);              //GEANT particle type

	gener->AddGenerator(electron,"electron",1);
	gener->AddGenerator(pion,"pion",1);
	gener->AddGenerator(muon,"muon",1);
	gener->AddGenerator(kaon,"kaon",1);
	gener->AddGenerator(proton,"proton",1);
	gGener=gener;
      }
      break;
    default: break;
    }
    return gGener;
}

void ProcessEnvironmentVars()
{
    // Run type
    if (gSystem->Getenv("CONFIG_RUN_TYPE")) {
      for (Int_t iRun = 0; iRun < kRunMax; iRun++) {
	if (strcmp(gSystem->Getenv("CONFIG_RUN_TYPE"), pprRunName[iRun])==0) {
	  srun = (PprRun_t)iRun;
	  cout<<"Run type set to "<<pprRunName[iRun]<<endl;
	}
      }
    }

    // Random Number seed
    if (gSystem->Getenv("CONFIG_SEED")) {
      sseed = atoi(gSystem->Getenv("CONFIG_SEED"));
    }
}
